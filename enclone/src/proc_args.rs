// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

use crate::proc_args2::*;
use crate::proc_args3::*;
use crate::proc_args_check::*;
use enclone_core::defs::*;
use enclone_core::testlist::*;
use io_utils::*;
use perf_stats::*;
use regex::Regex;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::{env, time::Instant};
use string_utils::*;
use vector_utils::*;

// Process arguments.

pub fn proc_args(mut ctl: &mut EncloneControl, args: &Vec<String>) {
    // Knobs.

    let heur = ClonotypeHeuristics {
        max_diffs: 50,
        ref_v_trim: 15,
        ref_j_trim: 15,
    };
    ctl.heur = heur;

    // Form the combined set of command-line arguments and "command-line" arguments
    // implied by environment variables.

    let targs = Instant::now();
    let mut args = args.clone();
    let mut args2 = Vec::<String>::new();
    args2.push(args[0].clone());
    for (key, value) in env::vars() {
        if key.starts_with("ENCLONE_") {
            args2.push(format!("{}={}", key.after("ENCLONE_"), value));
        }
    }
    for i in 1..args.len() {
        args2.push(args[i].clone());
    }
    args = args2;

    // Test for internal run.

    for (key, value) in env::vars() {
        if (key == "HOST" || key == "HOSTNAME") && value.ends_with(".fuzzplex.com") {
            ctl.gen_opt.internal_run = true;
        }
    }
    for i in 1..args.len() {
        if args[i] == "FORCE_EXTERNAL".to_string() {
            ctl.gen_opt.internal_run = false;
        }
    }
    if ctl.gen_opt.internal_run {
        ctl.gen_opt.current_ref = true; // not sure this is right
        ctl.gen_opt.pre = vec![
            format!("/mnt/assembly/vdj/current{}", TEST_FILES_VERSION),
            format!("enclone/test/inputs"),
        ];
    }

    // Set up general options.

    ctl.gen_opt.h5_pre = true;
    ctl.gen_opt.min_cells_exact = 1;
    ctl.gen_opt.min_chains_exact = 1;
    ctl.gen_opt.exact = None;
    for i in 1..args.len() {
        if args[i].starts_with("PRE=") {
            let pre = args[i].after("PRE=").split(',').collect::<Vec<&str>>();
            ctl.gen_opt.pre.clear();
            for x in pre.iter() {
                ctl.gen_opt.pre.push(x.to_string());
            }
        }
    }
    ctl.gen_opt.full_counts = true;
    ctl.silent = true;

    // Set up clonotyping control parameters.

    ctl.clono_filt_opt.ncells_low = 1;
    ctl.clono_filt_opt.ncells_high = 1_000_000_000;
    ctl.clono_filt_opt.min_umi = 0;
    ctl.clono_filt_opt.max_chains = 1000000;
    ctl.clono_filt_opt.qual_filter = true;
    ctl.clono_filt_opt.weak_chains = true;
    ctl.clono_filt_opt.weak_onesies = true;
    ctl.clono_filt_opt.weak_foursies = true;
    ctl.clono_filt_opt.bc_dup = true;
    ctl.clono_filt_opt.max_datasets = 1000000000;
    ctl.clono_filt_opt.umi_filt = true;
    ctl.clono_filt_opt.umi_ratio_filt = true;

    ctl.clono_print_opt.amino = vec![
        "cdr3".to_string(),
        "var".to_string(),
        "share".to_string(),
        "donor".to_string(),
    ];
    ctl.clono_print_opt.cvars = vec!["u".to_string(), "const".to_string(), "notes".to_string()];
    ctl.clono_print_opt.lvars = vec!["datasets".to_string(), "n".to_string()];

    ctl.clono_group_opt.min_group = 1;

    ctl.allele_alg_opt.min_mult = 4;
    ctl.allele_alg_opt.min_alt = 4;

    ctl.join_alg_opt.max_score = 1_000_000.0;
    ctl.join_alg_opt.merge_onesies = true; // should just kill this as an option

    ctl.join_print_opt.pfreq = 1_000_000_000;
    ctl.join_print_opt.quiet = true;

    ctl.parseable_opt.pchains = 4;

    ctl.onesie_mult = 10_000;

    // Pretest for consistency amongst TCR, BCR, GEX and META.  Also preparse GEX.

    let (mut have_tcr, mut have_bcr) = (false, false);
    let mut have_gex = false;
    let mut have_meta = false;
    let mut gex = String::new();
    let mut bc = String::new();
    let mut metas = Vec::<String>::new();
    let mut xcrs = Vec::<String>::new();
    for i in 1..args.len() {
        if args[i].starts_with("BI=") {
            have_bcr = true;
            have_gex = true;
        } else if args[i].starts_with("TCR=") {
            have_tcr = true;
        } else if args[i].starts_with("BCR=") {
            have_bcr = true;
        } else if args[i].starts_with("GEX=") {
            have_gex = true;
        } else if args[i].starts_with("META=") {
            have_meta = true;
        }
        if args[i].starts_with("GEX=") {
            gex = args[i].after("GEX=").to_string();
        }
        if args[i].starts_with("BC=") {
            bc = args[i].after("BC=").to_string();
        }
        if is_simple_arg(&args[i], "MARK_STATS") {
            ctl.gen_opt.mark_stats = true;
        }
        if is_simple_arg(&args[i], "MARKED_B") {
            ctl.clono_filt_opt.marked_b = true;
        }
    }
    if have_meta && (have_tcr || have_bcr || have_gex || bc.len() > 0) {
        eprintln!("\nIf META is specified, then none of TCR, BCR, GEX or BC can be specified.\n");
        std::process::exit(1);
    }
    if have_tcr && have_bcr {
        eprintln!("\nKindly please do not specify both TCR and BCR.\n");
        std::process::exit(1);
    }
    let mut using_plot = false;

    // Preprocess BI argument.

    if ctl.gen_opt.internal_run {
        for i in 1..args.len() {
            if args[i].starts_with("BI=") {
                let n = args[i].after("BI=");
                if n != "m1" {
                    if !n.parse::<usize>().is_ok() || n.force_usize() < 1 || n.force_usize() > 13 {
                        eprintln!("\nBI=n only works if 1 <= n <= 13, or n = m1.\n");
                        std::process::exit(1);
                    }
                }
                let mut args2 = Vec::<String>::new();
                for j in 0..i {
                    args2.push(args[j].clone());
                }
                let f = include_str!["enclone.testdata.bcr.gex"];
                let mut found = false;
                for s in f.lines() {
                    if s == format!("DONOR={}", n) {
                        found = true;
                    } else if found && s.starts_with("DONOR=") {
                        break;
                    }
                    if found {
                        if s.starts_with("BCR=") || s.starts_with("GEX=") {
                            args2.push(s.to_string());
                        }
                        if s.starts_with("GEX=") {
                            gex = s.after("GEX=").to_string();
                        }
                        if s == "SPECIES=mouse" {
                            args2.push("MOUSE".to_string());
                        }
                    }
                }
                for j in i + 1..args.len() {
                    args2.push(args[j].clone());
                }
                args = args2;
                break;
            }
        }
    }

    // Define arguments that set something to true.

    let mut simple_set = vec![
        ("ACCEPT_INCONSISTENT", &mut ctl.gen_opt.accept_inconsistent),
        ("ACCEPT_REUSE", &mut ctl.gen_opt.accept_reuse),
        ("ANN", &mut ctl.join_print_opt.ann),
        ("ANN0", &mut ctl.join_print_opt.ann0),
        ("BARCODES", &mut ctl.clono_print_opt.barcodes),
        ("BASELINE", &mut ctl.gen_opt.baseline),
        ("BCJOIN", &mut ctl.join_alg_opt.bcjoin),
        ("CDIFF", &mut ctl.clono_filt_opt.cdiff),
        ("CHAIN_BRIEF", &mut ctl.clono_print_opt.chain_brief),
        ("CON", &mut ctl.allele_print_opt.con),
        ("CON_CON", &mut ctl.gen_opt.con_con),
        ("CON_TRACE", &mut ctl.allele_print_opt.con_trace),
        ("CURRENT_REF", &mut ctl.gen_opt.current_ref),
        ("DEBUG_TABLE_PRINTING", &mut ctl.debug_table_printing),
        ("DEL", &mut ctl.clono_filt_opt.del),
        ("DESCRIP", &mut ctl.gen_opt.descrip),
        ("EASY", &mut ctl.join_alg_opt.easy),
        ("ECHO", &mut ctl.gen_opt.echo),
        ("EXP", &mut ctl.gen_opt.exp),
        ("FORCE", &mut ctl.force),
        ("FULL_SEQC", &mut ctl.clono_print_opt.full_seqc),
        ("GRAPH", &mut ctl.gen_opt.graph),
        ("GROUP_HEAVY_CDR3", &mut ctl.clono_group_opt.heavy_cdr3_aa),
        ("GROUP_VJ_REFNAME", &mut ctl.clono_group_opt.vj_refname),
        ("HAVE_ONESIE", &mut ctl.clono_filt_opt.have_onesie),
        ("HEAVY_CHAIN_REUSE", &mut ctl.gen_opt.heavy_chain_reuse),
        ("IMGT", &mut ctl.gen_opt.imgt),
        ("IMGT_FIX", &mut ctl.gen_opt.imgt_fix),
        ("INDELS", &mut ctl.gen_opt.indels),
        ("INSERTIONS", &mut ctl.gen_opt.insertions),
        ("JC1", &mut ctl.gen_opt.jc1),
        ("KEEP_IMPROPER", &mut ctl.merge_all_impropers),
        ("MARKED", &mut ctl.clono_filt_opt.marked),
        ("MEAN", &mut ctl.clono_print_opt.mean),
        ("MIX_DONORS", &mut ctl.clono_filt_opt.donor),
        ("MOUSE", &mut ctl.gen_opt.mouse),
        ("NCELL", &mut ctl.gen_opt.ncell),
        ("NCROSS", &mut ctl.clono_filt_opt.ncross),
        ("NGEX", &mut ctl.clono_filt_opt.ngex),
        ("NGRAPH_FILTER", &mut ctl.gen_opt.ngraph_filter),
        ("NGROUP", &mut ctl.gen_opt.ngroup),
        ("NOPRINT", &mut ctl.gen_opt.noprint),
        ("NOTE_SIMPLE", &mut ctl.clono_print_opt.note_simple),
        ("NPLAIN", &mut ctl.pretty),
        ("NWHITEF", &mut ctl.gen_opt.nwhitef),
        ("NWARN", &mut ctl.gen_opt.nwarn),
        ("PCELL", &mut ctl.parseable_opt.pbarcode),
        ("PER_CELL", &mut ctl.clono_print_opt.bu),
        ("PROTECT_BADS", &mut ctl.clono_filt_opt.protect_bads),
        ("RE", &mut ctl.gen_opt.reannotate),
        ("REUSE", &mut ctl.gen_opt.reuse),
        ("SEQC", &mut ctl.clono_print_opt.seqc),
        ("SHOW_BC", &mut ctl.join_print_opt.show_bc),
        ("STABLE_DOC", &mut ctl.gen_opt.stable_doc),
        ("SUM", &mut ctl.clono_print_opt.sum),
        ("SUMMARY", &mut ctl.gen_opt.summary),
        ("SUMMARY_CLEAN", &mut ctl.gen_opt.summary_clean),
        ("SUMMARY_CSV", &mut ctl.gen_opt.summary_csv),
        ("TOY", &mut ctl.toy),
        ("UMI_FILT_MARK", &mut ctl.clono_filt_opt.umi_filt_mark),
        (
            "UMI_RATIO_FILT_MARK",
            &mut ctl.clono_filt_opt.umi_ratio_filt_mark,
        ),
        ("UTR_CON", &mut ctl.gen_opt.utr_con),
        ("VDUP", &mut ctl.clono_filt_opt.vdup),
        ("WEAK", &mut ctl.gen_opt.weak),
        ("WHITEF", &mut ctl.clono_filt_opt.whitef),
    ];

    // Define arguments that set something to a usize.

    let usize_set = [
        ("CHAINS_EXACT", &mut ctl.gen_opt.chains_exact),
        ("MAX_DATASETS", &mut ctl.clono_filt_opt.max_datasets),
        ("MIN_ALT", &mut ctl.allele_alg_opt.min_alt),
        ("MIN_CELLS_EXACT", &mut ctl.gen_opt.min_cells_exact),
        ("MIN_CHAINS_EXACT", &mut ctl.gen_opt.min_chains_exact),
        ("MIN_DATASETS", &mut ctl.clono_filt_opt.min_datasets),
        ("MIN_EXACTS", &mut ctl.clono_filt_opt.min_exacts),
        ("MIN_GROUP", &mut ctl.clono_group_opt.min_group),
        ("MIN_MULT", &mut ctl.allele_alg_opt.min_mult),
        ("MIN_UMI", &mut ctl.clono_filt_opt.min_umi),
        ("ONESIE_MULT", &mut ctl.onesie_mult),
        ("PCHAINS", &mut ctl.parseable_opt.pchains),
    ];

    // Traverse arguments.

    'args_loop: for i in 1..args.len() {
        let mut arg = args[i].to_string();

        // Strip out certain quoted expressions.

        if arg.contains("=\"") && arg.ends_with("\"") {
            let mut quotes = 0;
            for c in arg.chars() {
                if c == '\"' {
                    quotes += 1;
                }
            }
            if quotes == 2 {
                arg = format!("{}={}", arg.before("="), arg.between("\"", "\""));
            }
        }

        // Check for weird case that might arise if testing code is screwed up.

        if arg.len() == 0 {
            eprintln!(
                "\nYou've passed a null argument to enclone.  Normally that isn't \
                 possible.\nPlease take a detailed look at how you're invoking enclone.\n"
            );
            std::process::exit(1);
        }

        // Process simple set arguments.

        for j in 0..simple_set.len() {
            if arg == simple_set[j].0.to_string() {
                *(simple_set[j].1) = true;
                continue 'args_loop;
            }
        }

        // Process usize args.

        for j in 0..usize_set.len() {
            if is_usize_arg(&arg, &usize_set[j].0) {
                *(usize_set[j].1) = arg.after(&format!("{}=", usize_set[j].0)).force_usize();
                continue 'args_loop;
            }
        }

        // Process the argument.

        if is_simple_arg(&arg, "SEQ") {
            ctl.join_print_opt.seq = true;

        // Not movable.
        } else if is_simple_arg(&arg, "H5") {
            ctl.gen_opt.force_h5 = true;
        } else if arg == "LEGEND" {
            ctl.gen_opt.use_legend = true;
        } else if is_usize_arg(&arg, "REQUIRED_FPS") {
            ctl.gen_opt.required_fps = Some(arg.after("REQUIRED_FPS=").force_usize());
        } else if is_usize_arg(&arg, "EXACT") {
            ctl.gen_opt.exact = Some(arg.after("EXACT=").force_usize());
        } else if is_usize_arg(&arg, "MIN_CHAINS") {
            ctl.clono_filt_opt.min_chains = arg.after("MIN_CHAINS=").force_usize();
        } else if is_usize_arg(&arg, "MAX_CHAINS") {
            ctl.clono_filt_opt.max_chains = arg.after("MAX_CHAINS=").force_usize();
        } else if arg == "PRINT_CPU" {
        } else if arg == "PRINT_CPU_INFO" {

            // Other.
        } else if is_simple_arg(&arg, "DUMP_INTERNAL_IDS") {
        } else if is_simple_arg(&arg, "COMP") {
        } else if is_simple_arg(&arg, "COMP2") {
        } else if is_simple_arg(&arg, "LONG_HELP") {
        } else if is_simple_arg(&arg, "FAIL_ONLY=true") {
            ctl.clono_filt_opt.fail_only = true;
        } else if arg.starts_with("BI=") {
            continue;
        } else if is_simple_arg(&arg, "MARK_STATS") {
        } else if arg == "HTML" || arg.starts_with("HTML=") {
        } else if arg == "SVG" {
        } else if arg.starts_with("LEGEND=") {
            let x = parse_csv(&arg.after("LEGEND="));
            if x.len() == 0 || x.len() % 2 != 0 {
                eprintln!("\nValue of LEGEND doesn't make sense.\n");
                std::process::exit(1);
            }
            ctl.gen_opt.use_legend = true;
            for i in 0..x.len() / 2 {
                ctl.gen_opt
                    .legend
                    .push((x[2 * i].clone(), x[2 * i + 1].clone()));
            }
        } else if is_simple_arg(&arg, "NH5") {
            ctl.gen_opt.force_h5 = false;
        } else if is_simple_arg(&arg, "H5_SLICE") {
            ctl.gen_opt.h5_pre = false;
        } else if is_simple_arg(&arg, "CTRLC") {
        } else if is_simple_arg(&arg, "CELLRANGER") {
        } else if is_simple_arg(&arg, "FORCE_EXTERNAL") {
        } else if arg.starts_with("BINARY=") {
            ctl.gen_opt.binary = arg.after("BINARY=").to_string();
        } else if arg.starts_with("PROTO=") {
            ctl.gen_opt.proto = arg.after("PROTO=").to_string();
        } else if is_simple_arg(&arg, "PRINT_FAILED_JOINS") {
            ctl.join_print_opt.quiet = false;
        } else if arg.starts_with("BARCODE=") {
            let bcs = arg.after("BARCODE=").split(',').collect::<Vec<&str>>();
            let mut x = Vec::<String>::new();
            for j in 0..bcs.len() {
                if !bcs[j].contains('-') {
                    eprintln!(
                        "\nValue for a barcode in BARCODE argument is invalid, must contain -.\n"
                    );
                    std::process::exit(1);
                }
                x.push(bcs[j].to_string());
            }
            ctl.clono_filt_opt.barcode = x;
        } else if is_simple_arg(&arg, "NUMI") {
            ctl.clono_filt_opt.umi_filt = false;
        } else if is_simple_arg(&arg, "NUMI_RATIO") {
            ctl.clono_filt_opt.umi_ratio_filt = false;
        } else if is_simple_arg(&arg, "MARKED_B") {
        } else if is_simple_arg(&arg, "NWEAK_CHAINS") {
            ctl.clono_filt_opt.weak_chains = false;
        } else if is_simple_arg(&arg, "NWEAK_ONESIES") {
            ctl.clono_filt_opt.weak_onesies = false;
        } else if is_simple_arg(&arg, "NFOURSIE_KILL") {
            ctl.clono_filt_opt.weak_foursies = false;
        } else if is_simple_arg(&arg, "NBC_DUP") {
            ctl.clono_filt_opt.bc_dup = false;
        } else if arg.starts_with("F=") {
            let filt = arg.after("F=").to_string();
            ctl.clono_filt_opt.bounds.push(LinearCondition::new(&filt));
        } else if arg.starts_with("SCAN=") {
            let mut x = arg.after("SCAN=").to_string();
            x = x.replace(" ", "").to_string();
            let x = x.split(',').collect::<Vec<&str>>();
            if x.len() != 3 {
                eprintln!("\nArgument to SCAN must have three components.\n");
                std::process::exit(1);
            }
            ctl.gen_opt.gene_scan_test = Some(LinearCondition::new(&x[0]));
            ctl.gen_opt.gene_scan_control = Some(LinearCondition::new(&x[1]));
            let threshold = LinearCondition::new(&x[2]);
            for i in 0..threshold.var.len() {
                if threshold.var[i] != "t".to_string() && threshold.var[i] != "c".to_string() {
                    eprintln!("\nIllegal variable in threshold for scan.\n");
                    std::process::exit(1);
                }
            }
            ctl.gen_opt.gene_scan_threshold = Some(threshold);
        } else if arg.starts_with("PLOT=") {
            using_plot = true;
            let x = arg.after("PLOT=").split(',').collect::<Vec<&str>>();
            if x.is_empty() {
                eprintln!("\nArgument to PLOT is invalid.\n");
                std::process::exit(1);
            }
            ctl.gen_opt.plot_file = x[0].to_string();
            for j in 1..x.len() {
                if !x[j].contains("->") {
                    eprintln!("\nArgument to PLOT is invalid.\n");
                    std::process::exit(1);
                }
                ctl.gen_opt
                    .sample_color_map
                    .insert(x[j].before("->").to_string(), x[j].after("->").to_string());
            }
        } else if arg.starts_with("PLOT_BY_ISOTYPE=") {
            ctl.gen_opt.plot_by_isotype = true;
            ctl.gen_opt.plot_file = arg.after("PLOT_BY_ISOTYPE=").to_string();
            if ctl.gen_opt.plot_file.is_empty() {
                eprintln!("\nFilename value needs to be supplied to PLOT_BY_ISOTYPE.\n");
                std::process::exit(1);
            }
        } else if arg.starts_with("PLOT_BY_MARK=") {
            ctl.gen_opt.plot_by_mark = true;
            ctl.gen_opt.plot_file = arg.after("PLOT_BY_MARK=").to_string();
            if ctl.gen_opt.plot_file.is_empty() {
                eprintln!("\nFilename value needs to be supplied to PLOT_BY_MARK.\n");
                std::process::exit(1);
            }
        } else if arg.starts_with("EMAIL=") {
        } else if arg.starts_with("REF=") {
            ctl.gen_opt.refname = arg.after("REF=").to_string();
        } else if is_simple_arg(&arg, "NSILENT") {
            ctl.silent = false;
        } else if is_simple_arg(&arg, "FAIL_ONLY=false") {
            ctl.clono_filt_opt.fail_only = false;
        } else if is_simple_arg(&arg, "NQUAL") {
            ctl.clono_filt_opt.qual_filter = false;
        } else if is_simple_arg(&arg, "NOPAGER") {
        } else if arg.starts_with("POUT=") {
            ctl.parseable_opt.pout = arg.after("POUT=").to_string();
        } else if arg.starts_with("DONOR_REF_FILE=") {
            ctl.gen_opt.dref_file = arg.after("DONOR_REF_FILE=").to_string();
        } else if arg.starts_with("EXT=") {
            ctl.gen_opt.ext = arg.after("EXT=").to_string();
        } else if arg.starts_with("TRACE_BARCODE=") {
            ctl.gen_opt.trace_barcode = arg.after("TRACE_BARCODE=").to_string();
        } else if is_usize_arg(&arg, "MAX_CORES") {
            let nthreads = arg.after("MAX_CORES=").force_usize();
            let _ = rayon::ThreadPoolBuilder::new()
                .num_threads(nthreads)
                .build_global();
        } else if arg.starts_with("PCOLS=") {
            ctl.parseable_opt.pcols.clear();
            let p = arg.after("PCOLS=").split(',').collect::<Vec<&str>>();
            for i in 0..p.len() {
                let mut x = p[i].to_string();
                x = x.replace("_sum", "_Σ");
                x = x.replace("_mean", "_μ");
                ctl.parseable_opt.pcols.push(x.to_string());
                ctl.parseable_opt.pcols_sort = ctl.parseable_opt.pcols.clone();
                ctl.parseable_opt.pcols_sortx = ctl.parseable_opt.pcols.clone();
                for j in 0..ctl.parseable_opt.pcols_sortx.len() {
                    if ctl.parseable_opt.pcols_sortx[j].contains(":") {
                        ctl.parseable_opt.pcols_sortx[j] =
                            ctl.parseable_opt.pcols_sortx[j].before(":").to_string();
                    }
                }
                unique_sort(&mut ctl.parseable_opt.pcols_sort);
                unique_sort(&mut ctl.parseable_opt.pcols_sortx);
            }
        } else if is_simple_arg(&arg, "PLAIN") {
        } else if is_simple_arg(&arg, "NOPRETTY") {
        } else if arg.starts_with("VJ=") {
            ctl.clono_filt_opt.vj = arg.after("VJ=").as_bytes().to_vec();
            for c in ctl.clono_filt_opt.vj.iter() {
                if !(*c == b'A' || *c == b'C' || *c == b'G' || *c == b'T') {
                    eprintln!("\nIllegal value for VJ, must be over alphabet ACGT.\n");
                    std::process::exit(1);
                }
            }
        } else if arg.starts_with("AMINO=") {
            ctl.clono_print_opt.amino.clear();
            for x in arg.after("AMINO=").split(',').collect::<Vec<&str>>() {
                if x != "" {
                    ctl.clono_print_opt.amino.push(x.to_string());
                }
            }
            for x in ctl.clono_print_opt.amino.iter() {
                let mut ok = false;
                if *x == "cdr3" || *x == "var" || *x == "share" || *x == "donor" || *x == "donorn" {
                    ok = true;
                } else if x.contains('-') {
                    let (start, stop) = (x.before("-"), x.after("-"));
                    if start.parse::<usize>().is_ok() && stop.parse::<usize>().is_ok() {
                        if start.force_usize() <= stop.force_usize() {
                            ok = true;
                        }
                    }
                }
                if !ok {
                    eprintln!(
                        "\nUnrecognized variable {} for AMINO.  Please type \
                         \"enclone help amino\".\n",
                        x
                    );
                    std::process::exit(1);
                }
            }
        } else if arg.starts_with("CVARS=") {
            ctl.clono_print_opt.cvars.clear();
            for x in arg.after("CVARS=").split(',').collect::<Vec<&str>>() {
                if x.len() > 0 {
                    ctl.clono_print_opt.cvars.push(x.to_string());
                }
            }
            for x in ctl.clono_print_opt.cvars.iter_mut() {
                *x = x.replace("_sum", "_Σ");
                *x = x.replace("_mean", "_μ");
            }
        } else if arg.starts_with("CVARSP=") {
            for x in arg.after("CVARSP=").split(',').collect::<Vec<&str>>() {
                if x.len() > 0 {
                    ctl.clono_print_opt.cvars.push(x.to_string());
                }
            }
            for x in ctl.clono_print_opt.cvars.iter_mut() {
                *x = x.replace("_sum", "_Σ");
                *x = x.replace("_mean", "_μ");
            }
        } else if arg.starts_with("LVARS=") {
            ctl.clono_print_opt.lvars.clear();
            for x in arg.after("LVARS=").split(',').collect::<Vec<&str>>() {
                ctl.clono_print_opt.lvars.push(x.to_string());
            }
            for x in ctl.clono_print_opt.lvars.iter_mut() {
                *x = x.replace("_sum", "_Σ");
                *x = x.replace("_mean", "_μ");
            }
        } else if arg.starts_with("LVARSP=") {
            let lvarsp = arg.after("LVARSP=").split(',').collect::<Vec<&str>>();
            for x in lvarsp {
                ctl.clono_print_opt.lvars.push(x.to_string());
            }
            for x in ctl.clono_print_opt.lvars.iter_mut() {
                *x = x.replace("_sum", "_Σ");
                *x = x.replace("_mean", "_μ");
            }
        } else if is_f64_arg(&arg, "MAX_SCORE") {
            ctl.join_alg_opt.max_score = arg.after("MAX_SCORE=").force_f64();
        } else if arg.starts_with("EXFASTA=") {
            ctl.gen_opt.fasta = arg.after("EXFASTA=").to_string();
        } else if arg.starts_with("FASTA=") {
            ctl.gen_opt.fasta_filename = arg.after("FASTA=").to_string();
        } else if arg.starts_with("FASTA_AA=") {
            ctl.gen_opt.fasta_aa_filename = arg.after("FASTA_AA=").to_string();
        } else if arg.starts_with("CDR3=") {
            let reg = Regex::new(&format!("^{}$", arg.after("CDR3=")));
            if !reg.is_ok() {
                eprintln!(
                    "\nYour CDR3 value {} could not be parsed as a regular expression.\n",
                    arg.after("CDR3=")
                );
                std::process::exit(1);
            }
            ctl.clono_filt_opt.cdr3 = Some(reg.unwrap());
        } else if arg.starts_with("GEX=") {
        } else if arg.starts_with("BC=") {
        } else if is_usize_arg(&arg, "CHAINS") {
            ctl.clono_filt_opt.min_chains = arg.after("CHAINS=").force_usize();
            ctl.clono_filt_opt.max_chains = arg.after("CHAINS=").force_usize();
        } else if arg.starts_with("SEG=") {
            let fields = arg.after("SEG=").split('|').collect::<Vec<&str>>();
            for x in fields.iter() {
                ctl.clono_filt_opt.seg.push(x.to_string());
            }
            ctl.clono_filt_opt.seg.sort();
        } else if arg.starts_with("SEGN=") {
            let fields = arg.after("SEGN=").split('|').collect::<Vec<&str>>();
            for x in fields.iter() {
                if !x.parse::<i32>().is_ok() {
                    eprintln!("\nInvalid argument to SEGN.\n");
                    std::process::exit(1);
                }
                ctl.clono_filt_opt.segn.push(x.to_string());
            }
            ctl.clono_filt_opt.segn.sort();
        } else if arg.starts_with("PRE=") {
        } else if is_usize_arg(&arg, "MIN_CELLS") {
            ctl.clono_filt_opt.ncells_low = arg.after("MIN_CELLS=").force_usize();
        } else if is_usize_arg(&arg, "MAX_CELLS") {
            ctl.clono_filt_opt.ncells_high = arg.after("MAX_CELLS=").force_usize();
        } else if is_usize_arg(&arg, "CELLS") {
            ctl.clono_filt_opt.ncells_low = arg.after("CELLS=").force_usize();
            ctl.clono_filt_opt.ncells_high = ctl.clono_filt_opt.ncells_low;
        } else if is_usize_arg(&arg, "PFREQ") {
            ctl.join_print_opt.pfreq = arg.after("PFREQ=").force_usize();
        } else if arg.starts_with("HAPS=") {
            // done above
        } else if arg.starts_with("META=") {
            let f = arg.after("META=");
            metas.push(f.to_string());
        } else if arg.starts_with("TCR=")
            || arg.starts_with("BCR=")
            || (arg.len() > 0 && arg.as_bytes()[0] >= b'0' && arg.as_bytes()[0] <= b'9')
        {
            xcrs.push(arg.to_string());
        } else {
            eprintln!("\nUnrecognized argument {}.\n", arg);
            std::process::exit(1);
        }
    }
    check_cvars(&ctl);
    if metas.len() > 0 {
        let f = &metas[metas.len() - 1];
        let f = get_path_fail(&f, &ctl, "META");
        proc_meta(&f, &mut ctl);
    }
    if xcrs.len() > 0 {
        let arg = &xcrs[xcrs.len() - 1];
        proc_xcr(&arg, &gex, &bc, have_gex, &mut ctl);
    }
    for i in 0..ctl.sample_info.n() {
        let (mut cells_cr, mut rpc_cr) = (None, None);
        if ctl.gen_opt.internal_run {
            let p = &ctl.sample_info.dataset_path[i];
            let mut f = format!("{}/metrics_summary_csv.csv", p);
            if !path_exists(&f) {
                f = format!("{}/metrics_summary.csv", p);
            }
            if path_exists(&f) {
                let f = open_for_read![&f];
                let mut count = 0;
                let (mut cells_field, mut rpc_field) = (None, None);
                for line in f.lines() {
                    count += 1;
                    let s = line.unwrap();
                    let fields = parse_csv(&s);
                    for (i, x) in fields.iter().enumerate() {
                        if count == 1 {
                            if *x == "Estimated Number of Cells" {
                                cells_field = Some(i);
                            } else if *x == "Mean Read Pairs per Cell" {
                                rpc_field = Some(i);
                            }
                        } else if count == 2 {
                            if Some(i) == cells_field {
                                let mut n = x.to_string();
                                if n.contains("\"") {
                                    n = n.between("\"", "\"").to_string();
                                }
                                n = n.replace(",", "");
                                cells_cr = Some(n.force_usize());
                            } else if Some(i) == rpc_field {
                                let mut n = x.to_string();
                                if n.contains("\"") {
                                    n = n.between("\"", "\"").to_string();
                                }
                                n = n.replace(",", "");
                                rpc_cr = Some(n.force_usize());
                            }
                        }
                    }
                }
            }
        }
        ctl.sample_info.cells_cellranger.push(cells_cr);
        ctl.sample_info
            .mean_read_pairs_per_cell_cellranger
            .push(rpc_cr);
    }
    if ctl.gen_opt.plot_by_isotype {
        if using_plot || ctl.gen_opt.use_legend {
            eprintln!("\nPLOT_BY_ISOTYPE cannot be used with PLOT or LEGEND.\n");
            std::process::exit(1);
        }
        if !have_bcr {
            eprintln!("\nPLOT_BY_ISOTYPE can only be used with BCR data.\n");
            std::process::exit(1);
        }
        if ctl.gen_opt.plot_by_mark {
            eprintln!("\nPLOT_BY_ISOTYPE and PLOT_BY_MARK cannot be used together.\n");
            std::process::exit(1);
        }
    }
    if ctl.gen_opt.plot_by_mark {
        if using_plot || ctl.gen_opt.use_legend {
            eprintln!("\nPLOT_BY_MARK cannot be used with PLOT or LEGEND.\n");
            std::process::exit(1);
        }
    }
    if ctl.parseable_opt.pbarcode && ctl.parseable_opt.pout.len() == 0 {
        eprintln!("\nIt does not make sense to specify PCELL unless POUT is also specified.\n");
        std::process::exit(1);
    }
    if ctl.sample_info.n() == 0 {
        eprintln!("\nNo TCR or BCR data have been specified.\n");
        std::process::exit(1);
    }
    let mut donors = Vec::<String>::new();
    let mut samples = Vec::<String>::new();
    let mut tags = Vec::<String>::new();
    let mut sample_for_bc = Vec::<String>::new();
    let mut donor_for_bc = Vec::<String>::new();
    for i in 0..ctl.sample_info.n() {
        for x in ctl.sample_info.sample_for_bc[i].iter() {
            samples.push(x.1.clone());
            sample_for_bc.push(x.1.clone());
        }
        for x in ctl.sample_info.donor_for_bc[i].iter() {
            donors.push(x.1.clone());
            donor_for_bc.push(x.1.clone());
        }
        for x in ctl.sample_info.tag[i].iter() {
            tags.push((x.1).clone());
        }
        donors.push(ctl.sample_info.donor_id[i].clone());
        samples.push(ctl.sample_info.sample_id[i].clone());
    }
    unique_sort(&mut donors);
    unique_sort(&mut samples);
    unique_sort(&mut tags);
    unique_sort(&mut sample_for_bc);
    unique_sort(&mut donor_for_bc);
    ctl.sample_info.donors = donors.len();
    ctl.sample_info.dataset_list = ctl.sample_info.dataset_id.clone();
    unique_sort(&mut ctl.sample_info.dataset_list);
    ctl.sample_info.sample_list = samples.clone();
    ctl.sample_info.donor_list = donors.clone();
    ctl.sample_info.tag_list = tags;
    for i in 0..ctl.sample_info.donor_for_bc.len() {
        if ctl.sample_info.donor_for_bc[i].len() > 0 {
            ctl.clono_filt_opt.donor = true;
        }
    }
    if ctl.comp2 {
        println!("\n-- used {:.2} seconds processing args", elapsed(&targs));
    }
    proc_args_tail(&mut ctl, &args);

    // Check for invalid variables in linear conditions.

    for i in 0..ctl.clono_filt_opt.bounds.len() {
        ctl.clono_filt_opt.bounds[i].require_valid_variables(&ctl);
    }
    if ctl.gen_opt.gene_scan_test.is_some() {
        ctl.gen_opt
            .gene_scan_test
            .as_ref()
            .unwrap()
            .require_valid_variables(&ctl);
        ctl.gen_opt
            .gene_scan_control
            .as_ref()
            .unwrap()
            .require_valid_variables(&ctl);
    }
}