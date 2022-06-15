// Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
//
// See also public_light_chain_analysis.rs.
//
// Analyze light chains.  Supply a single file of data, with one line per cell, and fields
// including donors_cell,v_name1,v_name2,dref,cdr3_aa1,clonotype_ncells,const1.
//
// Data from:
//
// enclone BCR=@test BUILT_IN CHAINS_EXACT=2 CHAINS=2 NOPRINT POUT=stdout PCELL ECHOC
//         PCOLS=donors_cell,v_name1,v_name2,dref,cdr3_aa1,clonotype_ncells,const1,hcomp
//         > per_cell_stuff
//
// Optional arguments:
//
// SHOW -- Instead print data for pairs of cells from the same donor at 100% identity with
// dref1 > 0 and dref2 > 0 and having the same light chain gene.
// For this, the order of output lines is nondeterministic.  Designed for use with J option.
//
// J -- Require different J genes rather than different V genes.
//
// JPLUS -- J, and also require that there are at least three positions in the last 25 bases of
// the different J gene reference sequences at which the reference sequences differ and the two
// cells both agree with their assigned J gene reference.
//
// SOLO -- reduce to one cell per clonotype.
//
// SAME -- instead require same heavy chain gene, and different light chain genes.

use enclone_core::hcat;
use io_utils::*;
use pretty_trace::PrettyTrace;
use rayon::prelude::*;
use std::collections::HashMap;
use std::env;
use std::io::{BufRead, Write};
use string_utils::strme;
use string_utils::TextUtils;
use tables::print_tabular_vbox;
use vector_utils::erase_if;

fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();
    let f = open_for_read![&args[1]];
    let mut show = false;
    // let mut use_j = false;
    // let mut jplus = false;
    let mut solo = false;
    let mut opt_same = false;
    for i in 2..args.len() {
        if args[i] == "SHOW" {
            show = true;
        /*
        } else if args[i] == "J" {
            use_j = true;
        */
        } else if args[i] == "SOLO" {
            solo = true;
        } else if args[i] == "SAME" {
            opt_same = true;
        /*
        } else if args[i] == "JPLUS" {
            jplus = true;
            use_j = true;
        */
        } else {
            eprintln!("\nIllegal argument.\n");
            std::process::exit(1);
        }
    }
    let mut first = true;
    let mut tof = HashMap::<String, usize>::new();

    #[derive(Eq, Ord, PartialEq, PartialOrd)]
    struct CellData {
        donor: String,
        cdr3_len1: usize,
        cdr3_aa1: Vec<u8>,
        v_name1: String,
        v_name2: String,
        dref: usize,
        dataset: String,
        barcode: String,
        j_name1: String,
        fwr4_dna_ref1: Vec<u8>,
        fwr4_dna1: Vec<u8>,
        v_name2_orig: String,
        j_name2: String,
        cdr3_aa2: Vec<u8>,
        leader_dna_ref2: Vec<u8>,
        leader_dna2: Vec<u8>,
        fwr1_dna_ref2: Vec<u8>,
        fwr1_dna2: Vec<u8>,
        fwr2_dna_ref2: Vec<u8>,
        fwr2_dna2: Vec<u8>,
        fwr3_dna_ref2: Vec<u8>,
        fwr3_dna2: Vec<u8>,
        fwr4_dna_ref2: Vec<u8>,
        fwr4_dna2: Vec<u8>,
        cdr1_dna_ref2: Vec<u8>,
        cdr1_dna2: Vec<u8>,
        cdr2_dna_ref2: Vec<u8>,
        cdr2_dna2: Vec<u8>,
        cdr3_dna2: Vec<u8>,
        cdr3_dna1: Vec<u8>,

    }

    let mut data = Vec::new();
    let mut clonotype = Vec::<usize>::new();
    for line in f.lines() {
        let s = line.unwrap();
        if s.starts_with("#") {
            continue;
        }
        let fields = s.split(',').collect::<Vec<&str>>();
        if first {
            for i in 0..fields.len() {
                tof.insert(fields[i].to_string(), i);
            }
            first = false;
        } else {
            data.push( CellData {
                donor: fields[tof["donors_cell"]].to_string(),
                cdr3_len1: fields[tof["cdr3_aa1"]].len(),
                cdr3_aa1: fields[tof["cdr3_aa1"]].to_string().as_bytes().to_vec(),
                v_name1: fields[tof["v_name1"]].to_string(),
                v_name2: fields[tof["v_name2"]].to_string(),
                dref: fields[tof["dref"]].force_usize(),
                dataset: fields[tof["datasets_cell"]].to_string(),
                barcode: fields[tof["barcode"]].to_string(),
                j_name1: fields[tof["j_name1"]].to_string(),
                fwr4_dna_ref1: fields[tof["fwr4_dna_ref1"]].to_string().as_bytes().to_vec(),
                fwr4_dna1: fields[tof["fwr4_dna1"]].to_string().as_bytes().to_vec(),
                v_name2_orig: fields[tof["v_name2"]].to_string(),
                j_name2: fields[tof["j_name2"]].to_string(),
                cdr3_aa2: fields[tof["cdr3_aa2"]].to_string().as_bytes().to_vec(),
                leader_dna_ref2: fields[tof["leader_dna_ref2"]].to_string().as_bytes().to_vec(),
                leader_dna2: fields[tof["leader_dna2"]].to_string().as_bytes().to_vec(),
                fwr1_dna_ref2: fields[tof["fwr1_dna_ref2"]].to_string().as_bytes().to_vec(),
                fwr1_dna2: fields[tof["fwr1_dna2"]].to_string().as_bytes().to_vec(),
                fwr2_dna_ref2: fields[tof["fwr2_dna_ref2"]].to_string().as_bytes().to_vec(),
                fwr2_dna2: fields[tof["fwr2_dna2"]].to_string().as_bytes().to_vec(),
                fwr3_dna_ref2: fields[tof["fwr3_dna_ref2"]].to_string().as_bytes().to_vec(),
                fwr3_dna2: fields[tof["fwr3_dna2"]].to_string().as_bytes().to_vec(),
                fwr4_dna_ref2: fields[tof["fwr4_dna_ref2"]].to_string().as_bytes().to_vec(),
                fwr4_dna2: fields[tof["fwr4_dna2"]].to_string().as_bytes().to_vec(),
                cdr1_dna_ref2: fields[tof["cdr1_dna_ref2"]].to_string().as_bytes().to_vec(),
                cdr1_dna2: fields[tof["cdr1_dna2"]].to_string().as_bytes().to_vec(),
                cdr2_dna_ref2: fields[tof["cdr2_dna_ref2"]].to_string().as_bytes().to_vec(),
                cdr2_dna2: fields[tof["cdr2_dna2"]].to_string().as_bytes().to_vec(),
                cdr3_dna2: fields[tof["cdr3_dna2"]].to_string().as_bytes().to_vec(),
                cdr3_dna1: fields[tof["cdr3_dna1"]].to_string().as_bytes().to_vec(),
            });
            clonotype.push(fields[tof["group_id"]].force_usize());
        }
    }

    // In solo case, reduce to one cell per clonotype.

    if solo {
        let mut to_delete = vec![false; data.len()];
        for i in 0..clonotype.len() {
            if i > 0 && clonotype[i] == clonotype[i - 1] {
                to_delete[i] = true;
            }
        }
        erase_if(&mut data, &to_delete);
        erase_if(&mut clonotype, &to_delete);
    }

    // Sort.

    data.sort();

    // Replace paralogous gene names.

    for i in 0..data.len() {
        data[i].v_name2 = data[i].v_name2.replace("D", "");
    }

    // Define groups based on equal donor and CDR3H length.
    // Plus placeholder for results, see next.

    let mut bounds = Vec::new();
    let mut i = 0;
    while i < data.len() {
        let mut j = i + 1;
        while j < data.len() {
            if data[j].donor != data[i].donor || data[j].cdr3_len1 != data[i].cdr3_len1 {
                break;
            }
            j += 1;
        }
        bounds.push((i, j, vec![vec![(0, 0, 0, 0); 11]; 5], Vec::new()));
        i = j;
    }

    // Results = for each percent identity, rounded down:
    // 1. count for equal light chain gene names and dref1 = 0 and dref2 = 0
    // 2. count for unequal light chain gene names and dref1 = 0 and dref2 = 0
    // 3. count for equal light chain gene names and dref1 > 0 and dref2 > 0
    // 4. count for unequal light chain gene names and dref1 > 0 and dref2 > 0.
    //
    // Make one pass for each donor.

    bounds.par_iter_mut().for_each(|res| {
        let i = res.0;
        let j = res.1;
        let log = &mut res.3;
        let d = data[i].donor.after("d").force_usize() - 1;
        for k1 in i..j {
            for k2 in k1 + 1..j {

                // For now requiring same heavy chain gene.

                if data[k1].v_name1 != data[k2].v_name1 {
                    continue;
                }

                // ================================================================================

                // Test for validated.  Work in progress, several steps.

                let mut validated = false;
                let d1 = &data[k1];
                let d2 = &data[k2];
                const MIN_VAL: usize = 3;

                // Test if light chain CDR3s have different lengths.

                if d1.cdr3_aa2.len() != d2.cdr3_aa2.len() {
                    validated = true;
                }

                // Test if light chain FWR3s have different lengths.
                // Seems a bit fishy.

                if !validated {
                    if d1.fwr3_dna2.len() != d2.fwr3_dna2.len() {
                        validated = true;
                    }
                }

                // Test FWR4 on heavy chain.

                if !validated {
                    let mut supp = 0;
                    let (mut sup1, mut sup2) = (0, 0);
                    // not sure I trust the following:
                    if d1.fwr4_dna1.len() != d2.fwr4_dna1.len()
                        || d1.fwr4_dna_ref1.len() != d2.fwr4_dna_ref1.len() {
                        validated = true;
                    } else {
                        for i in 0..d1.fwr4_dna1.len() {
                            if d1.fwr4_dna_ref1[i] != d2.fwr4_dna_ref1[i] {
                                if d1.fwr4_dna1[i] == d1.fwr4_dna_ref1[i]
                                    && d2.fwr4_dna1[i] == d2.fwr4_dna_ref1[i] {
                                    supp += 1;
                                } else if d1.fwr4_dna1[i] == d1.fwr4_dna_ref1[i]
                                    && d2.fwr4_dna1[i] == d1.fwr4_dna_ref1[i] {
                                    sup1 += 1;
                                } else if d1.fwr4_dna1[i] == d2.fwr4_dna_ref1[i]
                                    && d2.fwr4_dna1[i] == d2.fwr4_dna_ref1[i] {
                                    sup2 += 1;
                                }
                            }
                        }
                        if supp >= MIN_VAL && sup1 == 0 && sup2 == 0 {
                            validated = true;
                        }
                    }
                }

                // Test FWR4 on light chain.

                if !validated {
                    let mut supp = 0;
                    let (mut sup1, mut sup2) = (0, 0);
                    // not sure I trust the following:
                    if d1.fwr4_dna2.len() != d2.fwr4_dna2.len()
                        || d1.fwr4_dna_ref2.len() != d2.fwr4_dna_ref2.len() {
                        validated = true;
                    } else {
                        for i in 0..d1.fwr4_dna2.len() {
                            if d1.fwr4_dna_ref2[i] != d2.fwr4_dna_ref2[i] {
                                if d1.fwr4_dna2[i] == d1.fwr4_dna_ref2[i]
                                    && d2.fwr4_dna2[i] == d2.fwr4_dna_ref2[i] {
                                    supp += 1;
                                } else if d1.fwr4_dna2[i] == d1.fwr4_dna_ref2[i]
                                    && d2.fwr4_dna2[i] == d1.fwr4_dna_ref2[i] {
                                    sup1 += 1;
                                } else if d1.fwr4_dna2[i] == d2.fwr4_dna_ref2[i]
                                    && d2.fwr4_dna2[i] == d2.fwr4_dna_ref2[i] {
                                    sup2 += 1;
                                }
                            }
                        }
                        if supp >= MIN_VAL && sup1 == 0 && sup2 == 0 {
                            validated = true;
                        }
                    }
                }

                if !validated {
                    continue;
                }

                // ================================================================================

                // For now requiring some reference difference.  Using gene names
                // now but that's not optimal.

                if data[k1].j_name1 == data[k2].j_name1 
                    && data[k1].v_name2 == data[k2].v_name2 
                    && data[k1].j_name2 == data[k2].j_name2 
                    && data[k1].cdr3_aa2.len() == data[k2].cdr3_aa2.len()  {
                    continue;
                }

                // Compute stuff.

                let mut same = 0;
                for m in 0..data[k1].cdr3_aa1.len() {
                    if data[k1].cdr3_aa1[m] == data[k2].cdr3_aa1[m] {
                        same += 1;
                    }
                }
                let ident = 100.0 * same as f64 / data[k1].cdr3_aa1.len() as f64;
                let ident = ident.floor() as usize;
                let ident = ident / 10;
                let (dref1, dref2) = (data[k1].dref, data[k2].dref);
                let eq_light = data[k1].v_name2 == data[k2].v_name2;

                if dref1 == 0 && dref2 == 0 {
                    if eq_light {
                        res.2[d + 1][ident].0 += 1;
                        res.2[0][ident].0 += 1;
                    } else {
                        res.2[d + 1][ident].1 += 1;
                        res.2[0][ident].1 += 1;
                    }
                } else if dref1 > 0 && dref2 > 0 {
                    if eq_light {
                        res.2[d + 1][ident].2 += 1;
                        res.2[0][ident].2 += 1;
                        if show && ident == 10 {
                            fwriteln!(log, "\n{}", strme(&vec![b'='; 100]));
                            fwriteln!(log,
                                "\n{} {} {} {} {} {}\n{} {} {} {} {} {}",
                                data[k1].dataset, data[k1].barcode, 
                                data[k1].j_name1,
                                data[k1].v_name2_orig, data[k1].j_name2, 
                                data[k1].cdr3_aa2.len(),
                                data[k2].dataset, data[k2].barcode,
                                data[k2].j_name1,
                                data[k2].v_name2_orig, data[k2].j_name2, 
                                data[k2].cdr3_aa2.len(),
                            );

                            // Process CDR3 on heavy chain.

                            let mut refdiffs = 0;
                            fwriteln!(log, "\nCDR3-heavy");
                            for i in 0..d1.cdr3_dna1.len() {
                                if d1.cdr3_dna1[i] == d2.cdr3_dna1[i] {
                                    fwrite!(log, " ");
                                } else {
                                    fwrite!(log, "*");
                                    refdiffs += 1;
                                }
                            }
                            fwriteln!(log, "");
                            fwriteln!(log, "{} seq1", strme(&d1.cdr3_dna1));
                            fwriteln!(log, "{} seq2", strme(&d2.cdr3_dna1));

                            // Process FWR4 on heavy chain.

                            let mut supp = 0;
                            let mut sup1 = 0;
                            let mut sup2 = 0;
                            if d1.fwr4_dna_ref1 != d2.fwr4_dna_ref1 {
                                fwriteln!(log, "\nFWR4-heavy");
                                for i in 0..d1.fwr4_dna1.len() {
                                    if d1.fwr4_dna_ref1[i] == d2.fwr4_dna_ref1[i] {
                                        fwrite!(log, " ");
                                    } else {
                                        fwrite!(log, "*");
                                        refdiffs += 1;
                                    }
                                }
                                fwriteln!(log, "");
                                fwriteln!(log, "{} ref1", strme(&d1.fwr4_dna_ref1));
                                fwriteln!(log, "{} ref1", strme(&d2.fwr4_dna_ref1));
                                fwriteln!(log, "{} seq1", strme(&d1.fwr4_dna1));
                                fwriteln!(log, "{} seq2", strme(&d2.fwr4_dna1));
                                for i in 0..d1.fwr4_dna1.len() {
                                    if d1.fwr4_dna_ref1[i] != d2.fwr4_dna_ref1[i] {
                                        if d1.fwr4_dna1[i] == d1.fwr4_dna_ref1[i]
                                            && d2.fwr4_dna1[i] == d2.fwr4_dna_ref1[i] {
                                            supp += 1;
                                        } else if d1.fwr4_dna1[i] == d1.fwr4_dna_ref1[i]
                                            && d2.fwr4_dna1[i] == d1.fwr4_dna_ref1[i] {
                                            sup1 += 1;
                                        } else if d1.fwr4_dna1[i] == d2.fwr4_dna_ref1[i]
                                            && d2.fwr4_dna1[i] == d2.fwr4_dna_ref1[i] {
                                            sup2 += 1;
                                        }
                                    }
                                }
                            }

                            // Process leader

                            if d1.leader_dna_ref2 != d2.leader_dna_ref2 {
                                if d1.leader_dna_ref2.len() != d2.leader_dna_ref2.len() ||
                                   d1.leader_dna2.len() != d1.leader_dna_ref2.len() {
                                    fwriteln!(log, "\nleader: unequal lengths");
                                } else {
                                    fwriteln!(log, "\nleader");
                                    for i in 0..d1.leader_dna2.len() {
                                        if d1.leader_dna_ref2[i] == d2.leader_dna_ref2[i] {
                                            fwrite!(log, " ");
                                        } else {
                                            fwrite!(log, "*");
                                            refdiffs += 1;
                                        }
                                    }
                                    fwriteln!(log, "");
                                    fwriteln!(log, "{} ref2", strme(&d1.leader_dna_ref2));
                                    fwriteln!(log, "{} ref2", strme(&d2.leader_dna_ref2));
                                    fwriteln!(log, "{} seq1", strme(&d1.leader_dna2));
                                    fwriteln!(log, "{} seq2", strme(&d2.leader_dna2));
                                    for i in 0..d1.leader_dna2.len() {
                                        if d1.leader_dna_ref2[i] != d2.leader_dna_ref2[i] {
                                            if d1.leader_dna2[i] == d1.leader_dna_ref2[i]
                                                && d2.leader_dna2[i] == d2.leader_dna_ref2[i] {
                                                supp += 1;
                                            } else if d1.leader_dna2[i] == d1.leader_dna_ref2[i]
                                                && d2.leader_dna2[i] == d1.leader_dna_ref2[i] {
                                                sup1 += 1;
                                            } else if d1.leader_dna2[i] == d2.leader_dna_ref2[i]
                                                && d2.leader_dna2[i] == d2.leader_dna_ref2[i] {
                                                sup2 += 1;
                                            }
                                        }
                                    }
                                }
                            }

                            // Process FWR1.

                            if d1.fwr1_dna_ref2 != d2.fwr1_dna_ref2 {
                                fwriteln!(log, "\nFWR1");
                                for i in 0..d1.fwr1_dna2.len() {
                                    if d1.fwr1_dna_ref2[i] == d2.fwr1_dna_ref2[i] {
                                        fwrite!(log, " ");
                                    } else {
                                        fwrite!(log, "*");
                                        refdiffs += 1;
                                    }
                                }
                                fwriteln!(log, "");
                                fwriteln!(log, "{} ref2", strme(&d1.fwr1_dna_ref2));
                                fwriteln!(log, "{} ref2", strme(&d2.fwr1_dna_ref2));
                                fwriteln!(log, "{} seq1", strme(&d1.fwr1_dna2));
                                fwriteln!(log, "{} seq2", strme(&d2.fwr1_dna2));
                                for i in 0..d1.fwr1_dna2.len() {
                                    if d1.fwr1_dna_ref2[i] != d2.fwr1_dna_ref2[i] {
                                        if d1.fwr1_dna2[i] == d1.fwr1_dna_ref2[i]
                                            && d2.fwr1_dna2[i] == d2.fwr1_dna_ref2[i] {
                                            supp += 1;
                                        } else if d1.fwr1_dna2[i] == d1.fwr1_dna_ref2[i]
                                            && d2.fwr1_dna2[i] == d1.fwr1_dna_ref2[i] {
                                            sup1 += 1;
                                        } else if d1.fwr1_dna2[i] == d2.fwr1_dna_ref2[i]
                                            && d2.fwr1_dna2[i] == d2.fwr1_dna_ref2[i] {
                                            sup2 += 1;
                                        }
                                    }
                                }
                            }

                            // Process FWR2.

                            if d1.fwr2_dna_ref2 != d2.fwr2_dna_ref2 {
                                if d1.fwr2_dna_ref2.len() != d2.fwr2_dna_ref2.len() ||
                                   d1.fwr2_dna2.len() != d1.fwr2_dna_ref2.len() {
                                    fwriteln!(log, "\nFWR2: unequal lengths");
                                } else {
                                    fwriteln!(log, "\nFWR2");
                                    for i in 0..d1.fwr2_dna2.len() {
                                        if d1.fwr2_dna_ref2[i] == d2.fwr2_dna_ref2[i] {
                                            fwrite!(log, " ");
                                        } else {
                                            fwrite!(log, "*");
                                            refdiffs += 1;
                                        }
                                    }
                                    fwriteln!(log, "");
                                    fwriteln!(log, "{} ref2", strme(&d1.fwr2_dna_ref2));
                                    fwriteln!(log, "{} ref2", strme(&d2.fwr2_dna_ref2));
                                    fwriteln!(log, "{} seq1", strme(&d1.fwr2_dna2));
                                    fwriteln!(log, "{} seq2", strme(&d2.fwr2_dna2));
                                    for i in 0..d1.fwr2_dna2.len() {
                                        if d1.fwr2_dna_ref2[i] != d2.fwr2_dna_ref2[i] {
                                            if d1.fwr2_dna2[i] == d1.fwr2_dna_ref2[i]
                                                && d2.fwr2_dna2[i] == d2.fwr2_dna_ref2[i] {
                                                supp += 1;
                                            } else if d1.fwr2_dna2[i] == d1.fwr2_dna_ref2[i]
                                                && d2.fwr2_dna2[i] == d1.fwr2_dna_ref2[i] {
                                                sup1 += 1;
                                            } else if d1.fwr2_dna2[i] == d2.fwr2_dna_ref2[i]
                                                && d2.fwr2_dna2[i] == d2.fwr2_dna_ref2[i] {
                                                sup2 += 1;
                                            }
                                        }
                                    }
                                }
                            }

                            // Process FWR3.

                            if d1.fwr3_dna_ref2 != d2.fwr3_dna_ref2 {
                                if d1.fwr3_dna_ref2.len() != d2.fwr3_dna_ref2.len() ||
                                   d1.fwr3_dna2.len() != d1.fwr3_dna_ref2.len() {
                                    fwriteln!(log, "\nFWR3: unequal lengths");
                                } else {
                                    fwriteln!(log, "\nFWR3");
                                    for i in 0..d1.fwr3_dna2.len() {
                                        if d1.fwr3_dna_ref2[i] == d2.fwr3_dna_ref2[i] {
                                            fwrite!(log, " ");
                                        } else {
                                            fwrite!(log, "*");
                                            refdiffs += 1;
                                        }
                                    }
                                    fwriteln!(log, "");
                                    fwriteln!(log, "{} ref2", strme(&d1.fwr3_dna_ref2));
                                    fwriteln!(log, "{} ref2", strme(&d2.fwr3_dna_ref2));
                                    fwriteln!(log, "{} seq1", strme(&d1.fwr3_dna2));
                                    fwriteln!(log, "{} seq2", strme(&d2.fwr3_dna2));
                                    for i in 0..d1.fwr3_dna2.len() {
                                        if d1.fwr3_dna_ref2[i] != d2.fwr3_dna_ref2[i] {
                                            if d1.fwr3_dna2[i] == d1.fwr3_dna_ref2[i]
                                                && d2.fwr3_dna2[i] == d2.fwr3_dna_ref2[i] {
                                                supp += 1;
                                            } else if d1.fwr3_dna2[i] == d1.fwr3_dna_ref2[i]
                                                && d2.fwr3_dna2[i] == d1.fwr3_dna_ref2[i] {
                                                sup1 += 1;
                                            } else if d1.fwr3_dna2[i] == d2.fwr3_dna_ref2[i]
                                                && d2.fwr3_dna2[i] == d2.fwr3_dna_ref2[i] {
                                                sup2 += 1;
                                            }
                                        }
                                    }
                                }
                            }

                            // Process FWR4.

                            if d1.fwr4_dna_ref2 != d2.fwr4_dna_ref2 {
                                fwriteln!(log, "\nFWR4");
                                for i in 0..d1.fwr4_dna2.len() {
                                    if d1.fwr4_dna_ref2[i] == d2.fwr4_dna_ref2[i] {
                                        fwrite!(log, " ");
                                    } else {
                                        fwrite!(log, "*");
                                        refdiffs += 1;
                                    }
                                }
                                fwriteln!(log, "");
                                fwriteln!(log, "{} ref2", strme(&d1.fwr4_dna_ref2));
                                fwriteln!(log, "{} ref2", strme(&d2.fwr4_dna_ref2));
                                fwriteln!(log, "{} seq1", strme(&d1.fwr4_dna2));
                                fwriteln!(log, "{} seq2", strme(&d2.fwr4_dna2));
                                for i in 0..d1.fwr4_dna1.len() {
                                    if d1.fwr4_dna_ref1[i] != d2.fwr4_dna_ref1[i] {
                                        if d1.fwr4_dna1[i] == d1.fwr4_dna_ref1[i]
                                            && d2.fwr4_dna1[i] == d2.fwr4_dna_ref1[i] {
                                            supp += 1;
                                        } else if d1.fwr4_dna1[i] == d1.fwr4_dna_ref1[i]
                                            && d2.fwr4_dna1[i] == d1.fwr4_dna_ref1[i] {
                                            sup1 += 1;
                                        } else if d1.fwr4_dna1[i] == d2.fwr4_dna_ref1[i]
                                            && d2.fwr4_dna1[i] == d2.fwr4_dna_ref1[i] {
                                            sup2 += 1;
                                        }
                                    }
                                }
                                for i in 0..d1.fwr4_dna2.len() {
                                    if d1.fwr4_dna_ref2[i] != d2.fwr4_dna_ref2[i] {
                                        if d1.fwr4_dna2[i] == d1.fwr4_dna_ref2[i]
                                            && d2.fwr4_dna2[i] == d2.fwr4_dna_ref2[i] {
                                            supp += 1;
                                        } else if d1.fwr4_dna2[i] == d1.fwr4_dna_ref2[i]
                                            && d2.fwr4_dna2[i] == d1.fwr4_dna_ref2[i] {
                                            sup1 += 1;
                                        } else if d1.fwr4_dna2[i] == d2.fwr4_dna_ref2[i]
                                            && d2.fwr4_dna2[i] == d2.fwr4_dna_ref2[i] {
                                            sup2 += 1;
                                        }
                                    }
                                }
                            }

                            // Process CDR1.

                            if d1.cdr1_dna_ref2 != d2.cdr1_dna_ref2 {
                                if d1.cdr1_dna_ref2.len() != d2.cdr1_dna_ref2.len()
                                    || d1.cdr1_dna2.len() != d1.cdr1_dna_ref2.len() {
                                    fwriteln!(log, "\nCDR1: unequal lengths");
                                } else {
                                    fwriteln!(log, "\nCDR1");
                                    for i in 0..d1.cdr1_dna2.len() {
                                        if d1.cdr1_dna_ref2[i] == d2.cdr1_dna_ref2[i] {
                                            fwrite!(log, " ");
                                        } else {
                                            fwrite!(log, "*");
                                            refdiffs += 1;
                                        }
                                    }
                                    fwriteln!(log, "");
                                    fwriteln!(log, "{} ref2", strme(&d1.cdr1_dna_ref2));
                                    fwriteln!(log, "{} ref2", strme(&d2.cdr1_dna_ref2));
                                    fwriteln!(log, "{} seq1", strme(&d1.cdr1_dna2));
                                    fwriteln!(log, "{} seq2", strme(&d2.cdr1_dna2));
                                    for i in 0..d1.cdr1_dna2.len() {
                                        if d1.cdr1_dna_ref2[i] != d2.cdr1_dna_ref2[i] {
                                            if d1.cdr1_dna2[i] == d1.cdr1_dna_ref2[i]
                                                && d2.cdr1_dna2[i] == d2.cdr1_dna_ref2[i] {
                                                supp += 1;
                                            } else if d1.cdr1_dna2[i] == d1.cdr1_dna_ref2[i]
                                                && d2.cdr1_dna2[i] == d1.cdr1_dna_ref2[i] {
                                                sup1 += 1;
                                            } else if d1.cdr1_dna2[i] == d2.cdr1_dna_ref2[i]
                                                && d2.cdr1_dna2[i] == d2.cdr1_dna_ref2[i] {
                                                sup2 += 1;
                                            }
                                        }
                                    }
                                }
                            }

                            // Process CDR2.

                            if d1.cdr2_dna_ref2 != d2.cdr2_dna_ref2 {
                                if d1.cdr2_dna_ref2.len() != d2.cdr2_dna_ref2.len()
                                    || d1.cdr2_dna2.len() != d1.cdr2_dna_ref2.len() {
                                    fwriteln!(log, "\nCDR2: unequal lengths");
                                } else {
                                    fwriteln!(log, "\nCDR2");
                                    for i in 0..d1.cdr2_dna2.len() {
                                        if d1.cdr2_dna_ref2[i] == d2.cdr2_dna_ref2[i] {
                                            fwrite!(log, " ");
                                        } else {
                                            fwrite!(log, "*");
                                            refdiffs += 1;
                                        }
                                    }
                                    fwriteln!(log, "");
                                    fwriteln!(log, "{} ref2", strme(&d1.cdr2_dna_ref2));
                                    fwriteln!(log, "{} ref2", strme(&d2.cdr2_dna_ref2));
                                    fwriteln!(log, "{} seq1", strme(&d1.cdr2_dna2));
                                    fwriteln!(log, "{} seq2", strme(&d2.cdr2_dna2));
                                    for i in 0..d1.cdr2_dna2.len() {
                                        if d1.cdr2_dna_ref2[i] != d2.cdr2_dna_ref2[i] {
                                            if d1.cdr2_dna2[i] == d1.cdr2_dna_ref2[i]
                                                && d2.cdr2_dna2[i] == d2.cdr2_dna_ref2[i] {
                                                supp += 1;
                                            } else if d1.cdr2_dna2[i] == d1.cdr2_dna_ref2[i]
                                                && d2.cdr2_dna2[i] == d1.cdr2_dna_ref2[i] {
                                                sup1 += 1;
                                            } else if d1.cdr2_dna2[i] == d2.cdr2_dna_ref2[i]
                                                && d2.cdr2_dna2[i] == d2.cdr2_dna_ref2[i] {
                                                sup2 += 1;
                                            }
                                        }
                                    }
                                }
                            }

                            // Process CDR3.

                            if d1.cdr3_dna2.len() != d2.cdr3_dna2.len() {
                                fwriteln!(log, "\nCDR3: unequal lengths");
                            } else {
                                fwriteln!(log, "\nCDR3");
                                for i in 0..d1.cdr3_dna2.len() {
                                    if d1.cdr3_dna2[i] == d2.cdr3_dna2[i] {
                                        fwrite!(log, " ");
                                    } else {
                                        fwrite!(log, "*");
                                        refdiffs += 1;
                                    }
                                }
                                fwriteln!(log, "");
                                fwriteln!(log, "{} seq1", strme(&d1.cdr3_dna2));
                                fwriteln!(log, "{} seq2", strme(&d2.cdr3_dna2));
                            }

                            // Summarize. 

                            fwriteln!(log, "\nright = {supp}");
                            fwriteln!(log, "wrong1 = {sup1}");
                            fwriteln!(log, "wrong2 = {sup2}");
                            fwriteln!(log, "validated = {validated}");

                            let _refdiffs = refdiffs;

                            if !opt_same {
                                /*
                                for i in 0..n {
                                    if ref1[i] == ref2[i] {
                                        print!(" ");
                                    } else {
                                        print!("*");
                                        refdiffs += 1;
                                    }
                                }
                                let mut log = Vec::<u8>::new();
                                fwriteln!(log, "\n{} ref1", strme(&ref1));
                                fwriteln!(log, "{} ref2", strme(&ref2));
                                fwriteln!(log, "{} seq1", strme(&seq1));
                                fwriteln!(log, "{} seq2", strme(&seq2));
                                let mut supp = 0;
                                let mut sup1 = 0;
                                let mut sup2 = 0;
                                let mut other = 0;
                                for i in 0..n {
                                    if ref1[i] != ref2[i] {
                                        if seq1[i] == ref1[i] && seq2[i] == ref2[i] {
                                            supp += 1;
                                        } else if seq1[i] == ref1[i] && seq2[i] == ref1[i] {
                                            sup1 += 1;
                                        } else if seq1[i] == ref2[i] && seq2[i] == ref2[i] {
                                            sup2 += 1;
                                        } else {
                                            other += 1;
                                        }
                                    }
                                }
                                fwriteln!(log, "right = {supp}");
                                fwriteln!(log, "wrong1 = {sup1}");
                                fwriteln!(log, "wrong2 = {sup2}");
                                fwriteln!(log, "dunno = {other}");
                                fwrite!(log, "summary: ");
                                if refdiffs == 0 {
                                    fwriteln!(log, "no reference differences");
                                } else if supp == 1 && sup1 == 0 && sup2 == 0 {
                                    fwriteln!(log, "right == 1 and wrongs = 0");
                                } else if supp == 2 && sup1 == 0 && sup2 == 0 {
                                    fwriteln!(log, "right == 2 and wrongs = 0");
                                } else if supp == 3 && sup1 == 0 && sup2 == 0 {
                                    fwriteln!(log, "right == 3 and wrongs = 0");
                                } else if supp >= 4 && sup1 == 0 && sup2 == 0 {
                                    fwriteln!(log, "right >= 4 and wrongs = 0");
                                } else if supp == 1 && (sup1 > 0 || sup2 > 0) {
                                    fwriteln!(log, "right = 1 and at least one wrong > 0");
                                } else if supp == 2 && (sup1 > 0 || sup2 > 0) {
                                    fwriteln!(log, "right = 2 and at least one wrong > 0");
                                } else if supp == 3 && (sup1 > 0 || sup2 > 0) {
                                    fwriteln!(log, "right = 3 and at least one wrong > 0");
                                } else if supp >= 4 && (sup1 > 0 || sup2 > 0) {
                                    fwriteln!(log, "right >= 4 and at least one wrong > 0");
                                } else if supp == 0 && ((sup1 > 0) ^ (sup2 > 0)) {
                                    fwriteln!(log, "right = 0 and exactly one wrong > 0");
                                } else {
                                    fwriteln!(log, "other");
                                }
                                */
                            }
                        }
                    } else {
                        res.2[d + 1][ident].3 += 1;
                        res.2[0][ident].3 += 1;
                    }
                }
            }
        }
    });
    if show {
        for i in 0..bounds.len() {
            print!("{}", strme(&bounds[i].3));
        }
    }


    // Sum.

    let mut res = vec![vec![(0, 0, 0, 0); 11]; 5];
    for pass in 0..5 {
        for i in 0..bounds.len() {
            for j in 0..=10 {
                res[pass][j].0 += bounds[i].2[pass][j].0;
                res[pass][j].1 += bounds[i].2[pass][j].1;
                res[pass][j].2 += bounds[i].2[pass][j].2;
                res[pass][j].3 += bounds[i].2[pass][j].3;
            }
        }
    }

    // Print results.

    println!(
        "\nConsider two cells from the same donor that have the same CDR3H length,\n\
            and different heavy chain V genes."
    );
    println!("\nColumn 1: percent identity rounded down to nearest ten percent");
    println!("Column > 1: probability that light chain gene names are the same");
    let mut logs = Vec::<String>::new();
    for xpass in 1..=2 {
        let mut log = String::new();
        let mut rows = Vec::<Vec<String>>::new();
        let row = vec![
            "CDR3H-AA".to_string(),
            "log10(cell pairs)".to_string(),
            "all".to_string(),
            "d1".to_string(),
            "d2".to_string(),
            "d3".to_string(),
            "d4".to_string(),
        ];
        rows.push(row);
        for j in 0..=10 {
            let row = vec!["\\hline".to_string(); 7];
            rows.push(row);
            let mut row = vec![format!("{}%", 10 * j)];
            let n = if xpass == 1 {
                res[0][j].2 + res[0][j].3
            } else {
                res[0][j].0 + res[0][j].1
            };
            row.push(format!("{:.1}", (n as f64).log10()));
            for pass in 0..5 {
                if xpass == 1 {
                    let n = res[pass][j].2 + res[pass][j].3;
                    let nznz = 100.0 * res[pass][j].2 as f64 / n as f64;
                    row.push(format!("{nznz:.1}%"));
                } else {
                    let n = res[pass][j].0 + res[pass][j].1;
                    let nznz = 100.0 * res[pass][j].0 as f64 / n as f64;
                    row.push(format!("{nznz:.1}%"));
                }
            }
            rows.push(row);
        }
        print_tabular_vbox(&mut log, &rows, 0, &b"l|r|r|r|r|r|r".to_vec(), false, false);
        logs.push(log);
    }
    let mut logr = vec![Vec::<String>::new(); 2];
    for xpass in 0..2 {
        let r = logs[xpass].split('\n').map(str::to_owned).collect();
        logr[xpass] = r;
    }
    print!("\n both cells have dref > 0");
    print!("                               ");
    println!("both cells have dref = 0");
    let r = hcat(&logr[0], &logr[1], 3);
    for i in 0..r.len() {
        println!("{}", r[i]);
    }
}
