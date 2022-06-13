// Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
//
// See also public_light_chain_analysis.rs.
//
// Consider pairs of cells from the same donor, having the same CDRH3 lengths and same FWR1 
// lengths.  Compute supp, the support for the assertion that the two cells are from different
// clonotypes.  This is the number of positions for which the FWR1 references differ, and both
// cells agree with their reference there.
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
//
// SOLO: reduce to one cell per clonotype.

use enclone_core::hcat;
use io_utils::*;
use pretty_trace::PrettyTrace;
use rayon::prelude::*;
use std::cmp::max;
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
    let mut solo = false;
    for i in 2..args.len() {
        if args[i] == "SHOW" {
            show = true;
        } else if args[i] == "SOLO" {
            solo = true;
        } else {
            eprintln!("\nIllegal argument.\n");
            std::process::exit(1);
        }
    }

    #[derive(Eq, Ord, PartialEq, PartialOrd)]
    struct CellData {
        donor: String,
        cdr3_len1: usize,
        fwr1_dna_len1: usize,
        fwr2_dna_len1: usize,
        fwr3_dna_len1: usize,
        fwr4_dna_len1: usize,
        cdr3_aa1: Vec<u8>,
        v_name1: String,
        v_name2: String,
        dref: usize,
        dataset: String,
        barcode: String,
        j_name1: String,
        fwr1_dna_ref1: Vec<u8>,
        fwr1_dna1: Vec<u8>,
        fwr2_dna_ref1: Vec<u8>,
        fwr2_dna1: Vec<u8>,
        fwr3_dna_ref1: Vec<u8>,
        fwr3_dna1: Vec<u8>,
        fwr4_dna_ref1: Vec<u8>,
        fwr4_dna1: Vec<u8>,
    }

    let mut first = true;
    let mut tof = HashMap::<String, usize>::new();
    let mut data = Vec::<CellData>::new();
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
                fwr1_dna_len1: fields[tof["fwr1_dna1"]].to_string().len(),
                fwr2_dna_len1: fields[tof["fwr2_dna1"]].to_string().len(),
                fwr3_dna_len1: fields[tof["fwr3_dna1"]].to_string().len(),
                fwr4_dna_len1: fields[tof["fwr4_dna1"]].to_string().len(),
                cdr3_aa1: fields[tof["cdr3_aa1"]].to_string().as_bytes().to_vec(),
                v_name1: fields[tof["v_name1"]].to_string(),
                v_name2: fields[tof["v_name2"]].to_string(),
                dref: fields[tof["dref"]].force_usize(),
                dataset: fields[tof["datasets_cell"]].to_string(),
                barcode: fields[tof["barcode"]].to_string(),
                j_name1: fields[tof["j_name1"]].to_string(),
                fwr1_dna_ref1: fields[tof["fwr1_dna_ref1"]].to_string().as_bytes().to_vec(),
                fwr1_dna1: fields[tof["fwr1_dna1"]].to_string().as_bytes().to_vec(),
                fwr2_dna_ref1: fields[tof["fwr2_dna_ref1"]].to_string().as_bytes().to_vec(),
                fwr2_dna1: fields[tof["fwr2_dna1"]].to_string().as_bytes().to_vec(),
                fwr3_dna_ref1: fields[tof["fwr3_dna_ref1"]].to_string().as_bytes().to_vec(),
                fwr3_dna1: fields[tof["fwr3_dna1"]].to_string().as_bytes().to_vec(),
                fwr4_dna_ref1: fields[tof["fwr4_dna_ref1"]].to_string().as_bytes().to_vec(),
                fwr4_dna1: fields[tof["fwr4_dna1"]].to_string().as_bytes().to_vec(),
            });
            clonotype.push(fields[tof["group_id"]].force_usize());
        }
    }
    data.sort();

    // Require no indels in certain places.

    let mut to_delete = vec![false; data.len()];
    for i in 0..data.len() {
        if data[i].fwr1_dna1.len() != data[i].fwr1_dna_ref1.len() {
            to_delete[i] = true;
        }
        if data[i].fwr2_dna1.len() != data[i].fwr2_dna_ref1.len() {
            to_delete[i] = true;
        }
        if data[i].fwr3_dna1.len() != data[i].fwr3_dna_ref1.len() {
            to_delete[i] = true;
        }
        if data[i].fwr4_dna1.len() != data[i].fwr4_dna_ref1.len() {
            to_delete[i] = true;
        }
    }
    erase_if(&mut data, &to_delete);
    erase_if(&mut clonotype, &to_delete);

    // Replace paralogous light chain gene names.

    for i in 0..data.len() {
        data[i].v_name2 = data[i].v_name2.replace("D", "");
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

    // Define groups based on equal donor and CDR3H length.
    // Plus placeholder for results, see next.

    let mut bounds = Vec::new();
    let mut i = 0;
    while i < data.len() {
        let mut j = i + 1;
        while j < data.len() {
            if data[j].donor != data[i].donor || data[j].cdr3_len1 != data[i].cdr3_len1
                || data[j].fwr1_dna_len1 != data[i].fwr1_dna_len1
                || data[j].fwr2_dna_len1 != data[i].fwr2_dna_len1
                || data[j].fwr3_dna_len1 != data[i].fwr3_dna_len1
                || data[j].fwr4_dna_len1 != data[i].fwr4_dna_len1 {
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
                let d1 = &data[k1];
                let d2 = &data[k2];

                // Calculate framework support.

                let mut supp = 0;
                let mut sup1 = 0;
                let mut sup2 = 0;
                for i in 0..d1.fwr1_dna1.len() {
                    if d1.fwr1_dna_ref1[i] != d2.fwr1_dna_ref1[i] {
                        if d1.fwr1_dna1[i] == d1.fwr1_dna_ref1[i]
                            && d2.fwr1_dna1[i] == d2.fwr1_dna_ref1[i] {
                            supp += 1;
                        } else if d1.fwr1_dna1[i] == d1.fwr1_dna_ref1[i] 
                            && d2.fwr1_dna1[i] == d1.fwr1_dna_ref1[i] {
                            sup1 += 1;
                        } else if d1.fwr1_dna1[i] == d2.fwr1_dna_ref1[i] 
                            && d2.fwr1_dna1[i] == d2.fwr1_dna_ref1[i] {
                            sup2 += 1;
                        }
                    }
                }
                for i in 0..d1.fwr2_dna1.len() {
                    if d1.fwr2_dna_ref1[i] != d2.fwr2_dna_ref1[i] {
                        if d1.fwr2_dna1[i] == d1.fwr2_dna_ref1[i]
                            && d2.fwr2_dna1[i] == d2.fwr2_dna_ref1[i] {
                            supp += 1;
                        } else if d1.fwr2_dna1[i] == d1.fwr2_dna_ref1[i] 
                            && d2.fwr2_dna1[i] == d1.fwr2_dna_ref1[i] {
                            sup1 += 1;
                        } else if d1.fwr2_dna1[i] == d2.fwr2_dna_ref1[i] 
                            && d2.fwr2_dna1[i] == d2.fwr2_dna_ref1[i] {
                            sup2 += 1;
                        }
                    }
                }
                for i in 0..d1.fwr3_dna1.len() {
                    if d1.fwr3_dna_ref1[i] != d2.fwr3_dna_ref1[i] {
                        if d1.fwr3_dna1[i] == d1.fwr3_dna_ref1[i]
                            && d2.fwr3_dna1[i] == d2.fwr3_dna_ref1[i] {
                            supp += 1;
                        } else if d1.fwr3_dna1[i] == d1.fwr3_dna_ref1[i] 
                            && d2.fwr3_dna1[i] == d1.fwr3_dna_ref1[i] {
                            sup1 += 1;
                        } else if d1.fwr3_dna1[i] == d2.fwr3_dna_ref1[i] 
                            && d2.fwr3_dna1[i] == d2.fwr3_dna_ref1[i] {
                            sup2 += 1;
                        }
                    }
                }
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

                // Set floor on support.

                if supp < 3 || supp < 3 * max(sup1, sup2) {
                    continue;
                }

                // Compute stuff.

                let mut same = 0;
                for m in 0..d1.cdr3_aa1.len() {
                    if d1.cdr3_aa1[m] == d2.cdr3_aa1[m] {
                        same += 1;
                    }
                }
                let ident = 100.0 * same as f64 / data[k1].cdr3_aa1.len() as f64;
                let ident = ident.floor() as usize;
                let ident = ident / 10;
                let (dref1, dref2) = (d1.dref, d2.dref);
                let eq_light = d1.v_name2 == d2.v_name2;

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

                            // Start analysis of cell pair.

                            let mut refdiffs = 0;
                            fwriteln!(log, "\n{}", strme(&vec![b'='; 100]));
                            fwriteln!(log,
                                "\n{} {} {} {}",
                                data[k1].dataset, data[k1].barcode, 
                                data[k2].dataset, data[k2].barcode,
                            );

                            // Process FWR1.

                            fwriteln!(log, "\nFWR1");
                            for i in 0..d1.fwr1_dna1.len() {
                                if d1.fwr1_dna_ref1[i] == d2.fwr1_dna_ref1[i] {
                                    fwrite!(log, " ");
                                } else {
                                    fwrite!(log, "*");
                                    refdiffs += 1;
                                }
                            }
                            fwriteln!(log, "");
                            fwriteln!(log, "{} ref1", strme(&d1.fwr1_dna_ref1));
                            fwriteln!(log, "{} ref2", strme(&d2.fwr1_dna_ref1));
                            fwriteln!(log, "{} seq1", strme(&d1.fwr1_dna1));
                            fwriteln!(log, "{} seq2", strme(&d2.fwr1_dna1));

                            // Process FWR2.

                            fwriteln!(log, "\nFWR2");
                            for i in 0..d1.fwr2_dna1.len() {
                                if d1.fwr2_dna_ref1[i] == d2.fwr2_dna_ref1[i] {
                                    fwrite!(log, " ");
                                } else {
                                    fwrite!(log, "*");
                                    refdiffs += 1;
                                }
                            }
                            fwriteln!(log, "");
                            fwriteln!(log, "{} ref1", strme(&d1.fwr2_dna_ref1));
                            fwriteln!(log, "{} ref2", strme(&d2.fwr2_dna_ref1));
                            fwriteln!(log, "{} seq1", strme(&d1.fwr2_dna1));
                            fwriteln!(log, "{} seq2", strme(&d2.fwr2_dna1));

                            // Process FWR3.

                            fwriteln!(log, "\nFWR3");
                            for i in 0..d1.fwr3_dna1.len() {
                                if d1.fwr3_dna_ref1[i] == d2.fwr3_dna_ref1[i] {
                                    fwrite!(log, " ");
                                } else {
                                    fwrite!(log, "*");
                                    refdiffs += 1;
                                }
                            }
                            fwriteln!(log, "");
                            fwriteln!(log, "{} ref1", strme(&d1.fwr3_dna_ref1));
                            fwriteln!(log, "{} ref2", strme(&d2.fwr3_dna_ref1));
                            fwriteln!(log, "{} seq1", strme(&d1.fwr3_dna1));
                            fwriteln!(log, "{} seq2", strme(&d2.fwr3_dna1));

                            // Process FWR4.

                            fwriteln!(log, "\nFWR4");
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
                            fwriteln!(log, "{} ref2", strme(&d2.fwr4_dna_ref1));
                            fwriteln!(log, "{} seq1", strme(&d1.fwr4_dna1));
                            fwriteln!(log, "{} seq2", strme(&d2.fwr4_dna1));

                            // Summarize.

                            fwriteln!(log, "\nright = {supp}");
                            fwriteln!(log, "wrong1 = {sup1}");
                            fwriteln!(log, "wrong2 = {sup2}");
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
        std::process::exit(0);
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
