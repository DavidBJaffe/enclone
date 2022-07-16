// Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
//
// See also private_light_chain_analysis.rs.
//
// Analyze light chains.  Supply a single file of data, with one line per cell, and various fields.
//
// enclone BCR=@test BUILT_IN CHAINS_EXACT=2 CHAINS=2 NOPRINT POUT=stdout PCELL ECHOC
//         PCOLS=datasets_cell,donors_cell,v_name1,v_name2,dref,cdr3_aa1,clonotype_ncells,
//         const1,hcomp,jun_ins,jun_mat,jun_sub,d1_name1
//         > per_cell_stuff
//
// public_light_chain_analysis per_cell_stuff
//
// We allow per_cell_stuff to be replaced by a comma-separated list, in which case each element
// in the list is treated as having a single donor, unique to that list.
//
// Optional arguments:
// - FLOW: compute using flow classification of naive/memory rather than dref
// - FLOW_UNSWITCHED: same as flow but require memory to be sorted as unswitched
// - FLOW_SWITCHED: same as flow but require memory to be sorted as switched
// - NAIVE: compute just some stats about naive cells
// - NO_PARALOGS: do not make paralogs equivalent
// - REVERSE: reverse role of heavy and light
// - SOLO: use just one cell in a clonotype
// - SHOW: show information about each cell pair, and just for the 100% identity, equal
//         light chain case
// - MANY: work with many donors.
// - SWITCH_AS_DREF: if class switching occurs and dref = 0, change to dref = 1.
// - PERMUTE=n: permute the order of light chains n times.
// - FLIP_MEMORY=p: randomly flip the memory/naive state for p% of cells.

use enclone_core::test_def::test_donor_id;
use enclone_paper::public::*;
use io_utils::*;
use pretty_trace::PrettyTrace;
use rand::RngCore;
use rand::seq::SliceRandom;
use rand_chacha;
use rand_chacha::rand_core::SeedableRng;
use rayon::prelude::*;
use std::collections::HashMap;
use std::env;
use std::io::BufRead;
use std::mem::swap;
use string_utils::{add_commas, TextUtils};
use tables::*;
use vector_utils::{bin_position, erase_if, unique_sort};

const NAIVE: [usize; 40] = [
    1279049, 1279053, 1279057, 1279061, 1279065, 1279069, 1279073, 1279077, 1287144, 1287145,
    1287146, 1287147, 1287152, 1287153, 1287154, 1287155, 1287160, 1287161, 1287162, 1287163,
    1287168, 1287169, 1287170, 1287171, 1287176, 1287177, 1287178, 1287179, 1287184, 1287185,
    1287186, 1287187, 1287192, 1287193, 1287194, 1287195, 1287200, 1287201, 1287202, 1287203,
];
const PLASMABLAST: [usize; 6] = [1279052, 1279060, 1279068, 1279072, 1279076, 1279080];
const SWITCHED: [usize; 24] = [
    1279051, 1279055, 1279059, 1279063, 1279067, 1279071, 1279075, 1279079, 1287150, 1287151,
    1287158, 1287159, 1287166, 1287167, 1287174, 1287175, 1287182, 1287183, 1287190, 1287191,
    1287198, 1287199, 1287206, 1287207,
];
const UNSWITCHED: [usize; 24] = [
    1279050, 1279054, 1279058, 1279062, 1279066, 1279070, 1279074, 1279078, 1287148, 1287149,
    1287156, 1287157, 1287164, 1287165, 1287172, 1287173, 1287180, 1287181, 1287188, 1287189,
    1287196, 1287197, 1287204, 1287205,
];

fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();
    let mut opt_flow = false;
    let mut opt_naive = false;
    let mut opt_no_paralogs = false;
    let mut opt_reverse = false;
    let mut opt_solo = false;
    let mut opt_show = false;
    let mut opt_many = false;
    let mut flow_unswitched = false;
    let mut flow_switched = false;
    let mut switch_as_dref = false;
    let mut permute = 0;
    let mut flip_memory = 0.0;
    for i in 2..args.len() {
        if args[i] == "FLOW" {
            opt_flow = true;
        } else if args[i] == "FLOW_UNSWITCHED" {
            flow_unswitched = true;
        } else if args[i] == "FLOW_SWITCHED" {
            flow_switched = true;
        } else if args[i] == "NAIVE" {
            opt_naive = true;
        } else if args[i] == "NO_PARALOGS" {
            opt_no_paralogs = true;
        } else if args[i] == "REVERSE" {
            opt_reverse = true;
        } else if args[i] == "SOLO" {
            opt_solo = true;
        } else if args[i] == "SHOW" {
            opt_show = true;
        } else if args[i] == "MANY" {
            opt_many = true;
        } else if args[i].starts_with("PERMUTE=") {
            permute = args[i].after("PERMUTE=").force_usize();
        } else if args[i].starts_with("FLIP_MEMORY=") {
            flip_memory = args[i].after("FLIP_MEMORY=").force_f64();
        } else if args[i] == "SWITCH_AS_DREF" {
            switch_as_dref = true;
        }
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    #[derive(Clone, Eq, Ord, PartialEq, PartialOrd)]
    struct CellData {
        v_name1: String,
        cdr3_aa1_len: usize,
        cdr3_aa1: Vec<u8>,
        donor: String,
        v_name2: String,
        dref: usize,
        jun_mat: usize,
        jun_sub: usize,
        hcomp: usize,
        jun_ins: usize,
        dataset: String,
        d1_name1: String,
        j_name2: String,
        barcode: String,
        v_name2_orig: String,
        const1: String,
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Load data.

    let inputs = args[1].split(',').collect::<Vec<&str>>();
    let mut data = Vec::new();
    let mut clonotype = Vec::<usize>::new();
    for i in 0..inputs.len() {
        let f = open_for_read![&inputs[i]];
        let mut first = true;
        let mut tof = HashMap::<String, usize>::new();
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
                let mut donor = fields[tof["donors_cell"]].to_string();
                if inputs.len() > 1 {
                    donor = format!("d{}", i + 1);
                }
                if !opt_reverse {
                    data.push(CellData {
                        v_name1: fields[tof["v_name1"]].to_string(),
                        cdr3_aa1_len: fields[tof["cdr3_aa1"]].len(),
                        cdr3_aa1: fields[tof["cdr3_aa1"]].to_string().as_bytes().to_vec(),
                        donor: donor,
                        v_name2: fields[tof["v_name2"]].to_string(),
                        dref: fields[tof["dref"]].force_usize(),
                        jun_mat: fields[tof["jun_mat"]].force_usize(),
                        jun_sub: fields[tof["jun_sub"]].force_usize(),
                        hcomp: fields[tof["hcomp"]].force_usize(),
                        jun_ins: fields[tof["jun_ins"]].force_usize(),
                        dataset: fields[tof["datasets_cell"]].to_string(),
                        d1_name1: fields[tof["d1_name1"]].to_string(),
                        j_name2: fields[tof["j_name2"]].to_string(),
                        barcode: fields[tof["barcode"]].to_string(),
                        v_name2_orig: fields[tof["v_name2"]].to_string(),
                        const1: fields[tof["const1"]].to_string(),
                    });
                } else {
                    data.push(CellData {
                        v_name1: fields[tof["v_name2"]].to_string(),
                        cdr3_aa1_len: fields[tof["cdr3_aa2"]].len(),
                        cdr3_aa1: fields[tof["cdr3_aa2"]].to_string().as_bytes().to_vec(),
                        donor: donor,
                        v_name2: fields[tof["v_name1"]].to_string(),
                        dref: fields[tof["dref"]].force_usize(),
                        jun_mat: fields[tof["jun_mat"]].force_usize(),
                        jun_sub: fields[tof["jun_sub"]].force_usize(),
                        hcomp: fields[tof["hcomp"]].force_usize(),
                        jun_ins: fields[tof["jun_ins"]].force_usize(),
                        dataset: fields[tof["datasets_cell"]].to_string(),
                        d1_name1: fields[tof["d1_name1"]].to_string(),
                        j_name2: fields[tof["j_name1"]].to_string(),
                        barcode: fields[tof["barcode"]].to_string(),
                        v_name2_orig: fields[tof["v_name1"]].to_string(),
                        const1: fields[tof["const1"]].to_string(),
                    });
                }
                clonotype.push(fields[tof["group_id"]].force_usize());
            }
        }
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // In solo case, reduce to one cell per clonotype.

    if opt_solo {
        let mut to_delete = vec![false; data.len()];
        for i in 0..clonotype.len() {
            if i > 0 && clonotype[i] == clonotype[i - 1] {
                to_delete[i] = true;
            }
        }
        erase_if(&mut data, &to_delete);
    }

    // Implement SWITCH_AS_DREF.

    if switch_as_dref {
        for i in 0..data.len() {
            if data[i].dref == 0 {
                if data[i].const1.starts_with("IGHA")
                    || data[i].const1.starts_with("IGHE")
                    || data[i].const1.starts_with("IGHG")
                {
                    data[i].dref = 1;
                }
            }
        }
    }

    // Replace paralogs.

    if !opt_no_paralogs && !opt_reverse {
        for i in 0..data.len() {
            data[i].v_name2 = data[i].v_name2.replace("D", "");
        }
    }

    // Flip memory and naive cells.

    let mut randme = rand_chacha::ChaCha8Rng::seed_from_u64(123456789);
    if flip_memory > 0.0 {
        for i in 0..data.len() {
            let r = randme.next_u64();
            if r % 1_000_000 < ((flip_memory/100.0) * 1_000_000.0).round() as u64 {
                if data[i].dref > 0 {
                    data[i].dref = 0;
                } else {
                    data[i].dref = 1;
                }
            }
        }
    }

    // Permute.

    let mut reps = 1;
    if permute > 0 {
        reps = permute;
    }
    for _ in 0..reps {
        if permute > 0 {
            let mut data2 = data.clone();
            data2.shuffle(&mut randme);
            for i in 0..data.len() {
                data[i].v_name2 = data2[i].v_name2.clone();
                data[i].v_name2_orig = data2[i].v_name2_orig.clone();
                data[i].j_name2 = data2[i].j_name2.clone();
            }
        }

        // Sort.

        data.sort();

        // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
    
        // Compute number of public naive and memory cells.
    
        if permute == 0 {
            for pass in 1..=2 {
                let mut p = 0;
                let mut i = 0;
                while i < data.len() {
                    let mut j = i + 1;
                    while j < data.len() {
                        if data[j].v_name1 != data[i].v_name1 
                            || data[j].cdr3_aa1 != data[i].cdr3_aa1 {
                            break;
                        }
                        j += 1;
                    }
                    let mut dx = Vec::new();
                    for k in i..j {
                        if (pass == 1 && data[k].dref == 0) || (pass == 2 && data[k].dref > 0) {
                            dx.push(data[k].clone());
                        }
                    }
                    let mut donors = Vec::<String>::new();
                    for k in 0..dx.len() {
                        donors.push(dx[k].donor.clone());
                    }
                    unique_sort(&mut donors);
                    if donors.len() > 1 {
                        p += dx.len();
                    }
                    i = j;
                }
                if !opt_show {
                    if pass == 1 {
                        println!("\n{} public naive cells", p);
                    } else {
                        println!("{} public memory cells", p);
                    }
                }
            }
        }
    
        // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
    
        // Show the CDRH3 length distribution for naive cells.
    
        if !opt_show && !opt_many && permute == 0 {
            println!("\nCDRH3 length distribution for naive cells");
            let mut bins = vec![0; 100];
            let mut total = 0;
            let mut total_len = 0;
            for k in 0..data.len() {
                let dref = data[k].dref;
                if dref == 0 {
                    let len = data[k].cdr3_aa1_len;
                    bins[len / 5] += 1;
                    total += 1;
                    total_len += len;
                }
            }
            for i in 0..bins.len() {
                if bins[i] > 0 {
                    println!(
                        "{}-{} ==> {:.1}%",
                        5 * i,
                        5 * (i + 1),
                        100.0 * bins[i] as f64 / total as f64
                    );
                }
            }
            println!("mean CDRH3 length = {:.1}", total_len as f64 / total as f64);
        }
    
        // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
    
        // Compute naive fraction for each of the sort classes.
    
        let mut is_naive = vec![false; data.len()];
        let mut is_memory = vec![false; data.len()];
        let mut is_unswitched_flow = vec![false; data.len()];
        let mut is_switched_flow = vec![false; data.len()];
        if (opt_naive || opt_flow || flow_unswitched || flow_switched) && !opt_many 
            && permute == 0 {
            let mut all = Vec::<usize>::new();
            all.append(&mut NAIVE.to_vec());
            all.append(&mut UNSWITCHED.to_vec());
            all.append(&mut SWITCHED.to_vec());
            all.append(&mut PLASMABLAST.to_vec());
            all.sort();
            let mut naive = vec![(0, 0, 0, 0); 5];
            let mut unswitched = vec![(0, 0, 0, 0); 5];
            let mut switched = vec![(0, 0, 0, 0); 5];
            let mut plasmablast = vec![(0, 0, 0, 0); 5];
            let mut memory_subtotal = vec![(0, 0, 0, 0); 5];
            let mut unswitched_naive = vec![(0, 0, 0, 0); 5]; // THESE ARE MIXED!
            let mut switched_naive = vec![(0, 0, 0, 0); 5]; // THESE ARE MIXED!
            let mut total = vec![(0, 0, 0, 0); 5];
            let mut cells = vec![(0, 0); all.len()];
            for pass in 1..=2 {
                for i in 0..data.len() {
                    let c1 = &data[i].const1;
                    let is_switched =
                        c1.starts_with("IGHA") || c1.starts_with("IGHE") || c1.starts_with("IGHG");
                    let dref = data[i].dref;
                    let dataset = data[i].dataset.force_usize();
                    if pass == 1 && dataset.to_string().starts_with("128") {
                        continue;
                    }
                    if pass == 2 && dataset.to_string().starts_with("127") {
                        continue;
                    }
                    let donor_id = test_donor_id(dataset);
                    let p = bin_position(&all, &dataset) as usize;
                    cells[p].1 += 1;
                    if dref == 0 {
                        cells[p].0 += 1;
                    }
                    total[0].1 += 1;
                    total[donor_id].1 += 1;
                    if dref == 0 {
                        total[0].0 += 1;
                        total[donor_id].0 += 1;
                    }
                    if dref == 0 && is_switched {
                        total[0].2 += 1;
                        total[donor_id].2 += 1;
                    }
                    if dref > 0 && is_switched {
                        total[0].3 += 1;
                        total[donor_id].3 += 1;
                    }
                    if NAIVE.contains(&dataset) {
                        naive[0].1 += 1;
                        naive[donor_id].1 += 1;
                        is_naive[i] = true;
                        if dref == 0 {
                            naive[0].0 += 1;
                            naive[donor_id].0 += 1;
                        }
                        if dref == 0 && is_switched {
                            naive[0].2 += 1;
                            naive[donor_id].2 += 1;
                        }
                        if dref > 0 && is_switched {
                            naive[0].3 += 1;
                            naive[donor_id].3 += 1;
                        }
                    } else if UNSWITCHED.contains(&dataset) {
                        if pass == 1 || donor_id == 1 {
                            unswitched[0].1 += 1;
                            unswitched[donor_id].1 += 1;
                            memory_subtotal[0].1 += 1;
                            memory_subtotal[donor_id].1 += 1;
                            is_memory[i] = true;
                            is_unswitched_flow[i] = true;
                        } else {
                            unswitched_naive[0].1 += 1;
                            unswitched_naive[donor_id].1 += 1;
                        }
                        if dref == 0 {
                            if pass == 1 || donor_id == 1 {
                                unswitched[0].0 += 1;
                                unswitched[donor_id].0 += 1;
                                memory_subtotal[0].0 += 1;
                                memory_subtotal[donor_id].0 += 1;
                            } else {
                                unswitched_naive[0].0 += 1;
                                unswitched_naive[donor_id].0 += 1;
                            }
                        }
                        if dref == 0 && is_switched {
                            if pass == 1 || donor_id == 1 {
                                unswitched[0].2 += 1;
                                unswitched[donor_id].2 += 1;
                                memory_subtotal[0].2 += 1;
                                memory_subtotal[donor_id].2 += 1;
                            } else {
                                unswitched_naive[0].2 += 1;
                                unswitched_naive[donor_id].2 += 1;
                            }
                        }
                        if dref > 0 && is_switched {
                            if pass == 1 || donor_id == 1 {
                                unswitched[0].3 += 1;
                                unswitched[donor_id].3 += 1;
                                memory_subtotal[0].3 += 1;
                                memory_subtotal[donor_id].3 += 1;
                            } else {
                                unswitched_naive[0].3 += 1;
                                unswitched_naive[donor_id].3 += 1;
                            }
                        }
                    } else if SWITCHED.contains(&dataset) {
                        if pass == 1 {
                            switched[0].1 += 1;
                            switched[donor_id].1 += 1;
                            memory_subtotal[0].1 += 1;
                            memory_subtotal[donor_id].1 += 1;
                            is_memory[i] = true;
                            is_switched_flow[i] = true;
                        } else {
                            switched_naive[0].1 += 1;
                            switched_naive[donor_id].1 += 1;
                        }
                        if dref == 0 {
                            if pass == 1 {
                                switched[0].0 += 1;
                                switched[donor_id].0 += 1;
                                memory_subtotal[0].0 += 1;
                                memory_subtotal[donor_id].0 += 1;
                            } else {
                                switched_naive[0].0 += 1;
                                switched_naive[donor_id].0 += 1;
                            }
                        }
                        if dref == 0 && is_switched {
                            if pass == 1 {
                                switched[0].2 += 1;
                                switched[donor_id].2 += 1;
                                memory_subtotal[0].2 += 1;
                                memory_subtotal[donor_id].2 += 1;
                            } else {
                                switched_naive[0].2 += 1;
                                switched_naive[donor_id].2 += 1;
                            }
                        }
                        if dref > 0 && is_switched {
                            if pass == 1 {
                                switched[0].3 += 1;
                                switched[donor_id].3 += 1;
                                memory_subtotal[0].3 += 1;
                                memory_subtotal[donor_id].3 += 1;
                            } else {
                                switched_naive[0].3 += 1;
                                switched_naive[donor_id].3 += 1;
                            }
                        }
                    } else if PLASMABLAST.contains(&dataset) {
                        plasmablast[0].1 += 1;
                        plasmablast[donor_id].1 += 1;
                        memory_subtotal[0].1 += 1;
                        memory_subtotal[donor_id].1 += 1;
                        is_memory[i] = true;
                        if dref == 0 {
                            plasmablast[0].0 += 1;
                            plasmablast[donor_id].0 += 1;
                            memory_subtotal[0].0 += 1;
                            memory_subtotal[donor_id].0 += 1;
                        }
                        if dref == 0 && is_switched {
                            plasmablast[0].2 += 1;
                            plasmablast[donor_id].2 += 1;
                            memory_subtotal[0].2 += 1;
                            memory_subtotal[donor_id].2 += 1;
                        }
                        if dref > 0 && is_switched {
                            plasmablast[0].3 += 1;
                            plasmablast[donor_id].3 += 1;
                            memory_subtotal[0].3 += 1;
                            memory_subtotal[donor_id].3 += 1;
                        }
                    } else {
                        panic!("unclassified dataset");
                    }
                }
            }
    
            // Print tables.
    
            if opt_naive {
                let counts = [
                    &naive,
                    &unswitched,
                    &switched,
                    &plasmablast,
                    &memory_subtotal,
                    &unswitched_naive,
                    &switched_naive,
                    &total,
                ];
                let names = [
                    "naive",
                    "unswitched",
                    "switched",
                    "plasmablast",
                    "memory_subtotal",
                    "unswitched_naive",
                    "switched_naive",
                    "total",
                ];
                let row1 = vec![
                    "class".to_string(),
                    "all".to_string(),
                    "d1".to_string(),
                    "d2".to_string(),
                    "d3".to_string(),
                    "d4".to_string(),
                ];
                println!("\nall cells");
                let mut rows = vec![row1.clone()];
                for i in 0..counts.len() {
                    rows.push(vec!["\\hline".to_string(); 6]);
                    let mut row = vec![names[i].to_string()];
                    for j in 0..5 {
                        if counts[i][j].1 > 0 {
                            row.push(format!("{}", add_commas(counts[i][j].1)));
                        } else {
                            row.push(String::new());
                        }
                    }
                    rows.push(row);
                }
                let mut log = String::new();
                print_tabular_vbox(&mut log, &rows, 0, &b"l|r|r|r|r|r".to_vec(), false, false);
                println!("{}", log);
    
                println!("naive switched cells");
                let mut rows = vec![row1.clone()];
                for i in 0..counts.len() {
                    rows.push(vec!["\\hline".to_string(); 6]);
                    let mut row = vec![names[i].to_string()];
                    for j in 0..5 {
                        if counts[i][j].1 > 0 {
                            row.push(format!("{}", add_commas(counts[i][j].2),));
                        } else {
                            row.push(String::new());
                        }
                    }
                    rows.push(row);
                }
                let mut log = String::new();
                print_tabular_vbox(&mut log, &rows, 0, &b"l|r|r|r|r|r".to_vec(), false, false);
                println!("{}", log);
    
                println!("memory switched cells");
                let mut rows = vec![row1.clone()];
                for i in 0..counts.len() {
                    rows.push(vec!["\\hline".to_string(); 6]);
                    let mut row = vec![names[i].to_string()];
                    for j in 0..5 {
                        if counts[i][j].1 > 0 {
                            row.push(format!("{}", add_commas(counts[i][j].3),));
                        } else {
                            row.push(String::new());
                        }
                    }
                    rows.push(row);
                }
                let mut log = String::new();
                print_tabular_vbox(&mut log, &rows, 0, &b"l|r|r|r|r|r".to_vec(), false, false);
                println!("{}", log);
    
                println!("naive cell fractions");
                let mut rows = vec![row1.clone()];
                for i in 0..counts.len() {
                    rows.push(vec!["\\hline".to_string(); 6]);
                    let mut row = vec![names[i].to_string()];
                    for j in 0..5 {
                        if counts[i][j].1 > 0 {
                            row.push(format!(
                                "{:.1}%",
                                100.0 * counts[i][j].0 as f64 / counts[i][j].1 as f64
                            ));
                        } else {
                            row.push(String::new());
                        }
                    }
                    rows.push(row);
                }
                let mut log = String::new();
                print_tabular_vbox(&mut log, &rows, 0, &b"l|r|r|r|r|r".to_vec(), false, false);
                println!("{}", log);
    
                std::process::exit(0);
            }
        }
    
        // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
    
        // Compute DD fraction.
    
        let mut naive = (0, 0);
        let mut memory = (0, 0);
        if !opt_many && permute == 0 {
            for i in 0..data.len() {
                let dref = data[i].dref;
                let d = &data[i].d1_name1;
                if dref == 0 {
                    naive.1 += 1;
                    if d.contains(":") {
                        naive.0 += 1;
                    }
                } else {
                    memory.1 += 1;
                    if d.contains(":") {
                        memory.0 += 1;
                    }
                }
            }
            if !opt_show {
                println!(
                    "\nDD in naive cells = {} = {:.2}%",
                    naive.0,
                    100.0 * naive.0 as f64 / naive.1 as f64
                );
                println!(
                    "DD in memory cells = {} = {:.2}%",
                    memory.0,
                    100.0 * memory.0 as f64 / memory.1 as f64
                );
            }
        }
    
        // Compute DD stuff.
    
        let mut i = 0;
        let mut total = vec![0; 2];
        let mut dd_memory = 0;
        let mut dd_naive = 0;
        if !opt_many {
            println!("");
            while i < data.len() {
                let mut j = i + 1;
                while j < data.len() {
                    if data[j].v_name1 != data[i].v_name1 || data[j].cdr3_aa1 != data[i].cdr3_aa1 {
                        break;
                    }
                    j += 1;
                }
                let mut donors = Vec::<String>::new();
                for k in i..j {
                    donors.push(data[k].donor.clone());
                }
                unique_sort(&mut donors);
                if donors.len() > 1 {
                    for k in i..j {
                        let dref = data[k].dref;
                        if dref == 0 {
                            total[1] += 1;
                            if data[k].d1_name1.contains(":") {
                                dd_naive += 1;
                            }
                        } else {
                            total[0] += 1;
                            if data[k].d1_name1.contains(":") {
                                dd_memory += 1;
                            }
                        }
                    }
                }
                i = j;
            }
            if !opt_show && permute == 0 {
                println!(
                    "DD public memory cells = {} = {:.1}%",
                    dd_memory,
                    100.0 * dd_memory as f64 / total[0] as f64,
                );
                println!(
                    "DD public naive cells = {} = {:.1}%",
                    dd_naive,
                    100.0 * dd_naive as f64 / total[1] as f64,
                );
            }
        }
    
        // Define groups based on equal heavy chain gene names and CDRH3 length.
        // Plus placeholder for results, see next.
        // We track the number of cell pairs (first vector), and the number of cells 
        // (second vector).
    
        let mut bounds = Vec::new();
        let mut i = 0;
        while i < data.len() {
            // let j = next_diff12_9(&data, i as i32) as usize;
            let mut j = i + 1;
            while j < data.len() {
                if data[j].v_name1 != data[i].v_name1 
                    || data[j].cdr3_aa1_len != data[i].cdr3_aa1_len {
                    break;
                }
                j += 1;
            }
            bounds.push((
                i,
                j,
                vec![vec![(0, 0, 0, 0); 11]; 7],
                vec![vec![(0, 0); 11]; 7],
            ));
            i = j;
        }
    
        // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
    
        // Results = for each percent identity, rounded down:
        // 1. count for equal light chain gene names and dref1 = 0 and dref2 = 0
        // 2. count for unequal light chain gene names and dref1 = 0 and dref2 = 0
        // 3. count for equal light chain gene names and dref1 > 0 and dref2 > 0
        // 4. count for unequal light chain gene names and dref1 > 0 and dref2 > 0.
        //
        // Make one pass for all donors, and one pass each for each pair of donors.
    
        bounds.par_iter_mut().for_each(|res| {
            let i = res.0;
            let j = res.1;
            let mut cells = vec![vec![(Vec::new(), Vec::new()); 11]; 7];
            for k1 in i..j {
                for k2 in k1 + 1..j {
                    // Require different donors.
    
                    if data[k1].donor == data[k2].donor {
                        continue;
                    }
                    let (mut d1, mut d2) = (data[k1].donor.clone(), data[k2].donor.clone());
                    if d1 > d2 {
                        swap(&mut d1, &mut d2);
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
    
                    // Go through passes.
    
                    for pass in 0..7 {
                        // Require specific donors.
    
                        if pass == 0 {
                            if opt_show && ident == 10 && eq_light {
                                let mut comment = String::new();
                                if data[k1].j_name2 != data[k2].j_name2
                                    || data[k1].v_name2_orig != data[k2].v_name2_orig
                                {
                                    comment = " ***".to_string();
                                }
                                println!(
                                    "{} {} {} {} ==> {} {} {} {} {}",
                                    data[k1].dataset,
                                    data[k1].barcode,
                                    data[k1].v_name2_orig,
                                    data[k1].j_name2,
                                    data[k2].dataset,
                                    data[k2].barcode,
                                    data[k2].v_name2_orig,
                                    data[k2].j_name2,
                                    comment,
                                );
                            }
                        } else if pass == 1 {
                            if d1 != "d1" || d2 != "d2" {
                                continue;
                            }
                        } else if pass == 2 {
                            if d1 != "d1" || d2 != "d3" {
                                continue;
                            }
                        } else if pass == 3 {
                            if d1 != "d1" || d2 != "d4" {
                                continue;
                            }
                        } else if pass == 4 {
                            if d1 != "d2" || d2 != "d3" {
                                continue;
                            }
                        } else if pass == 5 {
                            if d1 != "d2" || d2 != "d4" {
                                continue;
                            }
                        } else {
                            if d1 != "d3" || d2 != "d4" {
                                continue;
                            }
                        }
    
                        // Add to results.
    
                        let mut naive = dref1 == 0 && dref2 == 0;
                        let mut memory = dref1 > 0 && dref2 > 0;
                        if opt_flow {
                            naive = is_naive[k1] && is_naive[k2];
                            memory = is_memory[k1] && is_memory[k2];
                        }
                        if flow_unswitched {
                            naive = is_naive[k1] && is_naive[k2];
                            memory = is_unswitched_flow[k1] && is_unswitched_flow[k2];
                        }
                        if flow_switched {
                            naive = is_naive[k1] && is_naive[k2];
                            memory = is_switched_flow[k1] && is_switched_flow[k2];
                        }
                        if naive {
                            cells[pass][ident].0.push(k1);
                            cells[pass][ident].0.push(k2);
                            if eq_light {
                                res.2[pass][ident].0 += 1;
                            } else {
                                res.2[pass][ident].1 += 1;
                            }
                        } else if memory {
                            cells[pass][ident].1.push(k1);
                            cells[pass][ident].1.push(k2);
                            if eq_light {
                                res.2[pass][ident].2 += 1;
                            } else {
                                res.2[pass][ident].3 += 1;
                            }
                        }
                    }
                }
            }
            for pass in 0..7 {
                for ident in 0..=10 {
                    unique_sort(&mut cells[pass][ident].0);
                    unique_sort(&mut cells[pass][ident].1);
                    res.3[pass][ident].0 = cells[pass][ident].0.len();
                    res.3[pass][ident].1 = cells[pass][ident].1.len();
                }
            }
        });
        if opt_show {
            std::process::exit(0);
        }
    
        // Sum.
    
        let mut res = vec![vec![(0, 0, 0, 0); 11]; 7];
        let mut res_cell = vec![vec![(0, 0); 11]; 7];
        for pass in 0..7 {
            for i in 0..bounds.len() {
                for j in 0..=10 {
                    res[pass][j].0 += bounds[i].2[pass][j].0;
                    res[pass][j].1 += bounds[i].2[pass][j].1;
                    res[pass][j].2 += bounds[i].2[pass][j].2;
                    res[pass][j].3 += bounds[i].2[pass][j].3;
                    res_cell[pass][j].0 += bounds[i].3[pass][j].0;
                    res_cell[pass][j].1 += bounds[i].3[pass][j].1;
                }
            }
        }
    
        // Print results.

        if permute > 0 {
            for xpass in 1..=2 {
                let pass = 0;
                if xpass == 1 {
                    print!("memory:");
                    for j in 0..=10 {
                        if j > 0 { 
                            print!(",");
                        }
                        let n = res[pass][j].2 + res[pass][j].3;
                        let nznz = 100.0 * res[pass][j].2 as f64 / n as f64;
                        print!("{nznz:.1}%");
                    }
                    println!("");
                } else {
                    print!("naive:");
                    for j in 0..=10 {
                        if j > 0 { 
                            print!(",");
                        }
                        let n = res[pass][j].0 + res[pass][j].1;
                        let nznz = 100.0 * res[pass][j].0 as f64 / n as f64;
                        print!("{nznz:.1}%");
                    }
                    println!("");
                }
            }
        } else {
            public_print_results(&res, &res_cell, opt_many);
        }
    }
}
