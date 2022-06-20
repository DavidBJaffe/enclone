// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

// Code split from public_light_chain_analysis to reduce file size.

use enclone_core::hcat;
use tables::*;

pub fn public_print_results(
    res: &Vec<Vec<(usize, usize, usize, usize)>>,
    res_cell: &Vec<Vec<(usize, usize)>>,
    opt_many: bool,
) {
    println!(
        "\nConsider two cells from different donors that have the same heavy chain gene name \
        and CDRH3 length."
    );
    println!("\nColumn 1: percent identity rounded down to nearest ten percent");
    println!("Column > 1: probability that light chain gene names are the same");
    let mut logs = Vec::<String>::new();
    for xpass in 1..=2 {
        let mut log = String::new();
        let mut rows = Vec::<Vec<String>>::new();
        let mut row = vec![
            "CDRH3-AA".to_string(),
            "log10(cell pairs)".to_string(),
            "any".to_string(),
        ];
        if !opt_many {
            row.append(&mut vec![
                "d1,d2".to_string(),
                "d1,d3".to_string(),
                "d1,d4".to_string(),
                "d2,d3".to_string(),
                "d2,d4".to_string(),
                "d3,d4".to_string(),
            ]);
        }
        rows.push(row);
        for j in 0..=10 {
            let mut cols = 9;
            if opt_many {
                cols = 3;
            }
            let row = vec!["\\hline".to_string(); cols];
            rows.push(row);
            let mut row = vec![format!("{}%", 10 * j)];
            let n = if xpass == 1 {
                res[0][j].2 + res[0][j].3
            } else {
                res[0][j].0 + res[0][j].1
            };
            row.push(format!("{:.1}", (n as f64).log10()));
            for pass in 0..7 {
                if opt_many && pass == 1 {
                    break;
                }
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
        let mut just = b"l|r|r|r|r|r|r|r|r".to_vec();
        if opt_many {
            just = b"l|r|r".to_vec();
        }
        print_tabular_vbox(&mut log, &rows, 0, &just, false, false);
        logs.push(log);
    }
    let mut logr = vec![Vec::<String>::new(); 2];
    for xpass in 0..2 {
        let r = logs[xpass].split('\n').map(str::to_owned).collect();
        logr[xpass] = r;
    }
    print!("\n both cells have dref > 0");
    if !opt_many {
        print!("                                                 ");
    } else {
        print!("             ");
    }
    println!("both cells have dref = 0");
    let r = hcat(&logr[0], &logr[1], 3);
    for i in 0..r.len() {
        println!("{}", r[i]);
    }

    // Print tables showing cell pair counts.

    if opt_many {
        std::process::exit(0);
    }
    println!("Cell pair counts:");
    let mut logs = Vec::<String>::new();
    for xpass in 1..=2 {
        let mut log = String::new();
        let mut rows = Vec::<Vec<String>>::new();
        let row = vec![
            "CDRH3-AA".to_string(),
            "any".to_string(),
            "d1,d2".to_string(),
            "d1,d3".to_string(),
            "d1,d4".to_string(),
            "d2,d3".to_string(),
            "d2,d4".to_string(),
            "d3,d4".to_string(),
        ];
        rows.push(row);
        for j in 0..=10 {
            let row = vec!["\\hline".to_string(); 8];
            rows.push(row);
            let mut row = vec![format!("{}%", 10 * j)];
            for pass in 0..7 {
                if xpass == 1 {
                    let n = res[pass][j].2 + res[pass][j].3;
                    row.push(format!("{:.1}", (n as f64).log10()));
                } else {
                    let n = res[pass][j].0 + res[pass][j].1;
                    row.push(format!("{:.1}", (n as f64).log10()));
                }
            }
            rows.push(row);
        }
        print_tabular_vbox(
            &mut log,
            &rows,
            0,
            &b"l|r|r|r|r|r|r|r".to_vec(),
            false,
            false,
        );
        logs.push(log);
    }
    let mut logr = vec![Vec::<String>::new(); 2];
    for xpass in 0..2 {
        let r = logs[xpass].split('\n').map(str::to_owned).collect();
        logr[xpass] = r;
    }
    print!("\n both cells have dref > 0");
    print!("                            ");
    println!("both cells have dref = 0");
    let r = hcat(&logr[0], &logr[1], 3);
    for i in 0..r.len() {
        println!("{}", r[i]);
    }

    // Print tables showing cell counts.

    println!("Cell counts:");
    let mut logs = Vec::<String>::new();
    for xpass in 1..=2 {
        let mut log = String::new();
        let mut rows = Vec::<Vec<String>>::new();
        let row = vec![
            "CDRH3-AA".to_string(),
            "any".to_string(),
            "d1,d2".to_string(),
            "d1,d3".to_string(),
            "d1,d4".to_string(),
            "d2,d3".to_string(),
            "d2,d4".to_string(),
            "d3,d4".to_string(),
        ];
        rows.push(row);
        for j in 0..=10 {
            let row = vec!["\\hline".to_string(); 8];
            rows.push(row);
            let mut row = vec![format!("{}%", 10 * j)];
            for pass in 0..7 {
                if xpass == 1 {
                    let n = res_cell[pass][j].1;
                    row.push(format!("{:.1}", (n as f64).log10()));
                } else {
                    let n = res_cell[pass][j].0;
                    row.push(format!("{:.1}", (n as f64).log10()));
                }
            }
            rows.push(row);
        }
        print_tabular_vbox(
            &mut log,
            &rows,
            0,
            &b"l|r|r|r|r|r|r|r".to_vec(),
            false,
            false,
        );
        logs.push(log);
    }
    let mut logr = vec![Vec::<String>::new(); 2];
    for xpass in 0..2 {
        let r = logs[xpass].split('\n').map(str::to_owned).collect();
        logr[xpass] = r;
    }
    print!("\n both cells have dref > 0");
    print!("                            ");
    println!("both cells have dref = 0");
    let r = hcat(&logr[0], &logr[1], 3);
    for i in 0..r.len() {
        println!("{}", r[i]);
    }
}
