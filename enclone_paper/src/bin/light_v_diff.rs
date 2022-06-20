// Copyright (c) 2022 10x Genomics, Inc. All rights reserved.
//
// Analyze the distance between light chain V genes.

use edit_distance::*;
use pretty_trace::*;
use string_utils::strme;
use string_utils::TextUtils;
use vdj_ann_ref::human_ref;

fn main() {
    PrettyTrace::new().on();

    // Get the human reference used by enclone.

    let mut refnames = Vec::<String>::new();
    let mut refs = Vec::<Vec<u8>>::new();
    let mut utr = Vec::<bool>::new();
    let href = human_ref();
    {
        for (i, line) in href.lines().enumerate() {
            if i % 2 == 0 {
                let n = line.between("|", " ").to_string();
                utr.push(line.contains("5'UTR"));
                refnames.push(n);
            } else {
                refs.push(line.as_bytes().to_vec());
            }
        }
    }

    // Pairwise compare.

    for i1 in 0..refnames.len() {
        if !refnames[i1].starts_with("IGKV") && !refnames[i1].starts_with("IGLV") {
            continue;
        }
        for i2 in i1 + 1..refnames.len() {
            if utr[i1] || utr[i2] {
                continue;
            }
            if !refnames[i2].starts_with("IGKV") && !refnames[i2].starts_with("IGLV") {
                continue;
            }
            if refnames[i1] == refnames[i2] {
                continue;
            }
            let seq1 = &refs[i1];
            let seq2 = &refs[i2];
            /*
            let n1 = refnames[i1].replace("D", "");
            let n2 = refnames[i2].replace("D", "");
            if n1 != n2 {
                continue;
            }
            */
            let dist = edit_distance(&strme(&seq1), &strme(&seq2));
            if dist <= 20 {
                println!("{} ==> {} ==> {}", refnames[i1], refnames[i2], dist);
            }
        }
    }
}
