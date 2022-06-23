// Given a file of junction nucleotide sequences, determine the distribution of dperfection
// values, defined as dperfect/junction-region-length.  See dperfect.rs.
//
// usage: dperfection junction-region-file
//
// How to create input file from per_cell_stuff, e.g. for d1:
//
// cat per_cell_stuff | cut -d, -f2,7,16 | grep d1 | grep ",0," | cut -d, -f3 > d1_cdr3_dna1_naive

use enclone_paper::dperfect::dperfect;
use io_utils::*;
use pretty_trace::*;
use std::env;
use std::io::BufRead;
use vdj_ann_ref::human_ref;

fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();
    let f = open_for_read![&args[1]];

    // Fetch junction regions.

    let mut juns = Vec::<Vec<u8>>::new();
    for line in f.lines() {
        let s = line.unwrap();
        juns.push(s.as_bytes().to_vec());
    }

    // Fetch D gene reference sequences.

    let mut ds = Vec::<Vec<u8>>::new();
    let href = human_ref();
    let mut hlines = Vec::<String>::new();
    for line in href.lines() {
        hlines.push(line.to_string());
    }
    for i in (0..hlines.len()).step_by(2) {
        if hlines[i].contains("|D-REGION|IG|") {
            ds.push(hlines[i + 1].as_bytes().to_vec());
        }
    }

    // Compute dperfection.

    let mut dps = Vec::<f64>::new();
    for i in 0..juns.len() {
        let dperf = dperfect(&juns[i], &ds);
        dps.push(dperf as f64 / juns[i].len() as f64);
    }

    // Print histogram.

    let mut buckets = vec![0; 11];
    for i in 0..dps.len() {
        let b = (dps[i] * 10.0).floor() as usize;
        buckets[b] += 1;
    }
    println!("\ndperfection%\tdata%");
    for i in 0..=10 {
        let label = if i < 10 { format!("{}-", 10 * i) } else { format!("100") };
        println!("{label}\t{:.1}", 100.0 * buckets[i] as f64 / dps.len() as f64);
    }
    println!("");
}
