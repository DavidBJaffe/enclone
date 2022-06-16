
// usage:
// map_donors infile outfile new-donor-name
//
// Read CSV infile, and change all the donors_cell entries to new-donor-name.  Write to outfile.

use io_utils::*;
use itertools::Itertools;
use pretty_trace::PrettyTrace;
use std::env;
use std::io::{BufRead, Write};

fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();
    let (infile, outfile) = (&args[1], &args[2]);
    let new_donor_name = &args[3];
    let f = open_for_read![&infile];
    let mut g = open_for_write_new![&outfile];
    let mut dcf = 1000000;
    let mut first = true;
    for line in f.lines() {
        let s = line.unwrap();
        if s.starts_with("#") {
            fwriteln!(g, "{}", s);
        } else {
            let mut fields: Vec<String> = s.split(',').map(str::to_owned).collect();
            if first {
                for j in 0..fields.len() {
                    if fields[j] == "donors_cell" {
                        dcf = j;
                    }
                }
                first = false;
            } else {
                fields[dcf] = new_donor_name.to_string();
            }
            fwriteln!(g, "{}", fields.iter().format(","));
        }
    }
}
