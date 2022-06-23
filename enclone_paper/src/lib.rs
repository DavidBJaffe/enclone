// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

pub mod dperfect;
pub mod public;

use string_utils::*;

pub fn unpack_aa_matrix_string(z: &str) -> Vec<Vec<f64>> {
    let mut m = Vec::<Vec<f64>>::new();
    for line in z.lines() {
        let mut s = line.to_string();
        let sb = s.replace(" ", "");
        if sb == "ACDEFGHIKLMNPQRSTVWY" {
            continue;
        }
        if s.len() > 2 && s.as_bytes()[0] >= b'A' {
            s = s[2..].to_string();
        }
        let fields = s.split(' ').collect::<Vec<&str>>();
        let mut row = Vec::<f64>::new();
        for i in 0..fields.len() {
            row.push(fields[i].force_f64());
        }
        m.push(row);
    }
    m
}
