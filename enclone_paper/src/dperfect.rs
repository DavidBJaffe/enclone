// Given a human heavy chain junction region, report the longest perfect match by a human D gene to
// the region, where the D gene can be placed within the region.

use std::cmp::max;

pub fn dperfect(x: &[u8], ds: &Vec<Vec<u8>>) -> usize {
    let mut best = 0;
    for i in 0..ds.len() {
        let d = &ds[i];
        if x.len() >= d.len() {
            for s in 0..=x.len() - d.len() {
                let mut m = 0;
                for j in 0..d.len() {
                    if x[s + j] == d[j] {
                        m += 1;
                        best = max(best, m);
                    } else {
                        m = 0;
                    }
                }
            }
        }
    }
    best
}
