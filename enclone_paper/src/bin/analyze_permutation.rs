// Analyze results of permutation analysis.
//
// Run with one argument, which is a file containing lines like:
//
// memory:4.8%,4.6%,4.6%,4.6%,4.6%,4.6%,4.6%,4.7%,4.6%,4.5%,4.3%
// naive:4.6%,4.6%,4.6%,4.6%,4.6%,4.6%,4.6%,4.6%,4.6%,4.3%,3.5%

use io_utils::*;
use pretty_trace::*;
use std::io::BufRead;
use string_utils::*;

fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = std::env::args().collect();
    let f = open_for_read![&args[1]];
    let mut memory = vec![Vec::<f64>::new(); 11];
    let mut naive = vec![Vec::<f64>::new(); 11];
    for line in f.lines() {
        let mut s = line.unwrap();
        if s.starts_with("memory:") || s.starts_with("naive:") {
            let is_memory = s.starts_with("memory:");
            s = s.replace("%", "");
            let vals = s.after(":").split(',').collect::<Vec<&str>>();
            for (i, v) in vals.iter().enumerate() {
                if is_memory {
                    memory[i].push(v.force_f64());
                } else {
                    naive[i].push(v.force_f64());
                }
            }
        }
    }
    let mut mean_memory = vec![0.0; 11];
    let mut mean_naive = vec![0.0; 11];
    for i in 0..=10 {
        assert_eq!(memory[i].len(), 1000);
        memory[i].sort_by(|a, b| a.partial_cmp(b).unwrap());
        naive[i].sort_by(|a, b| a.partial_cmp(b).unwrap());
        for j in 0..memory[i].len() {
            mean_memory[i] += memory[i][j];
        }
        mean_memory[i] /= memory[i].len() as f64;
        for j in 0..naive[i].len() {
            mean_naive[i] += naive[i][j];
        }
        mean_naive[i] /= naive[i].len() as f64;
    }
    println!("stat,0,10,20,30,40,50,60,70,80,90,100");
    print!("memory_min");
    for i in 0..=10 {
        print!(",{:.1}", memory[i][0]);
    }
    println!("");
    print!("memory_mean");
    for i in 0..=10 {
        print!(",{:.3}", mean_memory[i]);
    }
    println!("");
    print!("memory_max");
    for i in 0..=10 {
        print!(",{:.1}", memory[i].last().unwrap());
    }
    println!("");
    print!("naive_min");
    for i in 0..=10 {
        print!(",{:.1}", naive[i][0]);
    }
    println!("");
    print!("naive_mean");
    for i in 0..=10 {
        print!(",{:.3}", mean_naive[i]);
    }
    println!("");
    print!("naive_max");
    for i in 0..=10 {
        print!(",{:.1}", naive[i].last().unwrap());
    }
    println!("");
}
