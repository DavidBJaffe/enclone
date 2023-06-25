// Copyright (c) 2022 10x Genomics, Inc. All rights reserved.

pub mod align_to_vdj_ref;
pub mod blacklist;
pub mod opt_d;
pub mod packing;
pub mod pdb;
pub mod print_tools;

use string_utils::*;

pub fn expand_integer_ranges(x: &str) -> String {
    let mut tokens = Vec::<String>::new();
    let mut token = String::new();
    for c in x.chars() {
        if c == ',' || c == ':' || c == ';' {
            if !token.is_empty() {
                tokens.push(token.clone());
                token.clear();
            }
            tokens.push(c.to_string());
        } else {
            token.push(c);
        }
    }
    if !token.is_empty() {
        tokens.push(token);
    }
    let mut tokens2 = Vec::<String>::new();
    for i in 0..tokens.len() {
        if tokens[i].contains('-')
            && tokens[i].before("-").parse::<usize>().is_ok()
            && tokens[i].after("-").parse::<usize>().is_ok()
        {
            let n1 = tokens[i].before("-").force_usize();
            let n2 = tokens[i].after("-").force_usize();
            if n1 <= n2 {
                for n in n1..=n2 {
                    if n > n1 {
                        tokens2.push(",".to_string());
                    }
                    tokens2.push(format!("{}", n));
                }
                continue;
            }
        }
        tokens2.push(tokens[i].clone());
    }
    let mut y = String::new();
    for i in 0..tokens2.len() {
        y += &tokens2[i];
    }
    y
}
