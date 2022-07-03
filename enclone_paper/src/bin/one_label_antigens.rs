// Label antigens in input CSV file.
//
// usage:
// one label_antigens <filename>
//
// These two files are derived from the two respective papers.
//
// v_name1: temporary handling!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use debruijn::dna_string::DnaString;
use fasta_tools::*;
use io_utils::*;
use pretty_trace::*;
use std::collections::HashMap;
use std::env;
use std::io::BufRead;
use string_utils::*;
use vdj_ann::annotate::*;
use vdj_ann::refx::*;
use vdj_ann_ref::human_ref;

fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();
    let stuff = &args[1];

    // Define data structure.

    #[derive(Clone, Eq, Ord, PartialEq, PartialOrd)]
    struct CellData {
        v_name1: String,
        v_name2: String,
        cdr3_aa1: Vec<u8>,
    }

    // Load data.

    let mut data = Vec::<CellData>::new();
    let mut headers = Vec::<String>::new();
    let mut seqs = Vec::<DnaString>::new();
    let mut vgenes = Vec::<String>::new();
    let mut cdrh3 = Vec::<Vec<u8>>::new();
    let mut antigens = Vec::<String>::new();
    {
        let href = human_ref();
        let mut refdata = RefData::new();
        let ext_ref = String::new();
        make_vdj_ref_data_core(&mut refdata, &href, &ext_ref, false, true, None);
        let fasta_file = include_str!["../../data/phad_clonal_structure_2022.fasta"];
        read_fasta_contents_into_vec_dna_string_plus_headers(&fasta_file, &mut seqs, &mut headers);
        for i in 0..seqs.len() {
            antigens.push(headers[i].between("A-1 ", " specific").to_string());
            let x = &seqs[i];
            let mut ann = Vec::<(i32, i32, i32, i32, i32)>::new();
            annotate_seq(&x, &refdata, &mut ann, true, false, true);
            let mut found_v = false;
            for j in 0..ann.len() {
                let t = ann[j].2 as usize;
                if refdata.is_v(t) {
                    vgenes.push(refdata.name[t].clone());
                    found_v = true;
                    break;
                }
            }
            if !found_v {
                eprintln!("\nFailed to find V gene for sequence {}.\n", i + 1);
                std::process::exit(1);
            }
            let mut cdr3x = Vec::<(usize, Vec<u8>, usize, usize)>::new();
            get_cdr3_using_ann(&x, &refdata, &ann, &mut cdr3x);
            if cdr3x.len() != 1 {
                eprintln!("failed to find unique CDR3\n");
                std::process::exit(1);
            }
            cdrh3.push(cdr3x[0].1.clone());
        }
        let f = open_for_read![&stuff];
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
                data.push(CellData {
                    // v_name1: fields[tof["v_name1"]].to_string(),
                    v_name1: "IGHV3-9".to_string(),
                    v_name2: fields[tof["v_name2"]].to_string(),
                    cdr3_aa1: fields[tof["cdr3_aa1"]].to_string().as_bytes().to_vec(),
                });
            }
        }
    }

    // Search the data for V gene and CDRH3-AA matches.

    for i in 0..data.len() {
        for j in 0..seqs.len() {
            if data[i].v_name1 == vgenes[j] {
                if data[i].cdr3_aa1.len() == cdrh3[j].len() {
                    let mut diffs = 0;
                    for k in 0..cdrh3[j].len() {
                        if data[i].cdr3_aa1[k] != cdrh3[j][k] {
                            diffs += 1;
                        }
                    }
                    if diffs <= 4 {
                        println!("{}/{}/{}, diffs = {diffs}, {} cell {}, seq {}", 
                            data[i].v_name1, data[i].v_name2, strme(&data[i].cdr3_aa1),
                            i + 1, j + 1, antigens[j]);
                    }
                }
            }
        }
    }
}
