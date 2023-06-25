// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Find the optimal D segment, the runner up, and the delta between the scores.  This uses
// the donor V and J segments that are assigned to the clonotype.  Note that the optimal D
// segment may be null.  This is obvious from looking at data.

use crate::align_to_vdj_ref::{align_to_vdj_ref, match_bit_score, zero_one};
use enclone_proto::types::DonorReferenceItem;
use std::cmp::min;
use vdj_ann::refx::RefData;

pub fn vflank(_seq: &[u8], vref: &[u8]) -> usize {
    let mut flank = 13;
    if flank > vref.len() {
        flank = vref.len();
    }
    flank
}

pub fn jflank(seq: &[u8], jref: &[u8]) -> usize {
    let flank = 13;
    if flank > jref.len() {
        return jref.len();
    }

    // Let start be the first position on the J gene where there is a perfect match of length at
    // least five to the contig.

    const MATCHLEN: usize = 5;
    let mut start = 0;
    for i in 0..=jref.len() - MATCHLEN {
        let mut matchlen = 0;
        for j in 0..MATCHLEN {
            if seq.len() >= jref.len() && seq[seq.len() - jref.len() + i + j] != jref[i + j] {
                break;
            }
            matchlen += 1;
        }
        if matchlen == MATCHLEN {
            start = i;
            break;
        }
    }

    // Add start to the flank, so long as that's possible.

    min(flank + start, jref.len())
}

pub fn evaluate_d(
    tig: &[u8],
    vref: &[u8],
    seq_start: usize,
    ds: &Vec<usize>,
    jref: &[u8],
    refdata: &RefData,
    jscore_match: i32,
    jscore_mismatch: i32,
    jscore_gap_open: i32,
    jscore_gap_extend: i32,
    jscore_bits_multiplier: f64,
    new: bool,
) -> (Vec<bio_edit::alignment::AlignmentOperation>, f64) {
    // Start to build reference concatenation.  First append the V segment.

    let mut concat = Vec::<u8>::new();
    let vstart = vref.len() - vflank(tig, vref);
    let vref = vref[vstart..vref.len()].to_vec();
    concat.append(&mut vref.clone());

    // Append the D segment or segments.

    let mut dref = Vec::<u8>::new();
    let mut d2ref = Vec::<u8>::new();
    let mut drefname = String::new();
    for j in 0..ds.len() {
        let d = ds[j];
        if j == 0 {
            dref = refdata.refs[d].to_ascii_vec();
        } else if j == 1 {
            d2ref = refdata.refs[d].to_ascii_vec();
        }
        if j > 0 {
            drefname += ":";
        }
        drefname += &mut refdata.name[d].clone();
    }
    concat.append(&mut dref.clone());
    concat.append(&mut d2ref.clone());

    // Append the J segment.

    let jend = jflank(tig, jref);

    // Align the V..J sequence on the contig to the reference concatenation.

    let mut seq_end;
    if tig.len() < jref.len() - jend {
        seq_end = tig.len();
    } else {
        seq_end = tig.len() - (jref.len() - jend);
    }
    if seq_end <= seq_start as usize {
        seq_end = tig.len(); // bug fix for problem found by customer, couldn't reproduce internally
    }
    let seq = tig[seq_start as usize..seq_end].to_vec();
    let jref = jref[0..jend].to_vec();
    concat.append(&mut jref.clone());
    let (ops, count) = align_to_vdj_ref(
        &seq,
        &Vec::new(),
        &vref,
        &dref,
        &d2ref,
        &jref,
        &Vec::new(),
        &Vec::new(),
        &drefname,
        true,
        jscore_match,
        jscore_mismatch,
        jscore_gap_open,
        jscore_gap_extend,
        jscore_bits_multiplier,
        new,
    );
    (ops, count)
}

pub fn opt_d(
    v_ref_id: usize,                       // ex.share[mid].v_ref_id
    jref: &[u8],                           // reference J segment
    tig: &Vec<u8>,                         // ex.share[mid].seq_del
    annv: &Vec<(i32, i32, i32, i32, i32)>, // ex.share[mid].annv
    cdr3_aa: &str,                         // ex.share[mid].cdr3_aa
    refdata: &RefData,
    dref: &Vec<DonorReferenceItem>,
    scores: &mut Vec<f64>,
    dsx: &mut Vec<Vec<usize>>,
    jscore_match: i32,
    jscore_mismatch: i32,
    jscore_gap_open: i32,
    jscore_gap_extend: i32,
    jscore_bits_multiplier: f64,
    v_alt: Option<usize>,
    new: bool,
) {
    let mut comp = 1000000.0;

    // Go through every D segment, or possibly every concatenation of D segments.

    let mut todo = Vec::<Vec<usize>>::new();
    todo.push(vec![]);
    for i in refdata.ds.iter() {
        todo.push(vec![*i]);
    }
    let mut ds = Vec::<Vec<usize>>::new();
    let mut counts = Vec::<f64>::new();
    let mut good_d = Vec::<usize>::new();
    let mut vref = refdata.refs[v_ref_id].to_ascii_vec();
    if v_alt.is_some() {
        vref = dref[v_alt.unwrap()].nt_sequence.clone();
    }
    let vstart = vref.len() - vflank(tig, &vref);
    let mut seq_start = vstart as isize;
    // probably not exactly right
    if annv.len() > 1 {
        let q1 = annv[0].0 + annv[0].1;
        let q2 = annv[1].0;

        seq_start += q1 as isize - q2 as isize;
    }
    // Workaround for very rare anomaly that would cause a subsequent assert.
    if seq_start as usize > tig.len() {
        seq_start = tig.len() as isize;
    }
    const MIN_BITS_FOR_D2: f64 = 14.0;
    for di in 0..todo.len() {
        let (ops, count) = evaluate_d(
            tig,
            &vref,
            seq_start as usize,
            &todo[di],
            &jref,
            refdata,
            jscore_match,
            jscore_mismatch,
            jscore_gap_open,
            jscore_gap_extend,
            jscore_bits_multiplier,
            new,
        );
        counts.push(count);
        if !todo[di].is_empty() {
            let drefx = refdata.refs[todo[di][0]].to_ascii_vec();
            let vstart = vref.len() - vflank(tig, &vref);
            let vref = vref[vstart..vref.len()].to_vec();
            let zos = zero_one(&ops, vref.len(), vref.len() + drefx.len());
            let bits = match_bit_score(&zos);
            if bits >= MIN_BITS_FOR_D2 {
                good_d.push(todo[di][0]);
            }
        }
        ds.push(todo[di].clone());
        if count > comp {
            comp = count;
        }
    }
    if cdr3_aa.len() >= 20 {
        todo.clear();
        for i1 in good_d.iter() {
            for i2 in good_d.iter() {
                todo.push(vec![*i1, *i2]);
            }
        }
        for di in 0..todo.len() {
            let (_ops, count) = evaluate_d(
                tig,
                &vref,
                seq_start as usize,
                &todo[di],
                &jref,
                refdata,
                jscore_match,
                jscore_mismatch,
                jscore_gap_open,
                jscore_gap_extend,
                jscore_bits_multiplier,
                new,
            );
            counts.push(count);
            ds.push(todo[di].clone());
            if count > comp {
                comp = count;
            }
        }
    }

    // Reverse sort sync (counts, ds).

    let mut counts_ds = Vec::new();
    for i in 0..counts.len() {
        counts_ds.push((counts[i], ds[i].clone()));
    }
    counts_ds.sort_by(|a, b| b.partial_cmp(a).unwrap()); // reverse sort
    counts.clear();
    ds.clear();
    for i in 0..counts_ds.len() {
        counts.push(counts_ds[i].0);
        ds.push(counts_ds[i].1.clone());
    }

    // Done.

    *scores = counts;
    *dsx = ds;
}