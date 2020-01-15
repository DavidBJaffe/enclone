// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// This file provides functions to find alternate alleles and substitute them into references.

use vdj_ann::*;

use self::refx::*;
use crate::defs::*;
use debruijn::{dna_string::*, Mer};
use itertools::Itertools;
use perf_stats::*;
use stats_utils::*;
use std::cmp::*;
use std::time::Instant;
use vector_utils::*;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Output is {(donor, ref id, alt seq)}.

pub fn find_alleles(
    refdata: &RefData,
    ctl: &EncloneControl,
    exact_clonotypes: &Vec<ExactClonotype>,
) -> Vec<(usize, usize, DnaString)> {
    // Derive consensus sequences for alternate alleles of V segments.
    //
    // The priority of this algorithm is to reduce the likelihood of false positive joins.  It
    // has been tuned and optimized for this purpose, possibly at the expense of generating some
    // 'fake' alternate alleles.  However we do not actually know that this happens.  If we
    // wanted to test the true accuracy of allele finding, we would need whole-genome and VDJ
    // datasets from the same donors.
    //
    // This calculation has to be made separately for each donor, which means that in a certain
    // sense the algorithm is not blind to the truth data.  But it is the right thing to do.
    //
    // Alternate alleles might correspond to duplicated segments, which is fine, as
    // for purposes of this code that's functionally equivalent to bona fide alternate alleles.
    //
    // We do not attempt to calculate the last 15 bases of an alternate allele.  These
    // bases are just copied.  If we want to really know these bases we may need to
    // have actual genomic sequences.  That would also provide a control for these calculations.
    //
    // Tried to do this for J segments too but was unsuccessful.
    //
    // Limitations and to do items:
    // 1. Hypothetically we could make a library of alternate alleles and use that
    //    as a supplement or in place of calculating on the fly.
    // 2. Document as part of algorithm.
    // 3. Make alt_refs into a more efficient data structure.
    // 4. Speed up.

    let ta = Instant::now();
    let mut alt_refs = Vec::<(usize, usize, DnaString)>::new(); // (donor, ref id, alt seq)

    // Organize data by reference id.  Note that we ignore exact subclonotypes having four chains.

    let mut allxy =
        vec![Vec::<(usize, Vec<u8>, Vec<usize>, usize, usize)>::new(); refdata.refs.len()];
    for (m, x) in exact_clonotypes.iter().enumerate() {
        if x.share.len() >= 2 && x.share.len() <= 3 {
            for j in 0..x.share.len() {
                let y = &x.share[j];
                let id = y.v_ref_id;

                // Find the partner chains.

                let mut partner = Vec::<usize>::new();
                for ja in 0..x.share.len() {
                    if x.share[ja].left != x.share[j].left {
                        partner.push(x.share[ja].v_ref_id);
                    }
                }
                if !partner.is_empty() {
                    if y.seq_del.len() >= refdata.refs[id].len() - ctl.heur.ref_v_trim {
                        for l in 0..x.clones.len() {
                            let donor = ctl.sample_info.donor_index[x.clones[l][j].lena_index];
                            allxy[id].push((
                                donor,
                                y.seq_del.clone(),
                                partner.clone(),
                                m,
                                x.clones[l][j].lena_index,
                            ));
                        }
                    }
                }
            }
        }
    }

    // Process each reference id.

    for id in 0..refdata.refs.len() {
        if !refdata.is_v(id) {
            continue;
        }
        let mut allx = allxy[id].clone();

        // Divide by donor.

        allx.sort();
        let mut alls = Vec::<Vec<(usize, Vec<u8>, Vec<usize>, usize, usize)>>::new();
        let mut i = 0;
        while i < allx.len() {
            let j = next_diff1_5(&allx, i as i32) as usize;
            let mut all = Vec::<(usize, Vec<u8>, Vec<usize>, usize, usize)>::new();
            for k in i..j {
                all.push(allx[k].clone());
            }
            alls.push(all);
            i = j;
        }

        // Process donor by donor.

        for di in 0..alls.len() {
            // Data here are given by "all", the relevant entries of which are:
            // 1: V..J sequence for one chain of a given info entry
            // 2: the reference id(s) of the partner chain(s) -- possibly not used
            // 3: the index in exact_clonotypes
            // 4: the dataset id.

            let mut all = alls[di].clone();
            let donor_id = all[0].0;

            // If two entries have
            // * the same length CDR3 sequence and
            // * the same partner chain reference ids for both V and J  and
            // * the same length partner chain CDR3 sequence,
            // arbitrarily delete one of them.
            //
            // What we're trying to do here is reduce the incidence wherein multiple entries from
            // the same clonotype (which we have not yet computed) contribute to the count for a
            // putative alternate allele, resulting in alternate alleles that are not true
            // germline changes, but rather arise from SHM.
            //
            // In a better system we might first compute the clonotypes, then the alternate
            // alleles, however knowing the alternate alleles is in fact prerequisite to computing
            // the clonotypes.
            //
            // Note that a vulnerability of the algorithm is that if there is a very large
            // clonotype, then artifactual pairs arising from it could provide enough "evidence"
            // to create an alternate allele, which in fact should not exist.  This has not been
            // observed, but we haven't looked carefully.

            let mut to_delete = vec![false; all.len()];
            {
                let mut trace = Vec::<((usize, usize, usize, usize), usize)>::new();
                for i in 0..all.len() {
                    let u = all[i].3;
                    let ex = &exact_clonotypes[u];
                    for j1 in 0..ex.share.len() {
                        if ex.share[j1].seq_del == all[i].1 {
                            for j2 in 0..ex.share.len() {
                                let (s1, s2) = (&ex.share[j1], &ex.share[j2]);
                                if s2.left != s1.left {
                                    trace.push((
                                        (
                                            s1.cdr3_dna.len(),
                                            s2.cdr3_dna.len(),
                                            s2.v_ref_id,
                                            s2.j_ref_id,
                                        ),
                                        i,
                                    ));
                                }
                            }
                        }
                    }
                }
                trace.sort();
                let mut i = 0;
                while i < trace.len() {
                    let j = next_diff1_2(&trace, i as i32) as usize;
                    for k in i..j {
                        if !to_delete[trace[k].1] {
                            for l in i..j {
                                if l != k {
                                    to_delete[trace[l].1] = true;
                                }
                            }
                            break;
                        }
                    }
                    i = j;
                }
                erase_if(&mut all, &to_delete);
            }

            // Traverse the positions in the reference V segment.

            let mut ps = Vec::<usize>::new();
            for p in 0..refdata.refs[id].len() - ctl.heur.ref_v_trim {
                // Let bases = {(contig base, contig index, index of the other ref V segment)}.
                // The other ref V segments are there only for diagnostic purposes.

                let mut bases = Vec::<(u8, usize, Vec<usize>)>::new();
                for (i, x) in all.iter().enumerate() {
                    bases.push((x.1[p], i, x.2.clone()));
                }
                bases.sort();

                // Let freqs = {( number of contig bases,
                //                {contig indices}, {other ref V segment indices}, base )}.

                let mut freqs = Vec::<(usize, Vec<usize>, Vec<Vec<usize>>, u8)>::new();
                let mut i = 0;
                while i < bases.len() {
                    let j = next_diff1_3(&bases, i as i32) as usize;
                    let mut x = Vec::<usize>::new();
                    let mut y = Vec::<Vec<usize>>::new();
                    for m in i..j {
                        x.push(bases[m].1);
                        y.push(bases[m].2.clone());
                    }
                    freqs.push((j - i, x, y, bases[i].0));
                    i = j;
                }
                reverse_sort(&mut freqs);

                // If frequency of second most frequent base is high enough, save
                // the position in "ps".  Likewise if frequency
                // of first most frequent base is high enough and it's non-reference.
                //
                // Note that this doesn't allow for the case where the third most frequent base
                // is frequent enough.

                let x = refdata.refs[id].get(p);
                let c;
                if x == 0 {
                    c = b'A';
                } else if x == 1 {
                    c = b'C';
                } else if x == 2 {
                    c = b'G';
                } else {
                    c = b'T';
                }
                if (freqs.len() > 0 && freqs[0].0 >= ctl.allele_alg_opt.min_alt && freqs[0].3 != c)
                    || (freqs.len() > 1
                        && ctl.allele_alg_opt.min_mult * freqs[1].0 >= bases.len()
                        && freqs[1].0 >= ctl.allele_alg_opt.min_alt)
                {
                    ps.push(p);
                }
            }
            if ps.len() == 0 {
                continue;
            }
            if ctl.allele_print_opt.con_trace {
                println!(
                    "di = {}, id = {}, transcript = {}, ps = {}",
                    di,
                    id,
                    refdata.transcript[id],
                    ps.iter().format(",")
                );
            }

            // Define types = { ( (contig base at each position in ps), index of contig ) },
            // and sort.

            let mut types = Vec::<(Vec<u8>, usize)>::new();
            for loc in 0..all.len() {
                let mut t = Vec::<u8>::new();
                for i in 0..ps.len() {
                    let p = ps[i];
                    let base = all[loc].1[p];
                    t.push(base);
                }
                types.push((t, loc));
            }
            types.sort();

            // Traverse the types, grouping contigs that have an identical footprint at
            // the positions in ps.

            let mut keep = Vec::<(Vec<u8>, usize, f64, bool)>::new();
            let mut i = 0;
            let mut have_ref = false;
            while i < types.len() {
                let j = next_diff1_2(&types, i as i32) as usize;

                // Determine if the contigs equal reference at the positions in ps.

                let mut is_ref = true;
                for (z, p) in ps.iter().enumerate() {
                    let x = refdata.refs[id].get(*p);
                    let c;
                    if x == 0 {
                        c = b'A';
                    } else if x == 1 {
                        c = b'C';
                    } else if x == 2 {
                        c = b'G';
                    } else {
                        c = b'T';
                    }
                    if c != types[i].0[z] {
                        is_ref = false;
                    }
                }
                if (j - i >= ctl.allele_alg_opt.min_alt
                    && j - i >= types.len() / ctl.allele_alg_opt.min_mult)
                    || is_ref
                {
                    let mut q = Vec::<Vec<usize>>::new();
                    for k in i..j {
                        let m = types[k].1;
                        q.push(all[m].2.clone());
                    }
                    q.sort();
                    let (mut m, mut r) = (0, 0);
                    while r < q.len() {
                        let s = next_diff(&q, r);
                        m = max(m, s - r);
                        r = s;
                    }
                    let purity = 100.0 - percent_ratio(m, q.len());
                    keep.push((types[i].0.clone(), j - i, purity, is_ref));
                    if is_ref {
                        have_ref = true;
                    }
                }
                i = j;
            }
            if keep.len() > 1 || (keep.len() > 0 && !have_ref) {
                // Remove columns that are pure reference.

                let mut to_delete = vec![false; keep[0].0.len()];
                for i in 0..keep[0].0.len() {
                    let mut is_ref = true;
                    for j in 0..keep.len() {
                        let p = ps[i];
                        let x = refdata.refs[id].get(p);
                        let c;
                        if x == 0 {
                            c = b'A';
                        } else if x == 1 {
                            c = b'C';
                        } else if x == 2 {
                            c = b'G';
                        } else {
                            c = b'T';
                        }
                        if c != keep[j].0[i] {
                            is_ref = false;
                        }
                    }
                    if is_ref {
                        to_delete[i] = true;
                    }
                }
                erase_if(&mut ps, &to_delete);
                for j in 0..keep.len() {
                    erase_if(&mut keep[j].0, &to_delete);
                }
                let mut keep0 = Vec::<Vec<u8>>::new();
                for i in 0..keep.len() {
                    keep0.push(keep[i].0.clone());
                }
                keep0.sort();
                keep.sort_by(|a, b| a.partial_cmp(b).unwrap());

                // Print.

                if ctl.allele_print_opt.con {
                    println!("\nDONOR {}", donor_id + 1);
                    println!("{} = |{}| = {}", id, refdata.id[id], refdata.name[id]);
                    println!("ps = {}", ps.iter().format(","));
                    for x in keep.iter() {
                        let mut bases = String::new();
                        for z in x.0.iter() {
                            bases.push(*z as char);
                        }
                        print!("{} [{}] {:.1}", bases, x.1, x.2);
                        if x.3 {
                            print!(" (ref)");
                        }
                        println!("");
                    }
                }

                // Save alternate references.

                for x in keep.iter() {
                    if !x.3 {
                        let mut b = refdata.refs[id].clone();
                        for i in 0..ps.len() {
                            let c;
                            if x.0[i] == b'A' {
                                c = 0;
                            } else if x.0[i] == b'C' {
                                c = 1;
                            } else if x.0[i] == b'G' {
                                c = 2;
                            } else {
                                c = 3;
                            }
                            b.set_mut(ps[i], c);
                        }
                        alt_refs.push((donor_id, id, b));
                    }
                }

                // For each lena id from the donor, classify the alleles from that lena
                // id according to the classification just derived.  Print the matrix.

                let all = &alls[di];
                let ll = &ctl.sample_info.dataset_list[donor_id];
                let mut count = vec![vec![0; keep0.len()]; ll.len()];
                for i in 0..all.len() {
                    let lid = all[i].4;
                    let mut allele = Vec::<u8>::new();
                    for p in ps.iter() {
                        allele.push(all[i].1[*p]);
                    }
                    let m = bin_position(&keep0, &allele);
                    if m >= 0 {
                        let l = bin_position(
                            &ctl.sample_info.dataset_list[donor_id],
                            // &ctl.sample_info.dataset_id[lid].force_usize());
                            &lid,
                        );
                        count[l as usize][m as usize] += 1;
                    }
                }
                /*
                if ctl.allele_print_opt.con {
                    println!("MATRIX");
                    for l in 0..ll.len() {
                        println!( "{} ==> {}", ll[l], count[l].iter().format("\t") );
                    }
                }
                */
            }
        }
    }
    if ctl.comp {
        println!(
            "used {:.2} seconds used finding alt alleles, peak mem = {:.2} GB",
            elapsed(&ta),
            peak_mem_usage_gb()
        );
    }
    alt_refs.sort();
    alt_refs
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Update reference sequences for V segments by substituting in alt alleles if better.
// Computational performance dubious because of full alt_refs traversal.

pub fn sub_alts(
    ctl: &EncloneControl,
    alt_refs: &Vec<(usize, usize, DnaString)>,
    info: &mut Vec<CloneInfo>,
    exact_clonotypes: &mut Vec<ExactClonotype>,
) {
    let t = Instant::now();
    for i in 0..info.len() {
        for j in 0..info[i].vs.len() {
            if info[i].vs[j].len() - ctl.heur.ref_v_trim <= info[i].tigs[j].len() {
                let mut errs = 0;
                for l in 0..info[i].vs[j].len() - ctl.heur.ref_v_trim {
                    let x = info[i].vs[j].get(l);
                    let c;
                    if x == 0 {
                        c = b'A';
                    } else if x == 1 {
                        c = b'C';
                    } else if x == 2 {
                        c = b'G';
                    } else {
                        c = b'T';
                    }
                    if info[i].tigs[j][l] != c {
                        errs += 1;
                    }
                }
                let mut donors = Vec::<usize>::new();
                for lena in info[i].origin.iter() {
                    donors.push(ctl.sample_info.donor_index[*lena]);
                }
                for donor in donors {
                    for m in 0..alt_refs.len() {
                        if alt_refs[m].0 == donor && alt_refs[m].1 == info[i].vsids[j] {
                            if alt_refs[m].2.len() - ctl.heur.ref_v_trim <= info[i].tigs[j].len() {
                                let mut alt_errs = 0;
                                for l in 0..alt_refs[m].2.len() - ctl.heur.ref_v_trim {
                                    let x = alt_refs[m].2.get(l);
                                    let c;
                                    if x == 0 {
                                        c = b'A';
                                    } else if x == 1 {
                                        c = b'C';
                                    } else if x == 2 {
                                        c = b'G';
                                    } else {
                                        c = b'T';
                                    }
                                    if info[i].tigs[j][l] != c {
                                        alt_errs += 1;
                                    }
                                }
                                if alt_errs < errs {
                                    info[i].vs[j] = alt_refs[m].2.clone();
                                    info[i].dref[j] = Some(m); // not sure we're actually using this
                                    let ex = &mut exact_clonotypes[info[i].clonotype_id];
                                    for z in 0..ex.share.len() {
                                        if ex.share[z].seq == info[i].tigs[j] {
                                            ex.share[z].v_ref_id_donor = Some(m);
                                            ex.share[z].v_ref_id_donor_donor = Some(donor);
                                            let mut alts = 0;
                                            let mut mm = m;
                                            while mm >= 1 {
                                                mm -= 1;
                                                if alt_refs[mm].0 == donor
                                                    && alt_refs[mm].1 == alt_refs[m].1
                                                {
                                                    alts += 1;
                                                }
                                            }
                                            ex.share[z].v_ref_id_donor_alt_id = Some(alts);
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    if ctl.comp {
        println!(
            "used {:.2} seconds used substituting alt alleles, peak mem = {:.2} GB",
            elapsed(&t),
            peak_mem_usage_gb()
        );
    }
}