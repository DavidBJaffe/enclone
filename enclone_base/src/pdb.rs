// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// Tools to work with PDB files.

use amino::*;
use flate2::read::MultiGzDecoder;
use io_utils::*;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::process::Command;
use string_utils::*;
use vector_utils::*;

#[derive(Default, Clone)]
pub struct PdbStructure {
    pub chain_names: Vec<String>,
    pub chains: Vec<Vec<u8>>,
    pub atoms: Vec<PdbAtom>,
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// PdbAtom: for the first two fields, the identifiers are left justified and zero-byte filed.

#[derive(Default, Clone)]
pub struct PdbAtom {
    pub atom_name: [u8; 2], // e.g. C for carbon etc.
    pub atom_alt: [u8; 4],  // identifier for the atom location within the amino acid
    pub amino_acid: u8,     // e.g. L for leucine etc.
    pub chain: u16,         // zero-based chain number
    pub chain_pos: u16,     // zero-based index of amino acid on chain
    pub x: f32,             // x coordinate of atom
    pub y: f32,             // y coordinate of atom
    pub z: f32,             // y coordinate of atom
}

// PdbAtomoid same as PdbAtom except that chain_pos is a string, could be e.g. 82A.

#[derive(Default, Clone)]
pub struct PdbAtomoid {
    pub atom_name: [u8; 2], // e.g. C for carbon etc.
    pub atom_alt: [u8; 4],  // identifier for the atom location within the amino acid
    pub amino_acid: u8,     // e.g. L for leucine etc.
    pub chain: u16,         // zero-based chain number
    pub chain_pos: String,  // string representing index of amino acid on chain
    pub x: f32,             // x coordinate of atom
    pub y: f32,             // y coordinate of atom
    pub z: f32,             // y coordinate of atom
}

impl PdbStructure {
    pub fn fetch_atoms_range(&self, chain: usize, start: usize, stop: usize) -> Vec<[f32; 3]> {
        let mut u = Vec::<[f32; 3]>::new();
        for j in 0..self.atoms.len() {
            if self.atoms[j].chain as usize == chain {
                if (self.atoms[j].chain_pos as usize) >= start
                    && (self.atoms[j].chain_pos as usize) < stop
                {
                    let x_s = self.atoms[j].x;
                    let y_s = self.atoms[j].y;
                    let z_s = self.atoms[j].z;
                    u.push([x_s, y_s, z_s]);
                }
            }
        }
        return u;
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Read an atom from an alternately formatted "PDB" file.  In this format, there is no field for
// entity, which has to be inferred from the file as a whole.   We pass it as an argument.

pub fn parse_atom_pidgin(line: &str, chain_count: usize) -> Result<PdbAtomoid, String> {
    let fields = line.split_ascii_whitespace().collect::<Vec<&str>>();
    let atom0 = fields[11].as_bytes();
    let mut atom = [0 as u8; 2];
    atom[0] = atom0[0];
    if atom0.len() == 2 {
        atom[1] = atom0[1];
    }
    if atom0.len() > 2 {
        return Err(format!(
            "Error reading PDB file, atom name appears to have more than two characters in this \
                line:\n{line}",
        ));
    }
    let atom_alt0 = fields[2].as_bytes();
    if atom_alt0.len() > 4 {
        return Err(format!(
            "atom_alt0 = {} has length {}, but should be at most 4",
            fields[2],
            atom_alt0.len()
        ));
    }
    let mut atom_alt = [0 as u8; 4];
    for j in 0..atom_alt0.len() {
        atom_alt[j] = atom_alt0[j];
    }

    // Check for defective amino acid.  Possibly these cases could be rescued.

    if fields[3].len() != 3 || fields[3] == "UNK" {
        return Err(format!(
            "in pdb file, unknown amino acid {} in line {}",
            fields[3], line
        ));
    }

    // Create atom.

    let aa = aa3_to_aa(&fields[3].as_bytes());
    let pos_on_entity = fields[5].to_string();
    let entity = chain_count as u16;
    let xcoord = fields[6].force_f64() as f32;
    let ycoord = fields[7].force_f64() as f32;
    let zcoord = fields[8].force_f64() as f32;
    let m = PdbAtomoid {
        atom_name: atom,
        atom_alt: atom_alt,
        amino_acid: aa,
        chain: entity,
        chain_pos: pos_on_entity,
        x: xcoord,
        y: ycoord,
        z: zcoord,
    };
    Ok(m)
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn parse_atom(line: &str) -> Result<PdbAtom, String> {
    let fields = line.split_ascii_whitespace().collect::<Vec<&str>>();
    let atom0 = fields[2].as_bytes();
    let mut atom = [0 as u8; 2];
    atom[0] = atom0[0];
    if atom0.len() == 2 {
        atom[1] = atom0[1];
    }
    assert!(atom0.len() <= 2);
    let atom_alt0 = fields[3].as_bytes();
    if atom_alt0.len() > 4 {
        return Err(format!(
            "atom_alt0 = {} has length {}, but should be at most 4",
            fields[3],
            atom_alt0.len()
        ));
    }
    let mut atom_alt = [0 as u8; 4];
    for j in 0..atom_alt0.len() {
        atom_alt[j] = atom_alt0[j];
    }

    // Check for defective amino acid.  Possibly these cases could be rescued.

    if fields[5].len() != 3 || fields[5] == "UNK" {
        return Err(format!(
            "in pdb, unknown amino acid {} in line {}",
            fields[5], line
        ));
    }

    // Create atom.

    let aa = aa3_to_aa(&fields[5].as_bytes());
    let pos_on_entity = fields[8].force_usize() as u16;
    let entity = fields[7].force_usize() as u16;
    let xcoord = fields[10].force_f64() as f32;
    let ycoord = fields[11].force_f64() as f32;
    let zcoord = fields[12].force_f64() as f32;
    let m = PdbAtom {
        atom_name: atom,
        atom_alt: atom_alt,
        amino_acid: aa,
        chain: entity - 1,
        chain_pos: pos_on_entity - 1,
        x: xcoord,
        y: ycoord,
        z: zcoord,
    };
    Ok(m)
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn fetch_pdb_structure(name: &str) -> Result<PdbStructure, String> {
    fetch_pdb_structure_gen(&name, "antibody_sets/pdbs")
}

pub fn fetch_pdb_structure_gen(name: &str, dir: &str) -> Result<PdbStructure, String> {
    let mut lines = Vec::<String>::new();
    let opath = format!("{}/{}.gz", dir, name);
    if path_exists(&opath) {
        let gz = MultiGzDecoder::new(File::open(&opath).unwrap());
        let b = BufReader::new(gz);
        for line in b.lines() {
            let s = line.unwrap();
            lines.push(s);
        }
    } else {
        let x = Command::new("wget")
            .arg("-q")
            .arg("-O")
            .arg("-")
            .arg(&format!("https://files.rcsb.org/download/{}.cif", name))
            .output()
            .expect(&format!("failed to execute wget"));
        for line in strme(&x.stdout).lines() {
            lines.push(line.to_string());
        }
    }
    if lines.is_empty() {
        eprintln!("");
        eprintme!(name, dir);
        eprintln!("No lines read.\n");
        std::process::exit(1);
    }
    let mut chain_names = Vec::<String>::new();
    let mut chains = Vec::<Vec<u8>>::new();
    let mut ms = Vec::<PdbAtom>::new();
    let mut i = 0;
    while i < lines.len() {
        // Parse the chain names.  This is madness.

        if lines[i].starts_with("_entity.details") {
            let mut elines = Vec::<String>::new();
            let mut j = i + 1;
            while !lines[j].starts_with('#') {
                elines.push(lines[j].to_string());
                j += 1;
            }
            let mtypes = ["polymer", "non-polymer", "branched", "water", "?"];
            let mut count = 1;
            let mut starts = Vec::<usize>::new();
            for m in 0..elines.len() {
                let mut s = elines[m].clone();
                while s.contains("  ") {
                    s = s.replace("  ", " ");
                }
                if s.starts_with(&format!("{} ", count)) {
                    if s.after(" ").contains(" ") {
                        let mtype = s.between(" ", " ");
                        let mut known = false;
                        for u in 0..mtypes.len() {
                            if mtype == mtypes[u] {
                                known = true;
                            }
                        }
                        if !known {
                            eprintln!("unknown type {}", mtype);
                            std::process::exit(1);
                        }
                    }
                    starts.push(m);
                    count += 1;
                }
            }
            starts.push(elines.len());
            let mut flines = Vec::<String>::new();
            for z in 0..starts.len() - 1 {
                let mut s = String::new();
                for k in starts[z]..starts[z + 1] {
                    s += &elines[k];
                }
                flines.push(s);
            }
            for m in 0..flines.len() {
                let mut s = flines[m].replace("  ", " ");
                while s.contains("  ") {
                    s = s.replace("  ", " ");
                }
                if s.between(" ", " ") == "polymer" {
                    let mut t = s.after("polymer").to_string();
                    if t.contains("'") && t.after("'").contains("'") {
                        t = t.between("'", "'").to_string();
                    }
                    chain_names.push(t);
                }
            }
            i = j + 1;

        // Parse the chain sequences.
        } else if lines[i].starts_with("_entity_poly.pdbx_target_identifier") {
            let mut j = i + 1;
            while j < lines.len() && !lines[j].starts_with('#') {
                if lines[j].contains("'polypeptide(L)'") {
                    if lines[j].ends_with("? ") {
                        let s = lines[j]
                            .after("no no ")
                            .between(" ", " ")
                            .as_bytes()
                            .to_vec();
                        chains.push(s);
                        j += 1;
                        continue;
                    }
                    j += 1;
                    if lines[j].starts_with("#") {
                        break;
                    }
                    if !lines[j].contains(";") {
                        let s = Vec::new();
                        chains.push(s);
                        continue;
                    }

                    // Skip the first entry.  For some reason the second one works better.

                    for _ in 0..2 {
                        while j < lines.len() {
                            j += 1;
                            if lines[j].starts_with(";") {
                                break;
                            }
                        }
                    }

                    // Extract the second entry.

                    let mut s = lines[j].after(";").as_bytes().to_vec();
                    while j < lines.len() {
                        j += 1;
                        if lines[j].starts_with(";") {
                            break;
                        }
                        s.append(&mut lines[j].as_bytes().to_vec());
                    }
                    chains.push(s);
                } else {
                    j += 1;
                }
            }
            i = j;
        // Parse an alternate representation of the chain sequences, present in some PDB files.
        } else if lines[i].starts_with("_entity_poly.pdbx_seq_one_letter_code ") {
            let mut j = i + 1;
            while j < lines.len() && lines[j].starts_with(';') {
                if lines[j] == ";" {
                    break;
                }
                let mut s = lines[j].after(";").as_bytes().to_vec();
                while j < lines.len() {
                    j += 1;
                    if lines[j].starts_with(";") {
                        break;
                    }
                    s.append(&mut lines[j].as_bytes().to_vec());
                }
                chains.push(s);
            }
            i = j;
        } else if lines[i].starts_with("ATOM ") {
            let m = parse_atom(&lines[i]).unwrap();
            ms.push(m);
            i += 1;
        } else {
            i += 1;
        }
    }
    let mut to_delete = vec![false; chains.len()];
    for i in 0..chains.len() {
        if chains[i].len() == 0 {
            to_delete[i] = true;
        }
    }
    erase_if(&mut chains, &to_delete);
    erase_if(&mut chain_names, &to_delete);
    if chains.len() != chain_names.len() {
        return Err(format!(
            "in pdb {name}, number of chains = {} but number of chain names = {}",
            chains.len(),
            chain_names.len(),
        ));
    } else {
        for i in 0..chains.len() {
            loop {
                let mut found = false;
                for j in 0..chains[i].len() - 4 {
                    if chains[i][j] == b'(' && chains[i][j + 4] == b')' {
                        let s = stringme(&chains[i]);
                        let t = s.replace(&strme(&chains[i][j..j + 4]), "X");
                        chains[i] = t.as_bytes().to_vec();
                        found = true;
                        break;
                    }
                }
                if !found {
                    break;
                }
            }
            if chains[i].contains(&b'(') {
                eprintln!("\nProblem with {} in {}.\n", strme(&chains[i]), name);
                std::process::exit(1);
            }
        }
        Ok(PdbStructure {
            chain_names: chain_names,
            chains: chains,
            atoms: ms,
        })
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn fetch_pdb_structure_pidgin(name: &str, dir: &str) -> Result<PdbStructure, String> {
    let mut lines = Vec::<String>::new();
    let opath = format!("{}/{}.gz", dir, name);
    let ppath = format!("{}/{}", dir, name);
    if path_exists(&opath) {
        let gz = MultiGzDecoder::new(File::open(&opath).unwrap());
        let b = BufReader::new(gz);
        for line in b.lines() {
            let s = line.unwrap();
            lines.push(s);
        }
    } else if path_exists(&ppath) {
        let f = open_for_read![&ppath];
        for line in f.lines() {
            let s = line.unwrap();
            lines.push(s);
        }
    } else {
        return Err(format!("Unable to find PDB file {name}."));
    }
    if lines.is_empty() {
        eprintln!("");
        eprintme!(name, dir);
        eprintln!("No lines read.\n");
        std::process::exit(1);
    }
    let mut ms0 = Vec::<PdbAtomoid>::new();
    let mut chain_count = 0;
    for i in 0..lines.len() {
        if lines[i].starts_with("ATOM ") {
            let m = parse_atom_pidgin(&lines[i], chain_count).unwrap();
            ms0.push(m);
        } else if lines[i].starts_with("TER ") {
            chain_count += 1;
        } else if !lines[i].starts_with("END") {
            return Err(format!("Don't understand {}.", lines[i]));
        }
    }

    let mut pos = 0;
    let mut ms = Vec::<PdbAtom>::new();
    for i in 0..ms0.len() {
        if i > 0 && ms0[i].chain != ms0[i - 1].chain {
            pos = 0;
        } else if i > 0 && ms0[i].chain_pos != ms0[i - 1].chain_pos {
            pos += 1;
        }
        let m = PdbAtom {
            atom_name: ms0[i].atom_name,
            atom_alt: ms0[i].atom_alt,
            amino_acid: ms0[i].amino_acid,
            chain: ms0[i].chain,
            chain_pos: pos,
            x: ms0[i].x,
            y: ms0[i].y,
            z: ms0[i].z,
        };
        ms.push(m);
    }

    let mut chains = vec![Vec::<u8>::new(); chain_count];
    let c = [b'C', 0, 0, 0];
    for x in ms.iter() {
        if x.atom_alt == c {
            chains[x.chain as usize].push(x.amino_acid);
        }
    }
    Ok(PdbStructure {
        chain_names: Vec::new(),
        chains: chains,
        atoms: ms,
    })
}
