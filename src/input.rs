use crate::cli::Analyze;
use crate::invalid;
use crate::model::{AccessibilityMask, FileIdentity, Selection};
use crate::neutral::windows::position_windows;
use flate2::read::MultiGzDecoder;
use rsomics_common::Result;
use sha2::{Digest, Sha256};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

pub(crate) fn reader(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path)?;
    let input: Box<dyn Read> = if path.extension().is_some_and(|value| value == "gz") {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };
    Ok(Box::new(BufReader::with_capacity(1024 * 1024, input)))
}

pub(crate) fn identity(path: &Path) -> Result<FileIdentity> {
    let mut input = File::open(path)?;
    let mut hash = Sha256::new();
    let mut buffer = vec![0; 1024 * 1024];
    loop {
        let n = input.read(&mut buffer)?;
        if n == 0 {
            break;
        }
        hash.update(&buffer[..n]);
    }
    Ok(FileIdentity {
        size: input.metadata()?.len(),
        sha256: format!("{:x}", hash.finalize()),
    })
}

pub(crate) struct Input {
    pub reader: Box<dyn BufRead>,
    pub metadata: Vec<String>,
    pub selection: Selection,
    pub masks: Vec<AccessibilityMask>,
    pub files: BTreeMap<String, FileIdentity>,
}

pub(crate) fn prepare(args: &Analyze, budget: usize) -> Result<Input> {
    if args.population.len() != 2 || args.population[0] == args.population[1] {
        return Err(invalid(
            "exactly two distinct population labels are required",
        ));
    }
    if args.window_bp == 0 || args.jackknife_variants == 0 {
        return Err(invalid("window and jackknife sizes must be positive"));
    }
    let mut panel = HashMap::new();
    let mut seen = HashSet::new();
    let mut panel_reader = reader(&args.panel)?;
    let mut line = String::new();
    panel_reader.read_line(&mut line)?;
    let header: Vec<_> = line.split_whitespace().collect();
    if header.len() < 2 || header[..2] != ["sample", "pop"] {
        return Err(invalid("panel header must start with sample and pop"));
    }
    let columns = header.len();
    for line in panel_reader.lines() {
        let line = line?;
        let values: Vec<_> = line.split_whitespace().collect();
        if values.len() != columns || !seen.insert(values[0].to_owned()) {
            return Err(invalid("panel identity/cardinality differs"));
        }
        if let Some(population) = args.population.iter().position(|name| name == values[1]) {
            panel.insert(values[0].to_owned(), population);
        }
    }
    let mut input = reader(&args.vcf)?;
    let mut metadata = Vec::new();
    let samples = loop {
        line.clear();
        if input.read_line(&mut line)? == 0 {
            return Err(invalid("VCF #CHROM header is absent"));
        }
        let value = line.trim_end_matches(['\r', '\n']);
        if value.starts_with("##") {
            metadata.push(value.to_owned());
            continue;
        }
        let fields: Vec<_> = value.split('\t').collect();
        if fields.len() < 10
            || fields[..9]
                != [
                    "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
                ]
        {
            return Err(invalid("VCF sample header schema differs"));
        }
        break fields[9..]
            .iter()
            .map(|value| (*value).to_owned())
            .collect::<Vec<_>>();
    };
    if !metadata
        .iter()
        .any(|line| line.starts_with("##FORMAT=<ID=GT,Number=1,Type=String,"))
        || !metadata
            .iter()
            .any(|line| line.starts_with("##INFO=<ID=AA,Number=1,Type=String,"))
    {
        return Err(invalid("VCF requires the admitted GT/AA schemas"));
    }
    if samples.iter().collect::<HashSet<_>>().len() != samples.len() {
        return Err(invalid("duplicate VCF sample identity"));
    }
    let mut selection = Selection {
        names: [args.population[0].clone(), args.population[1].clone()],
        ids: Vec::new(),
        original_columns: Vec::new(),
        population_of_original: vec![None; samples.len()],
        population_indices: [Vec::new(), Vec::new()],
        original_samples: samples.len(),
    };
    for (column, sample) in samples.iter().enumerate() {
        if let Some(&population) = panel.get(sample) {
            selection.population_indices[population].push(selection.ids.len());
            selection.ids.push(sample.clone());
            selection.original_columns.push(column);
            selection.population_of_original[column] = Some(population);
        }
    }
    if selection.ids.len() != panel.len() || selection.population_indices.iter().any(Vec::is_empty)
    {
        return Err(invalid("a selected population/sample is absent from VCF"));
    }
    let mut requested = HashMap::<String, PathBuf>::new();
    for value in &args.mask {
        let (name, path) = value
            .split_once('=')
            .ok_or_else(|| invalid("--mask requires CONTIG=PATH"))?;
        if name.is_empty()
            || name == "."
            || name == ".."
            || name.contains(['/', '\\'])
            || requested.insert(name.to_owned(), path.into()).is_some()
        {
            return Err(invalid("duplicate or unsafe mask contig"));
        }
    }
    let mut files = BTreeMap::new();
    files.insert("vcf".into(), identity(&args.vcf)?);
    files.insert("panel".into(), identity(&args.panel)?);
    files.insert("reference_fai".into(), identity(&args.reference_fai)?);
    if let Some(path) = &args.index {
        files.insert("tbi".into(), identity(path)?);
    }
    let mut masks = Vec::new();
    let mut lengths = HashSet::new();
    let mut mask_bytes = 0usize;
    for record in reader(&args.reference_fai)?.lines() {
        let record = record?;
        let values: Vec<_> = record.split('\t').collect();
        if values.len() != 5 || !lengths.insert(values[0].to_owned()) {
            return Err(invalid("reference FAI schema/identity differs"));
        }
        let Some(path) = requested.remove(values[0]) else {
            continue;
        };
        let length: u64 = values[1]
            .parse()
            .map_err(|_| invalid("invalid FAI length"))?;
        if length < 2 {
            return Err(invalid(
                "connected workflow requires contig length at least two",
            ));
        }
        let words = usize::try_from(length.div_ceil(64))
            .map_err(|_| invalid("mask address range overflows"))?;
        mask_bytes = mask_bytes
            .checked_add(
                words
                    .checked_mul(8)
                    .ok_or_else(|| invalid("mask size overflows"))?,
            )
            .ok_or_else(|| invalid("mask size overflows"))?;
        if mask_bytes > budget / 4 {
            return Err(invalid("accessibility masks exceed the planning budget"));
        }
        if length
            .div_ceil(args.window_bp)
            .checked_mul(12 * 1024)
            .is_none_or(|value| value > budget as u64 / 3)
        {
            return Err(invalid("coordinate windows exceed the planning budget"));
        }
        let windows = position_windows(length, args.window_bp);
        let mut mask = AccessibilityMask {
            name: values[0].into(),
            length,
            bits: vec![0; words],
            codes: BTreeMap::new(),
            window_accessible: vec![0; windows.len()],
            windows,
        };
        let mut mask_reader = reader(&path)?;
        line.clear();
        mask_reader.read_line(&mut line)?;
        let header = line
            .trim_end()
            .strip_prefix('>')
            .and_then(|value| value.split_whitespace().next())
            .ok_or_else(|| invalid("mask FASTA header is absent"))?;
        if header != mask.name && header != format!("chr{}", mask.name) {
            return Err(invalid("mask FASTA contig differs"));
        }
        let mut position = 0usize;
        let mut counts = [0u64; 256];
        let mut window = 0usize;
        for row in mask_reader.lines() {
            let row = row?;
            if row.starts_with('>') {
                return Err(invalid("mask FASTA contains multiple contigs"));
            }
            for code in row.trim().bytes() {
                if !b"PNLHZQ0".contains(&code) || position as u64 >= length {
                    return Err(invalid("mask length/code differs"));
                }
                counts[code as usize] += 1;
                while position as u64 + 1 > mask.windows[window].stop {
                    window += 1;
                }
                if code == b'P' {
                    mask.bits[position / 64] |= 1 << (position % 64);
                    mask.window_accessible[window] += 1;
                }
                position += 1;
            }
        }
        if position as u64 != length {
            return Err(invalid("mask length differs from FAI"));
        }
        mask.codes = counts
            .iter()
            .enumerate()
            .filter(|(_, count)| **count > 0)
            .map(|(code, &count)| ((code as u8 as char).to_string(), count))
            .collect();
        files.insert(format!("mask:{}", mask.name), identity(&path)?);
        masks.push(mask);
    }
    if !requested.is_empty() || masks.is_empty() {
        return Err(invalid("mask contig is absent from reference FAI"));
    }
    Ok(Input {
        reader: input,
        metadata,
        selection,
        masks,
        files,
    })
}
