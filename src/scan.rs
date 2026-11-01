use crate::input::Input;
use crate::model::{
    AccessibilityMask, CountRow, DiploidCall, FileIdentity, LEGACY_PROFILE, LegacyIdentity,
    PopulationCounts, Selection,
};
use crate::{invalid, json_error};
use rayon::prelude::*;
use rsomics_common::Result;
use serde_json::{Value, json};
use sha2::{Digest, Sha256};
use std::collections::BTreeMap;
use std::io::{BufRead, BufWriter, Write};
use tempfile::NamedTempFile;

pub(crate) struct Scan {
    pub input: Input,
    pub spool: NamedTempFile,
    pub max_alt: usize,
    pub rows: u64,
    pub retained: Vec<usize>,
    pub original_rows: Vec<u64>,
    pub max_an: Vec<[u64; 2]>,
    pub window_rows: Vec<Vec<usize>>,
    pub filters: BTreeMap<String, u64>,
    pub gt_sha: String,
    pub record_sha: String,
    pub legacy: LegacyIdentity,
}

struct Parsed {
    row: CountRow,
    reason: &'static str,
    gt: Vec<u8>,
    identity: Vec<u8>,
    legacy: Option<Vec<u8>>,
}

fn genotype(value: &str, alt_count: usize) -> Result<DiploidCall> {
    let separators: Vec<_> = value
        .match_indices(['/', '|'])
        .map(|(index, _)| index)
        .collect();
    if separators.len() != 1 {
        return Err(invalid("GT is not diploid"));
    }
    let index = separators[0];
    let first = &value[..index];
    let second = &value[index + 1..];
    if first == "." && second == "." {
        return Ok(DiploidCall::Missing);
    }
    if first == "." || second == "." {
        return Err(invalid("partial-missing GT"));
    }
    let parse = |value: &str| -> Result<u16> {
        if value.is_empty() || !value.bytes().all(|value| value.is_ascii_digit()) {
            return Err(invalid("malformed genotype allele"));
        }
        let index: u16 = value
            .parse()
            .map_err(|_| invalid("genotype allele overflows int16 profile"))?;
        if index as usize > alt_count || index > 32767 {
            return Err(invalid("GT allele index exceeds ALT/int16 profile"));
        }
        Ok(index)
    };
    Ok(DiploidCall::Called(parse(first)?, parse(second)?))
}

fn ascii_json(value: &Value) -> Result<Vec<u8>> {
    let value = serde_json::to_string(value).map_err(json_error)?;
    let mut output = String::new();
    for character in value.chars() {
        if character.is_ascii() {
            output.push(character);
        } else {
            use std::fmt::Write;
            let mut encoded = [0; 2];
            for unit in character.encode_utf16(&mut encoded) {
                write!(&mut output, "\\u{unit:04x}").expect("writing into a String");
            }
        }
    }
    output.push('\n');
    Ok(output.into_bytes())
}

fn parse(
    line: &str,
    ordinal: u64,
    selection: &Selection,
    masks: &[AccessibilityMask],
) -> Result<Parsed> {
    let line = line.trim_end_matches(['\r', '\n']);
    let mut fields = line.split('\t');
    let basic: Vec<_> = fields.by_ref().take(9).collect();
    if basic.len() != 9 {
        return Err(invalid("VCF record has fewer than nine fixed columns"));
    }
    let contig = masks
        .iter()
        .position(|mask| mask.name == basic[0])
        .ok_or_else(|| invalid("VCF contig outside admitted masks"))?;
    let position: u64 = basic[1].parse().map_err(|_| invalid("invalid VCF POS"))?;
    if position == 0 || position > masks[contig].length {
        return Err(invalid("VCF coordinate outside admitted contig"));
    }
    let alts: Vec<_> = if basic[4] == "." {
        Vec::new()
    } else {
        basic[4].split(',').collect()
    };
    if alts.len() > 32767 || basic[3].is_empty() || alts.iter().any(|value| value.is_empty()) {
        return Err(invalid("invalid ALT width or empty VCF allele"));
    }
    let formats: Vec<_> = basic[8].split(':').collect();
    let gt_index = formats
        .iter()
        .position(|&value| value == "GT")
        .ok_or_else(|| invalid("FORMAT has no GT"))?;
    if formats
        .iter()
        .collect::<std::collections::HashSet<_>>()
        .len()
        != formats.len()
    {
        return Err(invalid("duplicate FORMAT field"));
    }
    let mut populations = std::array::from_fn(|_| PopulationCounts {
        copies: vec![0; alts.len() + 1],
        heterozygotes: vec![0; alts.len() + 1],
    });
    let mut gt = Vec::with_capacity(selection.ids.len() * 4);
    let mut selected_fields = Vec::new();
    let mut sample_count = 0;
    for (column, field) in fields.enumerate() {
        sample_count += 1;
        let population = selection
            .population_of_original
            .get(column)
            .copied()
            .flatten();
        let Some(population) = population else {
            continue;
        };
        if field.bytes().filter(|&value| value == b':').count() + 1 > formats.len() {
            return Err(invalid("sample FORMAT cardinality differs"));
        }
        let value = field
            .split(':')
            .nth(gt_index)
            .ok_or_else(|| invalid("sample GT field absent"))?;
        let calls = genotype(value, alts.len())?;
        match calls {
            DiploidCall::Missing => {
                gt.extend_from_slice(&(-1i16).to_le_bytes());
                gt.extend_from_slice(&(-1i16).to_le_bytes());
            }
            DiploidCall::Called(first, second) => {
                gt.extend_from_slice(&(first as i16).to_le_bytes());
                gt.extend_from_slice(&(second as i16).to_le_bytes());
                let counts = &mut populations[population];
                counts.copies[first as usize] += 1;
                counts.copies[second as usize] += 1;
                if first != second {
                    counts.heterozygotes[first as usize] += 1;
                    counts.heterozygotes[second as usize] += 1;
                }
            }
        }
        selected_fields.push(field);
    }
    if sample_count != selection.original_samples {
        return Err(invalid("VCF record sample cardinality differs"));
    }
    let mut aa = basic[7]
        .split(';')
        .filter_map(|value| value.strip_prefix("AA="));
    let raw_aa = aa.next().unwrap_or("");
    if aa.next().is_some() {
        return Err(invalid("duplicate INFO/AA"));
    }
    let passed = basic[6] == "PASS";
    let nucleotide = |value: &str| matches!(value, "A" | "C" | "G" | "T");
    let reason = if !passed {
        "not_pass"
    } else if !nucleotide(basic[3])
        || alts.is_empty()
        || !alts.iter().all(|value| nucleotide(value))
    {
        "not_literal_snv"
    } else if !masks[contig].accessible(position) {
        "inaccessible"
    } else {
        "retained"
    };
    let legacy = if reason == "retained"
        && alts.len() == 1
        && populations.iter().enumerate().all(|(pop, value)| {
            value
                .copies
                .iter()
                .map(|&value| value as usize)
                .sum::<usize>()
                == selection.population_indices[pop].len() * 2
        }) {
        let mut value = basic.join("\t");
        value.push('\t');
        value.push_str(&selected_fields.join("\t"));
        value.push('\n');
        Some(value.into_bytes())
    } else {
        None
    };
    let identity = ascii_json(&json!([
        ordinal, basic[0], position, basic[3], alts, passed, raw_aa
    ]))?;
    let mut alleles = vec![basic[3].to_owned()];
    alleles.extend(alts.into_iter().map(str::to_owned));
    Ok(Parsed {
        row: CountRow {
            ordinal,
            contig,
            position,
            alleles,
            raw_aa: raw_aa.to_owned(),
            populations,
        },
        reason,
        gt,
        identity,
        legacy,
    })
}

pub(crate) fn scan(mut input: Input, budget: usize, pool: &rayon::ThreadPool) -> Result<Scan> {
    let spool = NamedTempFile::new_in(std::env::temp_dir())?;
    let mut output = BufWriter::with_capacity(1024 * 1024, spool.reopen()?);
    let mut gt_hash = Sha256::new();
    let mut record_hash = Sha256::new();
    let mut legacy_hash = Sha256::new();
    let mut legacy_bytes = 0u64;
    let mut legacy_rows = 0u64;
    for metadata in &input.metadata {
        legacy_hash.update(metadata.as_bytes());
        legacy_hash.update(b"\n");
        legacy_bytes += metadata.len() as u64 + 1;
    }
    let header = format!(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n",
        input.selection.ids.join("\t")
    );
    legacy_hash.update(header.as_bytes());
    legacy_bytes += header.len() as u64;
    let mut rows = 0u64;
    let mut maximum = 1usize;
    let mut previous = None;
    let mut retained = vec![0; input.masks.len()];
    let mut original_rows = vec![0; input.masks.len()];
    let mut max_an = vec![[0; 2]; input.masks.len()];
    let mut window_rows: Vec<Vec<usize>> = input
        .masks
        .iter()
        .map(|mask| vec![0; mask.windows.len()])
        .collect();
    let mut filters: BTreeMap<String, u64> =
        ["retained", "not_pass", "not_literal_snv", "inaccessible"]
            .into_iter()
            .map(|value| (value.into(), 0))
            .collect();
    let batch_limit = (budget / 16).min(8 * 1024 * 1024);
    let mut done = false;
    while !done {
        let mut batch = Vec::new();
        let mut bytes = 0usize;
        while batch.len() < 8192 && bytes < batch_limit {
            let mut line = String::new();
            let n = input.reader.read_line(&mut line)?;
            if n == 0 {
                done = true;
                break;
            }
            if line.starts_with('#') {
                continue;
            }
            if n > budget / 8 {
                return Err(invalid("single VCF record exceeds planning budget"));
            }
            bytes += n;
            batch.push(line);
        }
        let selection = &input.selection;
        let masks = &input.masks;
        let parsed: Result<Vec<_>> = pool.install(|| {
            batch
                .par_iter()
                .enumerate()
                .map(|(index, line)| parse(line, rows + index as u64, selection, masks))
                .collect()
        });
        for parsed in parsed? {
            let row = parsed.row;
            let coordinate = (row.contig, row.position);
            if previous.is_some_and(|previous| coordinate < previous) {
                return Err(invalid("VCF contig/position order differs"));
            }
            previous = Some(coordinate);
            maximum = maximum.max(row.alleles.len() - 1);
            rows += 1;
            original_rows[row.contig] += 1;
            gt_hash.update(parsed.gt);
            record_hash.update(parsed.identity);
            *filters.get_mut(parsed.reason).expect("known exclusion") += 1;
            if let Some(value) = parsed.legacy {
                legacy_rows += 1;
                legacy_bytes += value.len() as u64;
                legacy_hash.update(value);
            }
            if parsed.reason == "retained" {
                let contig = row.contig;
                retained[contig] += 1;
                let window = input.masks[contig]
                    .windows
                    .partition_point(|window| window.stop < row.position);
                window_rows[contig][window] += 1;
                for (maximum, population) in max_an[contig].iter_mut().zip(&row.populations) {
                    *maximum = (*maximum).max(
                        population
                            .copies
                            .iter()
                            .map(|&value| u64::from(value))
                            .sum(),
                    );
                }
                serde_json::to_writer(&mut output, &row).map_err(json_error)?;
                output.write_all(b"\n")?;
            }
        }
    }
    output.flush()?;
    if rows == 0 || original_rows.contains(&0) || retained.contains(&0) {
        return Err(invalid(
            "complete workflow has an absent original/retained contig",
        ));
    }
    let gt_sha = format!("{:x}", gt_hash.finalize());
    let record_sha = format!("{:x}", record_hash.finalize());
    let mut files = BTreeMap::new();
    files.insert(
        "eligible.vcf".into(),
        FileIdentity {
            size: legacy_bytes,
            sha256: format!("{:x}", legacy_hash.finalize()),
        },
    );
    for pop in 0..2 {
        let value: String = input.selection.population_indices[pop]
            .iter()
            .map(|&index| format!("{}\n", input.selection.ids[index]))
            .collect();
        files.insert(
            format!("{}.txt", input.selection.names[pop]),
            FileIdentity {
                size: value.len() as u64,
                sha256: format!("{:x}", Sha256::digest(value.as_bytes())),
            },
        );
    }
    let legacy = LegacyIdentity {
        files,
        records: legacy_rows,
        input_gt_stream_sha256: gt_sha.clone(),
        original_record_profile: LEGACY_PROFILE,
        materialized: false,
    };
    Ok(Scan {
        input,
        spool,
        max_alt: maximum,
        rows,
        retained,
        original_rows,
        max_an,
        window_rows,
        filters,
        gt_sha,
        record_sha,
        legacy,
    })
}
