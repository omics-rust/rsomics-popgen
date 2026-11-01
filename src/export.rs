use crate::cli::Analyze;
use crate::input::identity;
use crate::legacy::Legacy;
use crate::model::{AccessibilityMask, CountRow, JACKKNIFE_PROFILE};
use crate::neutral::aggregate::{Aggregate, Totals};
use crate::neutral::diversity::Tajima;
use crate::neutral::sfs::Spectrum;
use crate::neutral::{self, SiteMetrics};
use crate::numeric::{array, general, scalar};
use crate::scan::Scan;
use crate::{invalid, json_error};
use flate2::{Compression, write::GzEncoder};
use rayon::prelude::*;
use rsomics_common::{AtomicFile, Result, write_atomic};
use serde_json::{Value, json};
use sha2::{Digest, Sha256};
use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

type Gzip = GzEncoder<BufWriter<File>>;

struct Contig {
    region: Aggregate,
    windows: Vec<Aggregate>,
    blocks: Vec<Aggregate>,
    spectrum: Spectrum,
    ancestry: BTreeMap<String, u64>,
    membership: Vec<[u64; 7]>,
    tail: Vec<u64>,
    observed: usize,
    writer: Gzip,
    transaction: AtomicFile,
}

fn gzip(path: &Path) -> Result<(Gzip, AtomicFile)> {
    let transaction = AtomicFile::new(path)?;
    let writer = GzEncoder::new(
        BufWriter::with_capacity(1024 * 1024, transaction.reopen()?),
        Compression::new(9),
    );
    Ok((writer, transaction))
}

pub(crate) fn save(path: &Path, value: &Value) -> Result<()> {
    write_atomic(path, |output| {
        serde_json::to_writer_pretty(&mut *output, value).map_err(json_error)?;
        output.write_all(b"\n")?;
        Ok(())
    })
}

fn site_text(row: &CountRow, value: &SiteMetrics, chrom: &str, width: usize) -> Result<String> {
    let mut fields = vec![
        row.ordinal.to_string(),
        chrom.to_owned(),
        row.position.to_string(),
        row.alleles[0].clone(),
        serde_json::to_string(&row.alleles[1..]).map_err(json_error)?,
        row.raw_aa.clone(),
        value.ancestry.0.to_string(),
        value.ancestry.1.into(),
    ];
    for pop in 0..2 {
        let mut counts = row.populations[pop].copies.clone();
        counts.resize(width, 0);
        let an: u64 = counts.iter().map(|&value| u64::from(value)).sum();
        fields.extend([
            an.to_string(),
            serde_json::to_string(&counts).map_err(json_error)?,
            general(value.within[pop], 17),
            if an < 2 { "True" } else { "False" }.into(),
        ]);
    }
    fields.extend([
        general(value.between, 17),
        if value.fixed { "True" } else { "False" }.into(),
        general(value.numerator, 17),
        general(value.denominator, 17),
    ]);
    for values in &value.wc {
        fields.push(serde_json::to_string(&array(values)).map_err(json_error)?);
    }
    let cells: Vec<_> = fields
        .into_iter()
        .map(|value| {
            if value.contains(['"', '\t', '\r', '\n']) {
                format!("\"{}\"", value.replace('"', "\"\""))
            } else {
                value
            }
        })
        .collect();
    Ok(cells.join("\t") + "\n")
}

fn window_text(
    mask: &AccessibilityMask,
    totals: &[Totals],
    names: &[String; 2],
    constants: &[Tajima; 2],
    counts: &[usize],
) -> String {
    let mut header: Vec<String> = ["start1", "stop1"].map(str::to_owned).into();
    for pop in names {
        for statistic in ["pi", "theta", "tajima_d"] {
            header.push(format!("{pop}_{statistic}"));
            if statistic != "tajima_d" {
                header.push(format!("{pop}_{statistic}_accessible_bp"));
            }
            header.push(format!("{pop}_{statistic}_records"));
        }
    }
    for statistic in ["dxy", "fixed_density", "hudson_fst", "wc_fst"] {
        header.push(statistic.into());
        if matches!(statistic, "dxy" | "fixed_density") {
            header.push(format!("{statistic}_accessible_bp"));
        }
        header.push(format!("{statistic}_records"));
    }
    header.extend(
        [
            "hudson_num_sum",
            "hudson_den_sum",
            "wc_a_sum",
            "wc_b_sum",
            "wc_c_sum",
        ]
        .map(str::to_owned),
    );
    let mut output = header.join("\t") + "\n";
    for (index, (window, total)) in mask.windows.iter().zip(totals).enumerate() {
        let accessible = mask.window_accessible[index];
        let bases = accessible as f64;
        let records = counts[index].to_string();
        let mut row = vec![window.start.to_string(), window.stop.to_string()];
        for (pop, constant) in constants.iter().enumerate() {
            row.extend([
                general(total.mpd[pop] / bases, 17),
                accessible.to_string(),
                records.clone(),
                general((total.segregating[pop] as f64 / constant.a1) / bases, 17),
                accessible.to_string(),
                records.clone(),
                general(constant.value(total.segregating[pop], total.mpd[pop]), 17),
                records.clone(),
            ]);
        }
        row.extend([
            general(total.dxy / bases, 17),
            accessible.to_string(),
            records.clone(),
            general(total.fixed as f64 / bases, 17),
            accessible.to_string(),
            records.clone(),
            general(total.hudson[0] / total.hudson[1], 17),
            records.clone(),
            general(
                total.wc[0] / ((total.wc[0] + total.wc[1]) + total.wc[2]),
                17,
            ),
            records,
        ]);
        row.extend(
            total
                .hudson
                .into_iter()
                .chain(total.wc)
                .map(|value| general(value, 17)),
        );
        output.push_str(&row.join("\t"));
        output.push('\n');
    }
    output
}

fn region(
    mask: &AccessibilityMask,
    totals: &Totals,
    names: &[String; 2],
    rows: usize,
    ancestry: &BTreeMap<String, u64>,
    constants: &[Tajima; 2],
) -> Value {
    let bases = mask.accessible_bases() as f64;
    let mut populations = serde_json::Map::new();
    for pop in 0..2 {
        populations.insert(names[pop].clone(), json!({"pi":scalar(totals.mpd[pop]/bases),
            "theta":scalar((totals.segregating[pop] as f64/constants[pop].a1)/bases),
            "tajima_d":scalar(constants[pop].value(totals.segregating[pop],totals.mpd[pop])),
            "segregating_sites":totals.segregating[pop],"n_max":totals.max_an[pop],"called_AN_min":totals.min_an[pop],
            "no_pair_sites":totals.no_pair[pop],"mpd_sum_fill0":scalar(totals.mpd[pop]),
            "theta_tajima_n_policy":"source n=max(called AN); no imputed invariant genotypes"}));
    }
    json!({"chrom":mask.name,"start":1,"stop":mask.length,"accessible_bp":mask.accessible_bases(),"records":rows,
        "populations":populations,"ancestral_status_counts":ancestry,"denominator":"project accessibility, not individual callability",
        "between":{"dxy":scalar(totals.dxy/bases),"fixed_differences":totals.fixed,"fixed_density":scalar(totals.fixed as f64/bases),
            "hudson_num_sum":scalar(totals.hudson[0]),"hudson_den_sum":scalar(totals.hudson[1]),
            "wc_a_sum":scalar(totals.wc[0]),"wc_b_sum":scalar(totals.wc[1]),"wc_c_sum":scalar(totals.wc[2])}})
}

fn scientific_identities(directory: &Path) -> Result<Value> {
    let mut values = serde_json::Map::new();
    for entry in fs::read_dir(directory)? {
        let path = entry?.path();
        if path.is_file() {
            values.insert(
                path.file_name()
                    .expect("file name")
                    .to_string_lossy()
                    .into(),
                serde_json::to_value(identity(&path)?).map_err(json_error)?,
            );
        }
    }
    Ok(Value::Object(values))
}

pub(crate) fn export(
    args: &Analyze,
    scan: &Scan,
    directory: &Path,
    pool: &rayon::ThreadPool,
    budget: usize,
) -> Result<Value> {
    let width = scan.max_alt + 1;
    let names = &scan.input.selection.names;
    let sizes = scan
        .input
        .selection
        .population_indices
        .each_ref()
        .map(Vec::len);
    let planned_aggregates: usize = scan.window_rows.iter().map(Vec::len).sum::<usize>()
        + scan
            .retained
            .iter()
            .map(|rows| rows / args.jackknife_variants)
            .sum::<usize>()
        + scan.retained.len();
    if planned_aggregates
        .checked_mul(12 * 1024)
        .is_none_or(|value| value > budget / 3)
    {
        return Err(invalid(
            "window/jackknife state exceeds the planning budget",
        ));
    }
    let mut contigs = Vec::new();
    for (index, mask) in scan.input.masks.iter().enumerate() {
        let count = scan.retained[index];
        let blocks = count / args.jackknife_variants;
        if blocks == 0 {
            return Err(invalid(
                "jackknife needs at least one full block; select a smaller --jackknife-variants",
            ));
        }
        let path = directory.join(&mask.name);
        fs::create_dir(&path)?;
        let (mut writer, transaction) = gzip(&path.join("sites.tsv.gz"))?;
        let mut header: Vec<_> = [
            "ordinal0",
            "chrom",
            "pos1",
            "ref",
            "alts",
            "raw_AA",
            "ancestral_allele_index",
            "ancestral_status",
        ]
        .map(str::to_owned)
        .into();
        for pop in names {
            for term in ["AN", "AC", "mpd", "no_pair"] {
                header.push(format!("{pop}_{term}"));
            }
        }
        header.extend(
            [
                "dxy_mpd",
                "fixed_difference",
                "hudson_num",
                "hudson_den",
                "wc_a",
                "wc_b",
                "wc_c",
            ]
            .map(str::to_owned),
        );
        writeln!(writer, "{}", header.join("\t"))?;
        contigs.push(Contig {
            region: Aggregate::new(count, width),
            windows: scan.window_rows[index]
                .iter()
                .map(|&rows| Aggregate::new(rows, width))
                .collect(),
            blocks: (0..blocks)
                .map(|_| Aggregate::new(args.jackknife_variants, width))
                .collect(),
            spectrum: Spectrum::new(names.clone(), sizes, budget)?,
            ancestry: BTreeMap::new(),
            membership: Vec::new(),
            tail: Vec::new(),
            observed: 0,
            writer,
            transaction,
        });
    }
    let mut legacy = if args.legacy_vcftools {
        Some(Legacy::new(directory, args.window_bp, sizes)?)
    } else {
        None
    };
    let mut reader = BufReader::new(scan.spool.reopen()?);
    let mut done = false;
    let mut ac_hash = Sha256::new();
    ac_hash.update(names[0].as_bytes());
    ac_hash.update([0]);
    let batch_limit = (budget / 16).min(8 * 1024 * 1024);
    while !done {
        let mut rows = Vec::new();
        let mut bytes = 0;
        while rows.len() < 8192 && bytes < batch_limit {
            let mut line = String::new();
            let n = reader.read_line(&mut line)?;
            if n == 0 {
                done = true;
                break;
            }
            bytes += n;
            rows.push(serde_json::from_str::<CountRow>(&line).map_err(json_error)?);
        }
        let masks = &scan.input.masks;
        let values: Result<Vec<_>> = pool.install(|| {
            rows.par_iter()
                .map(|row| {
                    let value = neutral::evaluate(row, width);
                    let text = site_text(row, &value, &masks[row.contig].name, width)?;
                    Ok((value, text))
                })
                .collect()
        });
        for (row, (value, text)) in rows.iter().zip(values?) {
            let mask = &scan.input.masks[row.contig];
            let contig = &mut contigs[row.contig];
            contig.writer.write_all(text.as_bytes())?;
            contig.region.push(row, &value);
            let window = mask
                .windows
                .partition_point(|window| window.stop < row.position);
            contig.windows[window].push(row, &value);
            contig.spectrum.push(row, value.ancestry.0);
            *contig.ancestry.entry(value.ancestry.1.into()).or_default() += 1;
            let block = contig.observed / args.jackknife_variants;
            if block < contig.blocks.len() {
                contig.blocks[block].push(row, &value);
                if contig.observed % args.jackknife_variants == 0 {
                    contig.membership.push([
                        block as u64,
                        contig.observed as u64,
                        (contig.observed + args.jackknife_variants) as u64,
                        row.ordinal,
                        row.ordinal,
                        row.position,
                        row.position,
                    ]);
                }
                let membership = &mut contig.membership[block];
                membership[4] = row.ordinal;
                membership[6] = row.position;
            } else {
                contig.tail.push(row.ordinal);
            }
            contig.observed += 1;
            if let Some(legacy) = &mut legacy {
                legacy.push(row, &mask.name)?;
            }
            for allele in 0..width {
                ac_hash.update(
                    i64::from(row.populations[0].copies.get(allele).copied().unwrap_or(0))
                        .to_le_bytes(),
                );
            }
        }
    }
    ac_hash.update(names[1].as_bytes());
    ac_hash.update([0]);
    for line in BufReader::new(scan.spool.reopen()?).lines() {
        let row: CountRow = serde_json::from_str(&line?).map_err(json_error)?;
        for allele in 0..width {
            ac_hash.update(
                i64::from(row.populations[1].copies.get(allele).copied().unwrap_or(0))
                    .to_le_bytes(),
            );
        }
    }
    if let Some(legacy) = legacy {
        legacy.finish()?;
    }
    let combined_ac_sha = format!("{:x}", ac_hash.finalize());
    let mut results = Vec::new();
    for (index, contig) in contigs.into_iter().enumerate() {
        let mask = &scan.input.masks[index];
        let path = directory.join(&mask.name);
        let mut output = contig.writer.finish()?;
        output.flush()?;
        drop(output);
        contig.transaction.commit()?;
        let constants = scan.max_an[index].map(Tajima::new);
        let totals = contig.region.finish();
        let windows: Vec<_> = contig.windows.into_iter().map(Aggregate::finish).collect();
        let (mut output, transaction) = gzip(&path.join("windows.tsv.gz"))?;
        output.write_all(
            window_text(mask, &windows, names, &constants, &scan.window_rows[index]).as_bytes(),
        )?;
        let mut output = output.finish()?;
        output.flush()?;
        drop(output);
        transaction.commit()?;
        let blocks: Vec<_> = contig.blocks.into_iter().map(Aggregate::finish).collect();
        let hudson_n: Vec<_> = blocks.iter().map(|block| block.hudson[0]).collect();
        let hudson_d: Vec<_> = blocks.iter().map(|block| block.hudson[1]).collect();
        let wc_n: Vec<_> = blocks.iter().map(|block| block.wc[0]).collect();
        let wc_d: Vec<_> = blocks.iter().map(|block| block.wc_den).collect();
        let jackknife_rows = blocks.len() * args.jackknife_variants;
        save(
            &path.join("fst-jackknife.json"),
            &json!({"profile":JACKKNIFE_PROFILE,"block_size":args.jackknife_variants,
            "point_rows":scan.retained[index],"jackknife_rows":jackknife_rows,"residual_tail_rows":contig.tail.len(),
            "residual_tail_index_range":[jackknife_rows,scan.retained[index]],"residual_tail_ordinals":contig.tail,
            "hudson":neutral::fst::jackknife(totals.hudson[0]/totals.hudson[1],&hudson_n,&hudson_d),
            "weir_cockerham":neutral::fst::jackknife(totals.wc[0]/((totals.wc[0]+totals.wc[1])+totals.wc[2]),&wc_n,&wc_d)}),
        )?;
        let (mut output, transaction) = gzip(&path.join("jackknife-membership.tsv.gz"))?;
        writeln!(
            output,
            "block0\tstart_index0\tstop_index_exclusive\tfirst_ordinal0\tlast_ordinal0\tfirst_pos1\tlast_pos1"
        )?;
        for row in contig.membership {
            writeln!(output, "{}", row.map(|value| value.to_string()).join("\t"))?;
        }
        let mut output = output.finish()?;
        output.flush()?;
        drop(output);
        transaction.commit()?;
        save(
            &path.join("region.json"),
            &region(
                mask,
                &totals,
                names,
                scan.retained[index],
                &contig.ancestry,
                &constants,
            ),
        )?;
        save(&path.join("sfs.json"), &contig.spectrum.finish()?)?;
        let ac_sha = if scan.input.masks.len() == 1 {
            combined_ac_sha.clone()
        } else {
            contig_ac_sha(scan, index, width)?
        };
        results.push(json!({"chrom":mask.name,"retained_rows":scan.retained[index],"window_rows":mask.windows.len(),"allele_counts_sha256":ac_sha,
            "allele_counts_identity":"population insertion order, UTF8 name NUL, full-width int64 little-endian variant/allele counts","output_files":scientific_identities(&path)?}));
    }
    Ok(
        json!({"contigs":results,"selected_gt_stream_sha256":scan.gt_sha,"record_stream_sha256":scan.record_sha}),
    )
}

fn contig_ac_sha(scan: &Scan, contig: usize, width: usize) -> Result<String> {
    let mut hash = Sha256::new();
    for pop in 0..2 {
        hash.update(scan.input.selection.names[pop].as_bytes());
        hash.update([0]);
        for line in BufReader::new(scan.spool.reopen()?).lines() {
            let row: CountRow = serde_json::from_str(&line?).map_err(json_error)?;
            if row.contig == contig {
                for allele in 0..width {
                    hash.update(
                        i64::from(
                            row.populations[pop]
                                .copies
                                .get(allele)
                                .copied()
                                .unwrap_or(0),
                        )
                        .to_le_bytes(),
                    );
                }
            }
        }
    }
    Ok(format!("{:x}", hash.finalize()))
}
