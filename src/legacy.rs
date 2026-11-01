use crate::model::CountRow;
use crate::numeric::general;
use rsomics_common::{AtomicFile, Result, RsomicsError};
use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
#[cfg(unix)]
use std::os::unix::fs::{DirBuilderExt, PermissionsExt};
use std::path::Path;

struct Table {
    output: BufWriter<File>,
    transaction: AtomicFile,
}

impl Table {
    fn new(path: &Path, header: &str) -> Result<Self> {
        let mut transaction = AtomicFile::new(path)?;
        #[cfg(unix)]
        transaction
            .file_mut()
            .set_permissions(fs::Permissions::from_mode(0o600))
            .map_err(RsomicsError::Io)?;
        let mut output = BufWriter::with_capacity(1024 * 1024, transaction.reopen()?);
        writeln!(output, "{header}").map_err(RsomicsError::Io)?;
        Ok(Self {
            output,
            transaction,
        })
    }

    fn finish(mut self) -> Result<AtomicFile> {
        self.output.flush().map_err(RsomicsError::Io)?;
        drop(self.output);
        Ok(self.transaction)
    }
}

#[derive(Default)]
struct PiBin {
    variants: u64,
    differences: u128,
    comparisons: u128,
}

#[derive(Default)]
struct TajimaBin {
    segregating: u64,
    pq_sum: f64,
}

#[derive(Default)]
struct FstBin {
    numerator: f64,
    denominator: f64,
    ratio_sum: f64,
    count: u64,
}

struct ContigBins {
    name: String,
    pi: BTreeMap<u64, PiBin>,
    tajima: BTreeMap<u64, TajimaBin>,
    fst: BTreeMap<u64, FstBin>,
}

impl ContigBins {
    fn new(name: &str) -> Self {
        Self {
            name: name.to_owned(),
            pi: BTreeMap::new(),
            tajima: BTreeMap::new(),
            fst: BTreeMap::new(),
        }
    }
}

struct TajimaCoefficients {
    n: f64,
    a1: f64,
    e1: f64,
    e2: f64,
}

impl TajimaCoefficients {
    fn new(n: u64) -> Self {
        let mut a1 = 0.0;
        let mut a2 = 0.0;
        for index in 1..n {
            a1 += 1.0 / index as f64;
            a2 += 1.0 / (u128::from(index) * u128::from(index)) as f64;
        }
        let nf = n as f64;
        let b1 = (nf + 1.0) / 3.0 / (nf - 1.0);
        let square_term = u128::from(n) * u128::from(n) + u128::from(n) + 3;
        let b2 = 2.0 * square_term as f64 / 9.0 / nf / (nf - 1.0);
        let c1 = b1 - (1.0 / a1);
        let c2 = b2 - ((nf + 2.0) / (a1 * nf)) + (a2 / a1 / a1);
        Self {
            n: nf,
            a1,
            e1: c1 / a1,
            e2: c2 / ((a1 * a1) + a2),
        }
    }

    fn evaluate(&self, bin: &TajimaBin) -> f64 {
        if bin.segregating == 0 {
            return f64::NAN;
        }
        let s = bin.segregating as f64;
        let pi = 2.0 * bin.pq_sum * self.n / (self.n - 1.0);
        let theta = s / self.a1;
        let variance = (self.e1 * s) + self.e2 * s * (s - 1.0);
        (pi - theta) / variance.sqrt()
    }
}

pub(crate) struct Legacy {
    directory: std::path::PathBuf,
    window_bp: u64,
    population_sizes: [usize; 2],
    pooled_chromosomes: u64,
    coefficients: TajimaCoefficients,
    site_pi: Table,
    site_fst: Table,
    contigs: Vec<ContigBins>,
}

impl Legacy {
    pub fn new(directory: &Path, window_bp: u64, population_sizes: [usize; 2]) -> Result<Self> {
        if window_bp == 0 || population_sizes.contains(&0) {
            return Err(RsomicsError::InvalidInput(
                "legacy reports require a positive window and two nonempty populations".into(),
            ));
        }
        let directory = directory.join("vcftools");
        let mut builder = fs::DirBuilder::new();
        builder.recursive(true);
        #[cfg(unix)]
        builder.mode(0o700);
        builder.create(&directory).map_err(RsomicsError::Io)?;
        let pooled_chromosomes = population_sizes
            .iter()
            .try_fold(0_u64, |total, &size| total.checked_add(size as u64))
            .and_then(|total| total.checked_mul(2))
            .ok_or_else(|| RsomicsError::InvalidInput("population size overflow".into()))?;
        Ok(Self {
            site_pi: Table::new(&directory.join("site-pi.sites.pi"), "CHROM\tPOS\tPI")?,
            site_fst: Table::new(
                &directory.join("weir-fst.weir.fst"),
                "CHROM\tPOS\tWEIR_AND_COCKERHAM_FST",
            )?,
            directory,
            window_bp,
            population_sizes,
            pooled_chromosomes,
            coefficients: TajimaCoefficients::new(pooled_chromosomes),
            contigs: Vec::new(),
        })
    }

    pub fn push(&mut self, row: &CountRow, chrom: &str) -> Result<()> {
        if row.alleles.len() != 2
            || row
                .populations
                .iter()
                .zip(self.population_sizes)
                .any(|(population, size)| {
                    population
                        .copies
                        .iter()
                        .map(|&value| u64::from(value))
                        .sum::<u64>()
                        != 2 * size as u64
                })
        {
            return Ok(());
        }
        let copies = [
            u128::from(row.populations[0].copies[0]) + u128::from(row.populations[1].copies[0]),
            u128::from(row.populations[0].copies[1]) + u128::from(row.populations[1].copies[1]),
        ];
        let chromosomes = u128::from(self.pooled_chromosomes);
        let comparisons = chromosomes * (chromosomes - 1);
        let mut differences = 0;
        for count in copies {
            differences += count * (chromosomes - count);
        }
        writeln!(
            self.site_pi.output,
            "{chrom}\t{}\t{}",
            row.position,
            general(differences as f64 / comparisons as f64, 6)
        )
        .map_err(RsomicsError::Io)?;
        let (numerator, denominator) = weir_components(row);
        let ratio = numerator / denominator;
        writeln!(
            self.site_fst.output,
            "{chrom}\t{}\t{}",
            row.position,
            general(ratio, 6)
        )
        .map_err(RsomicsError::Io)?;

        if self
            .contigs
            .last()
            .is_none_or(|contig| contig.name != chrom)
        {
            self.contigs.push(ContigBins::new(chrom));
        }
        let contig = self
            .contigs
            .last_mut()
            .expect("the current contig was inserted");
        let tajima_index = (row.position as f64 * (1.0 / self.window_bp as f64)) as u64;
        let tajima = contig.tajima.entry(tajima_index).or_default();
        let ref_frequency = copies[0] as f64 / self.pooled_chromosomes as f64;
        if ref_frequency > 0.0 && ref_frequency < 1.0 {
            tajima.segregating += 1;
            tajima.pq_sum += ref_frequency * (1.0 - ref_frequency);
        }

        let first = ((row.position as f64 - self.window_bp as f64) / self.window_bp as f64)
            .ceil()
            .max(0.0) as u64;
        let last = (row.position as f64 / self.window_bp as f64).ceil() as u64;
        for index in first..last {
            if differences > 0 {
                let bin = contig.pi.entry(index).or_default();
                bin.variants += 1;
                bin.differences += differences;
                bin.comparisons += comparisons;
            }
            if !ratio.is_nan() {
                let bin = contig.fst.entry(index).or_default();
                bin.numerator += numerator;
                bin.denominator += denominator;
                bin.ratio_sum += ratio;
                bin.count += 1;
            }
        }
        Ok(())
    }

    pub fn finish(self) -> Result<()> {
        let mut window_pi = Table::new(
            &self.directory.join("window-pi.windowed.pi"),
            "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tN_MONOMORPHIC\tPI",
        )?;
        let mut tajima = Table::new(
            &self.directory.join("tajima-d.Tajima.D"),
            "CHROM\tBIN_START\tN_SNPS\tTajimaD",
        )?;
        let mut window_fst = Table::new(
            &self.directory.join("window-weir-fst.windowed.weir.fst"),
            "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tWEIGHTED_FST\tMEAN_FST",
        )?;
        let n = u128::from(self.pooled_chromosomes);
        let invariant_comparisons = n * (n - 1);
        for contig in &self.contigs {
            for (&index, bin) in &contig.pi {
                let monomorphic = self.window_bp.checked_sub(bin.variants).ok_or_else(|| {
                    RsomicsError::InvalidInput(format!(
                        "legacy window {}:{} has more variant rows than bases",
                        contig.name,
                        index * self.window_bp + 1
                    ))
                })?;
                let denominator = bin.comparisons + u128::from(monomorphic) * invariant_comparisons;
                writeln!(
                    window_pi.output,
                    "{}\t{}\t{}\t{}\t{}\t{}",
                    contig.name,
                    index * self.window_bp + 1,
                    index * self.window_bp + self.window_bp,
                    bin.variants,
                    monomorphic,
                    general(bin.differences as f64 / denominator as f64, 6)
                )
                .map_err(RsomicsError::Io)?;
            }
            if let Some((&first, _)) = contig.tajima.iter().find(|(_, bin)| bin.segregating > 0) {
                let last = *contig
                    .tajima
                    .last_key_value()
                    .expect("a segregating bin exists")
                    .0;
                let empty = TajimaBin::default();
                for index in first..=last {
                    let bin = contig.tajima.get(&index).unwrap_or(&empty);
                    writeln!(
                        tajima.output,
                        "{}\t{}\t{}\t{}",
                        contig.name,
                        index * self.window_bp,
                        bin.segregating,
                        general(self.coefficients.evaluate(bin), 6)
                    )
                    .map_err(RsomicsError::Io)?;
                }
            }
            for (&index, bin) in &contig.fst {
                if bin.denominator != 0.0
                    && !bin.numerator.is_nan()
                    && !bin.denominator.is_nan()
                    && bin.count > 0
                {
                    writeln!(
                        window_fst.output,
                        "{}\t{}\t{}\t{}\t{}\t{}",
                        contig.name,
                        index * self.window_bp + 1,
                        index * self.window_bp + self.window_bp,
                        general(bin.count as f64, 6),
                        general(bin.numerator / bin.denominator, 6),
                        general(bin.ratio_sum / bin.count as f64, 6)
                    )
                    .map_err(RsomicsError::Io)?;
                }
            }
        }
        AtomicFile::commit_all(vec![
            self.site_pi.finish()?,
            self.site_fst.finish()?,
            window_pi.finish()?,
            tajima.finish()?,
            window_fst.finish()?,
        ])
    }
}

fn weir_components(row: &CountRow) -> (f64, f64) {
    let mut n = [0.0; 2];
    let mut frequencies = [[0.0; 2]; 2];
    let mut pbar = [0.0; 2];
    let mut hbar = [0.0; 2];
    let mut squared_sizes = 0.0;
    for (population_index, population) in row.populations.iter().enumerate() {
        for (allele, &heterozygotes) in population.heterozygotes.iter().take(2).enumerate() {
            let heterozygotes = heterozygotes as f64;
            let homozygotes = (population.copies[allele] as f64 - heterozygotes) / 2.0;
            n[population_index] += homozygotes + 0.5 * heterozygotes;
            frequencies[population_index][allele] = heterozygotes + 2.0 * homozygotes;
            pbar[allele] += frequencies[population_index][allele];
            hbar[allele] += heterozygotes;
        }
        for frequency in &mut frequencies[population_index] {
            *frequency /= 2.0 * n[population_index];
        }
        squared_sizes += n[population_index] * n[population_index];
    }
    let n_sum = n[0] + n[1];
    let nbar = n_sum / 2.0;
    let nc = (n_sum - (squared_sizes / n_sum)) / (2.0 - 1.0);
    let mut numerator = 0.0;
    let mut denominator = 0.0;
    for (allele, (pbar, hbar)) in pbar.iter_mut().zip(&mut hbar).enumerate() {
        *pbar /= n_sum * 2.0;
        *hbar /= n_sum;
        let mut variance = 0.0;
        for (size, frequencies) in n.iter().zip(&frequencies) {
            let difference = frequencies[allele] - *pbar;
            variance += *size * difference * difference;
        }
        variance /= (2.0 - 1.0) * nbar;
        let pq = *pbar * (1.0 - *pbar);
        let a = (variance - (pq - (((2.0 - 1.0) * variance) / 2.0) - (*hbar / 4.0)) / (nbar - 1.0))
            * nbar
            / nc;
        let b =
            (pq - (variance * (2.0 - 1.0) / 2.0) - *hbar * (((2.0 * nbar) - 1.0) / (4.0 * nbar)))
                * nbar
                / (nbar - 1.0);
        let c = *hbar / 2.0;
        if !a.is_nan() && !b.is_nan() && !c.is_nan() {
            numerator += a;
            denominator += a + b + c;
        }
    }
    (numerator, denominator)
}
