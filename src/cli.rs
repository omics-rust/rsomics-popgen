use clap::{Args, Parser, Subcommand};
use std::path::PathBuf;

#[derive(Parser)]
#[command(
    name = "rsomics-popgen",
    version,
    about = "Connected neutral variation and population differentiation"
)]
pub(crate) struct Cli {
    #[command(flatten)]
    pub output: rsomics_common::OutputArgs,
    #[command(flatten)]
    pub threads: rsomics_common::ThreadArgs,
    #[arg(long, global = true, default_value_t = 512)]
    pub memory_mib: usize,
    #[command(subcommand)]
    pub command: Command,
}

#[derive(Subcommand)]
pub(crate) enum Command {
    /// Analyze diversity, spectra and population differentiation together.
    Analyze(Analyze),
}

#[derive(Args)]
pub(crate) struct Analyze {
    #[arg(long)]
    pub vcf: PathBuf,
    #[arg(long)]
    pub panel: PathBuf,
    #[arg(long, num_args = 1, action = clap::ArgAction::Append, required = true)]
    pub population: Vec<String>,
    #[arg(long)]
    pub reference_fai: PathBuf,
    /// Strict project-accessibility FASTA, supplied as CONTIG=PATH.
    #[arg(long, required = true)]
    pub mask: Vec<String>,
    /// Retain an original VCF index identity without cropping the contig.
    #[arg(long)]
    pub index: Option<PathBuf>,
    #[arg(long, default_value_t = 100000)]
    pub window_bp: u64,
    #[arg(long, default_value_t = 10000)]
    pub jackknife_variants: usize,
    /// Also emit the distinct vcftools0.1.17 legacy estimator reports.
    #[arg(long)]
    pub legacy_vcftools: bool,
    #[arg(long)]
    pub out: PathBuf,
}
