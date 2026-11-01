#![forbid(unsafe_code)]

mod cli;
mod export;
mod input;
mod legacy;
mod model;
mod neutral;
mod numeric;
mod scan;

use cli::{Analyze, Cli, Command};
use rsomics_common::{Result, RsomicsError, ToolMeta};
use serde_json::{Value, json};
use std::collections::BTreeMap;
use std::fs;
use std::path::Path;
use std::time::Instant;

pub(crate) fn invalid(message: impl Into<String>) -> RsomicsError {
    RsomicsError::InvalidInput(message.into())
}
pub(crate) fn json_error(error: serde_json::Error) -> RsomicsError {
    invalid(error.to_string())
}

pub fn run_cli() -> std::process::ExitCode {
    let cli: Cli = rsomics_help::parse();
    rsomics_common::run(
        &cli.output,
        ToolMeta {
            name: env!("CARGO_PKG_NAME"),
            version: env!("CARGO_PKG_VERSION"),
        },
        || {
            let threads = cli.threads.requested().map_or_else(
                || std::thread::available_parallelism().map(usize::from),
                |value| Ok(value.get()),
            )?;
            let budget = cli
                .memory_mib
                .checked_mul(1024 * 1024)
                .ok_or_else(|| invalid("memory budget overflows"))?;
            if budget < 32 * 1024 * 1024 {
                return Err(invalid("--memory-mib must be at least 32"));
            }
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .map_err(|error| RsomicsError::ConfigError(error.to_string()))?;
            let Command::Analyze(args) = cli.command;
            analyze(&args, &pool, budget, cli.memory_mib)
        },
    )
}

fn analyze(
    args: &Analyze,
    pool: &rayon::ThreadPool,
    budget: usize,
    memory_mib: usize,
) -> Result<Value> {
    let started = Instant::now();
    if args.out.exists() || args.out.is_symlink() {
        return Err(invalid("output directory must not exist"));
    }
    let input = input::prepare(args, budget)?;
    let sizes = input.selection.population_indices.each_ref().map(Vec::len);
    let cells = (sizes[0] * 2 + 1)
        .checked_mul(sizes[1] * 2 + 1)
        .ok_or_else(|| invalid("SFS axes overflow"))?;
    let fixed_state = input
        .masks
        .len()
        .checked_mul(
            2 * 1024 * 1024
                + cells
                    .checked_mul(16)
                    .ok_or_else(|| invalid("SFS size overflows"))?,
        )
        .ok_or_else(|| invalid("contig state overflows"))?;
    if fixed_state > budget / 3 {
        return Err(invalid("contig exports/SFS exceed the planning budget"));
    }
    let scan = scan::scan(input, budget, pool)?;
    fs::create_dir(&args.out)?;
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        fs::set_permissions(&args.out, fs::Permissions::from_mode(0o700))?;
    }
    let result = export::export(args, &scan, &args.out, pool, budget)?;
    let names = &scan.input.selection.names;
    let mut contigs = serde_json::Map::new();
    for (index, mask) in scan.input.masks.iter().enumerate() {
        contigs.insert(mask.name.clone(),json!({"length":mask.length,"rows":scan.original_rows[index],"accessible_bp":mask.accessible_bases(),"mask_codes":mask.codes}));
    }
    let population_indices = BTreeMap::from([
        (
            names[0].clone(),
            scan.input.selection.population_indices[0].clone(),
        ),
        (
            names[1].clone(),
            scan.input.selection.population_indices[1].clone(),
        ),
    ]);
    export::save(
        &args.out.join("input-admission.json"),
        &json!({"schema":1,"original_samples":scan.input.selection.original_samples,
        "sample_ids":scan.input.selection.ids,"sample_indices":scan.input.selection.original_columns,"population_order":names,
        "population_indices":population_indices,"max_alt_number":scan.max_alt,"gt_shape":[scan.rows,scan.input.selection.ids.len() as u64,2],
        "filter_counts":scan.filters,"selected_gt_stream_sha256":scan.gt_sha,"record_stream_sha256":scan.record_sha,
        "contigs":contigs,"files":scan.input.files,"legacy":scan.legacy}),
    )?;
    export::save(
        &args.out.join("filter-counts.json"),
        &serde_json::to_value(&scan.filters).map_err(json_error)?,
    )?;
    export::save(&args.out.join("result.json"), &result)?;
    export::save(
        &args.out.join("run-profile.json"),
        &json!({"population_order":names,"window_bp":args.window_bp,"jackknife_variants":args.jackknife_variants,
        "legacy_vcftools":args.legacy_vcftools,"canonical":"scikit-allel-1.3.13-neutral","gt":"complete-diploid-or-fully-missing",
        "sites":"PASS-P-accessible-literal-uppercase-SNV","aa":"unique-nucleotide-before-pipe-no-fallback"}),
    )?;
    export::save(
        &args.out.join("legacy-provenance.json"),
        &json!({"backend":"native-rust","upstream":"vcftools-0.1.17","derivative":scan.legacy,
        "operations":if args.legacy_vcftools { vec!["site-pi","window-pi","tajima-d","weir-fst","window-weir-fst"] } else { Vec::<&str>::new() }}),
    )?;
    export::save(
        &args.out.join("request.json"),
        &json!({"argv":std::env::args_os().map(|value|value.to_string_lossy().into_owned()).collect::<Vec<_>>(),"cwd":std::env::current_dir()?.to_string_lossy()}),
    )?;
    export::save(
        &args.out.join("runtime.json"),
        &json!({"implementation":"native-rust","product":"rsomics-popgen","version":env!("CARGO_PKG_VERSION"),
        "target":std::env::consts::ARCH,"system":std::env::consts::OS,"threads":pool.current_num_threads(),"memory_mib":memory_mib,
        "planning_budget_is_not_process_rss_limit":true,"elapsed_seconds":started.elapsed().as_secs_f64(),"count_spool_bytes":scan.spool.as_file().metadata()?.len(),
        "executable":input::identity(&std::env::current_exe()?)?}),
    )?;
    let mut manifest = serde_json::Map::new();
    collect_manifest(&args.out, &args.out, &mut manifest)?;
    export::save(
        &args.out.join("output-manifest.json"),
        &Value::Object(manifest),
    )?;
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        fs::set_permissions(
            args.out.join("output-manifest.json"),
            fs::Permissions::from_mode(0o600),
        )?;
    }
    if !std::env::args_os().any(|value| value == "--json") {
        eprintln!(
            "Analyzed {} original records; retained {} neutral sites in {}",
            scan.rows,
            scan.retained.iter().sum::<usize>(),
            args.out.display()
        );
    }
    Ok(
        json!({"output":args.out,"original_records":scan.rows,"retained_records":scan.retained.iter().sum::<usize>(),"contigs":scan.input.masks.len(),"threads":pool.current_num_threads()}),
    )
}

fn collect_manifest(
    root: &Path,
    directory: &Path,
    output: &mut serde_json::Map<String, Value>,
) -> Result<()> {
    for entry in fs::read_dir(directory)? {
        let path = entry?.path();
        if path.is_symlink() {
            return Err(invalid("output tree unexpectedly contains a link"));
        }
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            fs::set_permissions(
                &path,
                fs::Permissions::from_mode(if path.is_dir() { 0o700 } else { 0o600 }),
            )?;
        }
        if path.is_dir() {
            collect_manifest(root, &path, output)?;
        } else {
            output.insert(
                path.strip_prefix(root)
                    .expect("output descendant")
                    .to_string_lossy()
                    .into_owned(),
                serde_json::to_value(input::identity(&path)?).map_err(json_error)?,
            );
        }
    }
    Ok(())
}
