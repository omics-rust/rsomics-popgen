# Neutral analysis input and output contract

`rsomics-popgen analyze` computes the connected neutral workflow for exactly
two distinct, nonempty populations. The canonical profile is
scikit-allel1.3.13-neutral; `--legacy-vcftools` adds five different native
vcftools0.1.17 report estimators. It does not invoke either upstream runtime.

## Inputs

VCF, panel and masks may be plain text or gzip files with a `.gz` suffix;
multi-member gzip is accepted. Input paths name regular readable files, not
stdin. All supplied input contents are recorded by size and SHA256.

| Argument | Required contract |
|---|---|
| `--vcf` | Ordered VCF with a standard sample header, `GT` FORMAT declaration `Number=1,Type=String`, and `AA` INFO declaration `Number=1,Type=String`. Selected calls are diploid or fully missing; partial missing, malformed/out-of-range alleles and unsorted coordinates fail. Phased/unphased separators are both accepted for this unphased workflow. |
| `--panel` | Whitespace-delimited header starting `sample pop`; additional columns are allowed and every row has the header width. Sample identities are unique. Every sample belonging to either selected population must occur in the VCF. Selected sample order follows the VCF, not the panel. |
| `--population` | Repeat exactly twice with distinct population labels. Each selected population must contain at least one VCF sample. Other panel populations are not selected. |
| `--reference-fai` | Five-column FAI with unique contig names and lengths of at least two bases. Supplied mask contigs are resolved in its order. |
| `--mask CONTIG=PATH` | Repeat for the declared contigs. A single-contig FASTA has exactly the FAI length; its header is `CONTIG` or `chrCONTIG`. Allowed uppercase codes are `P N L H Z Q 0`; only `P` is accessible. This is project accessibility, not individual callability. |
| `--index` | Optional original index file identity only. It neither crops input nor enables random access. |
| `--window-bp` | Positive inclusive-coordinate window width, default100000. |
| `--jackknife-variants` | Positive retained-variant block size, default10000. This is not a base-pair window. Each analyzed contig needs at least one complete block. |
| `--out` | A directory that must not already exist. No overwrite/append mode is offered. |

The complete declared contig profile needs retained sites, complete-call
biallelic SFS eligibility and matched-ancestry eligibility for unfolded spectra.
An unsupported empty required population fails rather than inventing a spectrum
or uncertainty result. Use an appropriate smaller positive jackknife size for
small valid data; the synthetic example uses3.

Canonical filtering is exclusive and ordered: not-PASS, not-literal-uppercase-
SNV, inaccessible, retained. Multiallelic literal SNVs are kept. Global ALT
width and the original record/selected-GT identities include excluded records.
AA uses the uppercase token before the first `|`, requires A/C/G/T and a unique
REF/ALT match. Missing, non-nucleotide, mismatching and ambiguous ancestry remain
distinct; no REF fallback is applied.

`--threads` selects native workers; omitted means available parallelism.
Scientific output order and arithmetic do not change with the worker count.
`--memory-mib` defaults512 and must be at least32. It constrains planning
allocations, not total process RSS; insufficient planning capacity is an error.
Private count-spool scratch follows TMPDIR and is not an exported genotype file.

## Outputs

Each declared contig has a directory containing:

| File | Meaning |
|---|---|
| `sites.tsv.gz` | Original ordinal/coordinates/alleles/AA classification, per-population allele copies and AN, MPD, Dxy, fixed-difference indicator, Hudson numerator/denominator and per-allele WC components. |
| `windows.tsv.gz` | Inclusive `start1`/`stop1` windows with statistic-specific record counts, accessible denominators and aggregated Hudson/WC components. Empty/zero-accessible windows remain present. |
| `region.json` | Full declared-contig summaries and population/accessibility/ancestry accounting. |
| `sfs.json` | Folded/unfolded1D and joint2D spectra, explicit dimensions, eligibility and ancestry accounting. |
| `fst-jackknife.json` | All-row point estimates, complete variant-count blocks, leave-one-block values, uncertainty and disclosed residual tail. |
| `jackknife-membership.tsv.gz` | Retained-row/original-ordinal membership in complete blocks or the exact excluded residual tail. |

Top-level files include `input-admission.json`, `filter-counts.json`,
`run-profile.json`, `result.json`, `request.json`, `runtime.json`,
`legacy-provenance.json` and `output-manifest.json`. They record selected
samples/order/populations, input digests, resolved policies, filtering,
runtime/budget/executable identity and every exported file's identity.

`output-manifest.json` is written last and does not hash itself. A nonzero exit
may leave partial reports, but never treat that directory as a completed run.
The manifest is not a substitute for checking the process exit code. Preserve
failed output for diagnosis; choose a new output directory for another run.

Canonical TSV cells use general17-significant-digit finite values. Array cells
use JSON text and may be quoted as TSV fields. Canonical JSON represents
nonfinite numbers as strings `"NaN"`, `"+Inf"`, `"-Inf"`; TSV uses numeric-style
`nan`/`inf` spellings. Do not coerce these to zero/null. NaN payload/sign and
zero sign are not cross-ISA scientific contracts; infinity sign and nonfinite
class are. Negative defined FST is not clamped.

## Estimator distinctions

Canonical pi/theta/Dxy densities use explicit P-accessible denominators.
Site MPD/Dxy can be undefined while their density aggregation applies the
upstream's separate fill-zero rule. Fixed difference means both populations
are called and share no supported allele; each need not be monomorphic.
Theta/Tajima use the maximum called AN for that complete contig call, not
an independently chosen per-window sample size. Tajima requires at least three
segregating sites. Inclusive final windows can exceed the nominal width by one
base under the upstream endpoint rule.

Hudson and Weir–Cockerham FST aggregate components; they are not averages of
per-site FST ratios. The all-row point estimate includes the jackknife residual
tail, whereas uncertainty uses only complete equal-size variant blocks. One
complete block can legitimately produce undefined leave-one-block/SE values.
The tail is disclosed, not silently discarded from the point estimate.

SFS uses complete-call biallelic eligibility and fixed twice-population-size
axes. Joint folding folds each marginal minor count independently. No
projection, scaling or imputed invariant-genome bin is offered.

With `--legacy-vcftools`, `vcftools/` additionally contains
`site-pi.sites.pi`, `window-pi.windowed.pi`, `tajima-d.Tajima.D`,
`weir-fst.weir.fst` and `window-weir-fst.windowed.weir.fst`. Their complete-case
biallelic derivative eligibility, pooled pi/Tajima cohort, bins, full-window
denominators and finite default6-significant-digit report values are distinct
from canonical outputs. These files must not be relabelled canonical estimates.
The exact would-be derivative VCF/population-list bytes are hashed with
`materialized=false`; the product does not claim those large files exist.

Help/version/completions use rsomics-help; machine result/errors use
rsomics-common via `--json`. That envelope is a command summary, not a second
copy of the scientific tables or selected human genotype data.
