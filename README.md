# rsomics-popgen

One Rust binary for a connected, two-population neutral-variation and
population-differentiation workflow. It reads selected diploid calls once,
shares compact allele counts across diversity, SFS and FST, and exports complete
reports with explicit filtering, accessibility and estimator provenance.

The implementation has accepted four-native-platform and full-chr22
correctness/performance evidence. This delivery revision still needs its own
package/native CI; publication and release assets are not claimed here.
It is not a complete scikit-allel or vcftools replacement.

## Install from source

Rust1.91 or newer is required. No Python, vcftools process or scientific FFI
runtime is used by the product.

```sh
git clone https://github.com/omics-rust/rsomics-popgen.git
cd rsomics-popgen
cargo install --locked --path .
rsomics-popgen --help
```

Owning-repository CI tests Linux and macOS natively on x86_64 and aarch64.
It prepares one native archive per platform and the source crate, with
SHA256SUMS and license notices, without publishing to a registry. A CI artifact
is a candidate delivery, not automatically a signed or published release.

## Analyze

```text
rsomics-popgen analyze --vcf input.vcf.gz --panel samples.tsv \
  --population YRI --population CEU --reference-fai reference.fai \
  --mask 22=strict_mask.fa.gz --window-bp 100000 \
  --jackknife-variants 10000 --legacy-vcftools --out results --threads 4
```

Supply exactly two distinct, nonempty populations. The VCF, panel, FAI and
strict-accessibility masks must describe the same declared contigs. See the
[input and output guide](docs/neutral-analysis.md) for schemas, supported inputs,
nonfinite values and the distinction between canonical and legacy estimates.
Invalid or partially missing calls fail nonzero; they are not repaired silently.

Canonical reports follow the pinned scikit-allel1.3.13 neutral profile: PASS, project-P-accessible literal uppercase SNVs including multiallelic SNVs; explicit accessibility denominators; site/region/inclusive-window diversity, theta, Tajima D, Dxy and disjoint-allele fixed-difference summaries; complete-call biallelic folded/unfolded1D/joint2D spectra; FST components and variant-count block jackknife. Undefined values remain NaN/Inf and the residual jackknife tail is disclosed separately from all-row point estimates. At least one full jackknife block and nonempty complete/polarized SFS eligibility are required for this complete workflow.

`--legacy-vcftools` adds five native reports with the distinct vcftools0.1.17 pooled-cohort estimators, bins and denominators. They are not renamed canonical statistics. The exact would-be biallelic complete-case VCF and population-list bytes are hashed without materializing another cohort VCF. No Python, vcftools process or scientific FFI backend is invoked.

Each contig exports sites, windows, region, SFS, FST jackknife and block
membership. Provenance, exact input identities, filtering, resolved profile and
an output manifest accompany the reports. `--json` uses the rsomics-common
result/error envelope; help uses rsomics-help. The output directory must not
already exist. `--memory-mib` controls planning allocations, not a process RSS
cap. Temporary work follows TMPDIR.

The source package includes a small synthetic casebook and pinned-source
goldens, not human cohort data or an upstream runtime. From the source root:

```sh
rsomics-popgen analyze \
  --vcf tests/fixtures/neutral/casebook.vcf \
  --panel tests/fixtures/neutral/panel.tsv \
  --population YRI --population CEU \
  --reference-fai tests/fixtures/neutral/reference.fai \
  --mask a=tests/fixtures/neutral/a.mask.fa \
  --mask b=tests/fixtures/neutral/b.mask.fa \
  --window-bp 10 --jackknife-variants 3 --legacy-vcftools \
  --threads 4 --out casebook-results
```

## Evidence and limits

The accepted implementation is
[`ff8e9513`](https://github.com/omics-rust/rsomics-world/blob/main/.autopilot/oracles/popgen-native-2026-11-01.json).
Its source-bound casebook passes debug/release tests on all four native target
classes in [run41265985145](https://github.com/omics-rust/rsomics-world/actions/runs/41265985145).
The later
[full-workflow acceptance](https://github.com/omics-rust/rsomics-world/blob/main/.autopilot/oracles/popgen-whole-real-2026-11-01.json)
adopts [run41272268312](https://github.com/omics-rust/rsomics-world/actions/runs/41272268312):
all1103547 chr22 VCF rows,207 selected persons,723783 canonical SNVs,514 windows,
full spectra/jackknife outputs and all five legacy reports. Canonical finite
comparison uses the fixed `max(2^-44,2^-40*abs(source))` rule; discrete and
nonfinite classes are exact. Legacy parsed decimals match their source profile.

On that Linux x86_64 AMD EPYC7763 host, one warmup plus five alternating pairs
at each CPU allowance measured the complete parse-to-output workflow:

| CPU allowance | Source median wall | Native median wall | Source process RSS median | Native process RSS median |
|---|---:|---:|---:|---:|
| 1 | 404.609s | 69.103s | 2272808960B | 34942976B |
| 4 | 383.532s | 45.453s | 2272841728B | 36777984B |

Source was CPython3.12.10/scikit-allel1.3.13/NumPy1.26.4/Dask2024.9.0 plus five
actual vcftools0.1.17 commands. Native avoids materializing their derivative
VCF while preserving exact virtual input identities and report semantics;
this is a whole-workflow benefit, not an isolated-kernel speed comparison.
RSS is process lifetime high-water; source RSS is the largest waited process,
not simultaneous aggregate memory. Cold cache and exclusive-host execution
are not claimed. Native1→4 wall improves1.520x while total CPU rises1.651x;
single-worker system CPU is a retained loss,5.571s versus source4.592s.
Results belong to the exact historical executable and profile, not every
machine or a newly built delivery binary without fresh qualification.

Selection, admixture, PBS, LD, haplotype statistics, ploidy models other than the declared diploid profile, SFS projection and scaled spectra are not implemented or advertised. See THIRD_PARTY.md for source provenance and retained notices.
