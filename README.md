# rsomics-popgen

Private, unpublished Rust candidate for a connected neutral-variation and population-differentiation workflow. This is not a complete scikit-allel or vcftools replacement and has no release/performance qualification yet.

```text
rsomics-popgen analyze --vcf input.vcf.gz --panel samples.tsv \
  --population YRI --population CEU --reference-fai reference.fai \
  --mask 22=strict_mask.fa.gz --window-bp 100000 \
  --jackknife-variants 10000 --legacy-vcftools --out results --threads 4
```

The checked diploid VCF/sample/population model is parsed once. Selected calls are reduced to a private compact allele/heterozygote-count spool, then shared by diversity, SFS and Hudson/Weir–Cockerham estimators. Original sample order, global ALT width, raw AA, complete-case eligibility and input/count identities are preserved. Invalid or partial calls fail nonzero; they are not repaired silently.

Canonical reports follow the pinned scikit-allel1.3.13 neutral profile: PASS, project-P-accessible literal uppercase SNVs including multiallelic SNVs; explicit accessibility denominators; site/region/inclusive-window diversity, theta, Tajima D, Dxy and disjoint-allele fixed-difference summaries; complete-call biallelic folded/unfolded1D/joint2D spectra; FST components and variant-count block jackknife. Undefined values remain NaN/Inf and the residual jackknife tail is disclosed separately from all-row point estimates. At least one full jackknife block and nonempty complete/polarized SFS eligibility are required for this complete workflow.

`--legacy-vcftools` adds five native reports with the distinct vcftools0.1.17 pooled-cohort estimators, bins and denominators. They are not renamed canonical statistics. The exact would-be biallelic complete-case VCF and population-list bytes are hashed without materializing another cohort VCF. No Python, vcftools process or scientific FFI backend is invoked.

Each contig exports sites, windows, region, SFS, FST jackknife and block membership. Provenance, exact input identities, filtering, resolved profile and an output manifest accompany the reports. `--json` uses the rsomics-common result/error envelope; help uses rsomics-help. The output directory must not already exist. `--memory-mib` controls planning allocations, not a promised process RSS cap. Temporary work follows TMPDIR; project execution uses `scripts/guarded-cargo.sh` and a dedicated external target directory.

The repository includes only the synthetic engineering casebook and its actual pinned-source goldens, not human cohort data or an upstream runtime. One connected integration test exercises the real CLI at1/4 workers, complete reports, virtual derivative identities and malformed-GT failures. Native four-platform, representative repeated performance and release review remain separate gates.

Selection, admixture, PBS, LD, haplotype statistics, ploidy models other than the declared diploid profile, SFS projection and scaled spectra are not implemented or advertised. See THIRD_PARTY.md for source provenance and retained notices.
