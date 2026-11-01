use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub(crate) struct PopulationCounts {
    pub copies: Vec<u32>,
    pub heterozygotes: Vec<u32>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub(crate) struct CountRow {
    pub ordinal: u64,
    pub contig: usize,
    pub position: u64,
    pub alleles: Vec<String>,
    pub raw_aa: String,
    pub populations: [PopulationCounts; 2],
}

#[derive(Clone, Copy, Debug)]
pub(crate) enum DiploidCall {
    Missing,
    Called(u16, u16),
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct CoordinateWindow {
    pub start: u64,
    pub stop: u64,
}

#[derive(Clone, Debug)]
pub(crate) struct Selection {
    pub names: [String; 2],
    pub ids: Vec<String>,
    pub original_columns: Vec<usize>,
    pub population_of_original: Vec<Option<usize>>,
    pub population_indices: [Vec<usize>; 2],
    pub original_samples: usize,
}

#[derive(Clone, Debug)]
pub(crate) struct AccessibilityMask {
    pub name: String,
    pub length: u64,
    pub bits: Vec<u64>,
    pub codes: BTreeMap<String, u64>,
    pub windows: Vec<CoordinateWindow>,
    pub window_accessible: Vec<u64>,
}

impl AccessibilityMask {
    pub fn accessible(&self, position: u64) -> bool {
        let index = (position - 1) as usize;
        self.bits[index / 64] & (1 << (index % 64)) != 0
    }

    pub fn accessible_bases(&self) -> u64 {
        self.codes.get("P").copied().unwrap_or(0)
    }
}

#[derive(Clone, Debug, Serialize)]
pub(crate) struct FileIdentity {
    pub size: u64,
    pub sha256: String,
}

#[derive(Clone, Debug, Serialize)]
pub(crate) struct LegacyIdentity {
    pub files: BTreeMap<String, FileIdentity>,
    pub records: u64,
    pub input_gt_stream_sha256: String,
    pub original_record_profile: &'static str,
    pub materialized: bool,
}

pub(crate) const LEGACY_PROFILE: &str = "PASS P-accessible literal biallelic SNV; both populations complete; original fields and GT preserved";
pub(crate) const SFS_PROFILE: &str = "supplied eligible VCF rows only; no projection; no invented invariant bin0; independently folded marginal minor counts";
pub(crate) const JACKKNIFE_PROFILE: &str = "source variant-count block jackknife; all-row point estimate; incomplete residual tail excluded from resampling; not demographic confidence interval";
