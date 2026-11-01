pub(crate) mod aggregate;
pub(crate) mod diversity;
pub(crate) mod fst;
pub(crate) mod sfs;
pub(crate) mod windows;

use crate::model::CountRow;

pub(crate) struct SiteMetrics {
    pub within: [f64; 2],
    pub between: f64,
    pub fixed: bool,
    pub numerator: f64,
    pub denominator: f64,
    pub wc: [Vec<f64>; 3],
    pub ancestry: (i32, &'static str),
}

pub(crate) fn ancestral(row: &CountRow) -> (i32, &'static str) {
    if row.raw_aa.is_empty() || row.raw_aa == "." {
        return (-1, "missing");
    }
    let token = row
        .raw_aa
        .split('|')
        .next()
        .expect("split has a first token")
        .to_ascii_uppercase();
    if !matches!(token.as_str(), "A" | "C" | "G" | "T") {
        return (-1, "non_nucleotide");
    }
    let mut matches = row
        .alleles
        .iter()
        .enumerate()
        .filter(|(_, value)| value.eq_ignore_ascii_case(&token));
    match (matches.next(), matches.next()) {
        (Some((index, _)), None) => (index as i32, "matched"),
        (None, _) => (-1, "allele_mismatch"),
        _ => (-1, "ambiguous"),
    }
}

pub(crate) fn evaluate(row: &CountRow, width: usize) -> SiteMetrics {
    let within = row
        .populations
        .each_ref()
        .map(|value| diversity::within(&value.copies));
    let between = diversity::between(&row.populations[0].copies, &row.populations[1].copies);
    SiteMetrics {
        within,
        between,
        fixed: diversity::fixed(&row.populations[0].copies, &row.populations[1].copies),
        numerator: between - (within[0] + within[1]) / 2.0,
        denominator: between,
        wc: fst::wc(row, width),
        ancestry: ancestral(row),
    }
}
