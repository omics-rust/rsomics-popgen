use crate::model::CountRow;
use crate::numeric::{array, scalar, source_sum};
use serde_json::{Value, json};

pub(crate) fn wc(row: &CountRow, width: usize) -> [Vec<f64>; 3] {
    let numbers = row.populations.each_ref().map(|pop| {
        pop.copies
            .iter()
            .map(|&value| u64::from(value))
            .sum::<u64>()
            / 2
    });
    let total = (numbers[0] + numbers[1]) as f64;
    let nbar = total / 2.0;
    let nc = total
        - (numbers[0] as f64 * numbers[0] as f64 + numbers[1] as f64 * numbers[1] as f64) / total;
    let mut components = std::array::from_fn(|_| Vec::with_capacity(width));
    for allele in 0..width {
        let copies = row
            .populations
            .each_ref()
            .map(|pop| f64::from(pop.copies.get(allele).copied().unwrap_or(0)));
        let pbar = (copies[0] + copies[1]) / (2.0 * total);
        let p = [
            copies[0] / (2.0 * numbers[0] as f64),
            copies[1] / (2.0 * numbers[1] as f64),
        ];
        let variance = (numbers[0] as f64 * (p[0] - pbar).powi(2)
            + numbers[1] as f64 * (p[1] - pbar).powi(2))
            / nbar;
        let heterozygotes: u64 = row
            .populations
            .iter()
            .map(|pop| u64::from(pop.heterozygotes.get(allele).copied().unwrap_or(0)))
            .sum();
        let hbar = heterozygotes as f64 / total;
        let a = (nbar / nc)
            * (variance
                - (1.0 / (nbar - 1.0)) * (pbar * (1.0 - pbar) - variance / 2.0 - hbar / 4.0));
        let b = (nbar / (nbar - 1.0))
            * (pbar * (1.0 - pbar) - variance / 2.0 - ((2.0 * nbar - 1.0) * hbar / (4.0 * nbar)));
        components[0].push(a);
        components[1].push(b);
        components[2].push(hbar / 2.0);
    }
    components
}

pub(crate) fn jackknife(point: f64, numerator: &[f64], denominator: &[f64]) -> Value {
    let count = numerator.len();
    let vb: Vec<_> = numerator
        .iter()
        .zip(denominator)
        .map(|(&n, &d)| n / d)
        .collect();
    let vj: Vec<_> = (0..count)
        .map(|omit| {
            if count == 1 {
                return f64::NAN;
            }
            let n: Vec<_> = numerator
                .iter()
                .enumerate()
                .map(|(index, &value)| if index == omit { 0.0 } else { value })
                .collect();
            let d: Vec<_> = denominator
                .iter()
                .enumerate()
                .map(|(index, &value)| if index == omit { 0.0 } else { value })
                .collect();
            source_sum(&n, false) / source_sum(&d, false)
        })
        .collect();
    let mean = source_sum(&vj, false) / count as f64;
    let squared: Vec<_> = vj.iter().map(|value| (value - mean).powi(2)).collect();
    let se = (((count - 1) as f64 / count as f64) * source_sum(&squared, false)).sqrt();
    json!({"point": scalar(point), "se": scalar(se), "vb": array(&vb), "vj": array(&vj)})
}
