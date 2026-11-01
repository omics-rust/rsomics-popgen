use crate::numeric::source_sum;

pub(crate) fn within(counts: &[u32]) -> f64 {
    let number: u64 = counts.iter().map(|&value| u64::from(value)).sum();
    if number < 2 {
        return f64::NAN;
    }
    let pairs = (number * (number - 1)) as f64 / 2.0;
    let same: Vec<f64> = counts
        .iter()
        .map(|&value| f64::from(value) * (f64::from(value) - 1.0) / 2.0)
        .collect();
    (pairs - source_sum(&same, false)) / pairs
}

pub(crate) fn between(first: &[u32], second: &[u32]) -> f64 {
    let n1: u64 = first.iter().map(|&value| u64::from(value)).sum();
    let n2: u64 = second.iter().map(|&value| u64::from(value)).sum();
    let pairs = n1 * n2;
    if pairs == 0 {
        return f64::NAN;
    }
    let same: u64 = first
        .iter()
        .zip(second)
        .map(|(&a, &b)| u64::from(a) * u64::from(b))
        .sum();
    (pairs - same) as f64 / pairs as f64
}

pub(crate) fn fixed(first: &[u32], second: &[u32]) -> bool {
    first.iter().any(|&value| value > 0)
        && second.iter().any(|&value| value > 0)
        && first.iter().zip(second).all(|(&a, &b)| a == 0 || b == 0)
}

#[derive(Clone, Copy)]
pub(crate) struct Tajima {
    pub a1: f64,
    e1: f64,
    e2: f64,
}

impl Tajima {
    pub fn new(number: u64) -> Self {
        let terms: Vec<_> = (1..number).map(|value| 1.0 / value as f64).collect();
        let a1 = source_sum(&terms, false);
        let terms: Vec<_> = (1..number)
            .map(|value| 1.0 / (value * value) as f64)
            .collect();
        let a2 = source_sum(&terms, false);
        let n = number as f64;
        let b1 = (n + 1.0) / (3.0 * (n - 1.0));
        let b2 = 2.0 * (n * n + n + 3.0) / (9.0 * n * (n - 1.0));
        let c1 = b1 - 1.0 / a1;
        let c2 = b2 - (n + 2.0) / (a1 * n) + a2 / (a1 * a1);
        Self {
            a1,
            e1: c1 / a1,
            e2: c2 / (a1 * a1 + a2),
        }
    }

    pub fn value(self, segregating: u64, mpd_sum: f64) -> f64 {
        if segregating < 3 {
            return f64::NAN;
        }
        let s = segregating as f64;
        (mpd_sum - s / self.a1) / (self.e1 * s + self.e2 * s * (s - 1.0)).sqrt()
    }
}
