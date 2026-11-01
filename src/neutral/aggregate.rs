use crate::model::CountRow;
use crate::neutral::SiteMetrics;
use crate::numeric::StreamingSum;

pub(crate) struct Aggregate {
    scalar: [StreamingSum; 5],
    wc: [StreamingSum; 3],
    wc_den: StreamingSum,
    segregating: [u64; 2],
    min_an: [u64; 2],
    max_an: [u64; 2],
    no_pair: [u64; 2],
    fixed: u64,
}

pub(crate) struct Totals {
    pub mpd: [f64; 2],
    pub dxy: f64,
    pub hudson: [f64; 2],
    pub wc: [f64; 3],
    pub wc_den: f64,
    pub segregating: [u64; 2],
    pub min_an: [u64; 2],
    pub max_an: [u64; 2],
    pub no_pair: [u64; 2],
    pub fixed: u64,
}

fn nan_zero(value: f64) -> f64 {
    if value.is_nan() { 0.0 } else { value }
}

impl Aggregate {
    pub fn new(rows: usize, width: usize) -> Self {
        Self {
            scalar: std::array::from_fn(|_| StreamingSum::new(rows)),
            wc: std::array::from_fn(|_| StreamingSum::new(rows * width)),
            wc_den: StreamingSum::new(rows * width),
            segregating: [0; 2],
            min_an: [u64::MAX; 2],
            max_an: [0; 2],
            no_pair: [0; 2],
            fixed: 0,
        }
    }

    pub fn push(&mut self, row: &CountRow, value: &SiteMetrics) {
        for pop in 0..2 {
            self.scalar[pop].push(nan_zero(value.within[pop]));
            let ac = &row.populations[pop].copies;
            let an: u64 = ac.iter().map(|&value| u64::from(value)).sum();
            self.segregating[pop] += u64::from(ac.iter().filter(|&&value| value > 0).count() > 1);
            self.min_an[pop] = self.min_an[pop].min(an);
            self.max_an[pop] = self.max_an[pop].max(an);
            self.no_pair[pop] += u64::from(an < 2);
        }
        self.scalar[2].push(nan_zero(value.between));
        self.scalar[3].push(nan_zero(value.numerator));
        self.scalar[4].push(nan_zero(value.denominator));
        for allele in 0..value.wc[0].len() {
            for component in 0..3 {
                self.wc[component].push(nan_zero(value.wc[component][allele]));
            }
            self.wc_den.push(nan_zero(
                (value.wc[0][allele] + value.wc[1][allele]) + value.wc[2][allele],
            ));
        }
        self.fixed += u64::from(value.fixed);
    }

    pub fn finish(self) -> Totals {
        let scalar = self.scalar.map(StreamingSum::finish);
        Totals {
            mpd: [scalar[0], scalar[1]],
            dxy: scalar[2],
            hudson: [scalar[3], scalar[4]],
            wc: self.wc.map(StreamingSum::finish),
            wc_den: self.wc_den.finish(),
            segregating: self.segregating,
            min_an: self.min_an,
            max_an: self.max_an,
            no_pair: self.no_pair,
            fixed: self.fixed,
        }
    }
}
