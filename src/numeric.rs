use serde_json::Value;

pub(crate) fn scalar(value: f64) -> Value {
    if value.is_finite() {
        Value::from(value)
    } else {
        Value::from(if value.is_nan() {
            "NaN"
        } else if value.is_sign_negative() {
            "-Inf"
        } else {
            "+Inf"
        })
    }
}

pub(crate) fn array(values: &[f64]) -> Value {
    Value::Array(values.iter().copied().map(scalar).collect())
}

pub(crate) fn general(value: f64, significant: usize) -> String {
    if value.is_nan() {
        return if value.is_sign_negative() {
            "-nan"
        } else {
            "nan"
        }
        .into();
    }
    if value.is_infinite() {
        return if value.is_sign_negative() {
            "-inf"
        } else {
            "inf"
        }
        .into();
    }
    if value == 0.0 {
        return if value.is_sign_negative() { "-0" } else { "0" }.into();
    }
    let precision = significant - 1;
    let scientific = format!("{value:.precision$e}");
    let (mantissa, exponent) = scientific
        .split_once('e')
        .expect("Rust scientific formatter emits an exponent");
    let exponent: i32 = exponent
        .parse()
        .expect("Rust scientific exponent is an integer");
    if (-4..significant as i32).contains(&exponent) {
        let decimals = (precision as i32 - exponent).max(0) as usize;
        let mut text = format!("{value:.decimals$}");
        if text.contains('.') {
            while text.ends_with('0') {
                text.pop();
            }
            if text.ends_with('.') {
                text.pop();
            }
        }
        text
    } else {
        format!(
            "{}e{}{:02}",
            mantissa.trim_end_matches('0').trim_end_matches('.'),
            if exponent < 0 { '-' } else { '+' },
            exponent.abs()
        )
    }
}

// NumPy1.26.4's 128-element/eight-lane tree fixes reduction order independently of workers.
pub(crate) fn source_sum(values: &[f64], nan_sum: bool) -> f64 {
    let get = |index: usize| {
        if nan_sum && values[index].is_nan() {
            0.0
        } else {
            values[index]
        }
    };
    if values.len() < 8 {
        let mut sum = -0.0;
        for index in 0..values.len() {
            sum += get(index);
        }
        sum
    } else if values.len() <= 128 {
        let mut lanes = std::array::from_fn::<_, 8, _>(get);
        let stop = values.len() - values.len() % 8;
        for index in (8..stop).step_by(8) {
            for (lane, value) in lanes.iter_mut().enumerate() {
                *value += get(index + lane);
            }
        }
        let mut sum = ((lanes[0] + lanes[1]) + (lanes[2] + lanes[3]))
            + ((lanes[4] + lanes[5]) + (lanes[6] + lanes[7]));
        for index in stop..values.len() {
            sum += get(index);
        }
        sum
    } else {
        let half = (values.len() / 2) & !7;
        source_sum(&values[..half], nan_sum) + source_sum(&values[half..], nan_sum)
    }
}

pub(crate) struct StreamingSum {
    remaining: usize,
    leaf: Vec<f64>,
    leaf_length: usize,
    stack: Vec<(usize, Option<f64>)>,
    result: Option<f64>,
}

impl StreamingSum {
    pub fn new(length: usize) -> Self {
        let mut value = Self {
            remaining: length,
            leaf: Vec::with_capacity(128),
            leaf_length: 0,
            stack: Vec::new(),
            result: (length == 0).then_some(-0.0),
        };
        value.descend(length);
        value
    }

    fn descend(&mut self, mut length: usize) {
        while length > 128 {
            let half = (length / 2) & !7;
            self.stack.push((length - half, None));
            length = half;
        }
        self.leaf_length = length;
    }

    pub fn push(&mut self, value: f64) {
        assert!(self.remaining > 0);
        self.remaining -= 1;
        self.leaf.push(value);
        if self.leaf.len() == self.leaf_length {
            let mut sum = source_sum(&self.leaf, false);
            self.leaf.clear();
            while let Some((right_length, left)) = self.stack.pop() {
                if let Some(left) = left {
                    sum += left;
                } else {
                    self.stack.push((right_length, Some(sum)));
                    self.descend(right_length);
                    return;
                }
            }
            self.result = Some(sum);
        }
    }

    pub fn finish(self) -> f64 {
        assert_eq!(self.remaining, 0);
        self.result.expect("complete planned reduction")
    }
}
