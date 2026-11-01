use flate2::read::MultiGzDecoder;
use serde_json::Value;
use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::Command;

fn fixture() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/neutral")
}

fn invoke(input: &Path, out: &Path, threads: usize) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_rsomics-popgen"))
        .arg("analyze")
        .arg("--vcf")
        .arg(input)
        .arg("--panel")
        .arg(fixture().join("panel.tsv"))
        .args(["--population", "YRI", "--population", "CEU"])
        .arg("--reference-fai")
        .arg(fixture().join("reference.fai"))
        .arg("--mask")
        .arg(format!("a={}", fixture().join("a.mask.fa").display()))
        .arg("--mask")
        .arg(format!("b={}", fixture().join("b.mask.fa").display()))
        .args([
            "--window-bp",
            "10",
            "--jackknife-variants",
            "3",
            "--legacy-vcftools",
        ])
        .arg("--threads")
        .arg(threads.to_string())
        .arg("--out")
        .arg(out)
        .output()
        .unwrap()
}

fn numeric(got: f64, expected: f64, legacy: bool) {
    if expected.is_nan() {
        assert!(got.is_nan());
        return;
    }
    if expected.is_infinite() {
        assert_eq!(got, expected);
        return;
    }
    assert!(got.is_finite());
    let bound = if legacy {
        0.0
    } else {
        (2.0_f64.powi(-44)).max(2.0_f64.powi(-40) * expected.abs())
    };
    assert!(
        (got - expected).abs() <= bound,
        "{got} != {expected}, bound {bound}"
    );
}

fn json_value(got: &Value, expected: &Value) {
    match expected {
        Value::Object(map) => {
            let actual = got.as_object().unwrap();
            assert_eq!(actual.len(), map.len());
            for (key, expected) in map {
                json_value(&actual[key], expected);
            }
        }
        Value::Array(values) => {
            let actual = got.as_array().unwrap();
            assert_eq!(actual.len(), values.len());
            for (a, b) in actual.iter().zip(values) {
                json_value(a, b);
            }
        }
        Value::Number(value) if value.is_f64() => {
            numeric(got.as_f64().unwrap(), value.as_f64().unwrap(), false)
        }
        _ => assert_eq!(got, expected),
    }
}

fn text(path: &Path) -> String {
    if path.extension().is_some_and(|value| value == "gz") {
        let mut value = String::new();
        MultiGzDecoder::new(fs::File::open(path).unwrap())
            .read_to_string(&mut value)
            .unwrap();
        value
    } else {
        fs::read_to_string(path).unwrap()
    }
}

fn table(got: &str, expected: &str, name: &str) {
    let mut a = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(got.as_bytes());
    let mut b = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(expected.as_bytes());
    let headers = b.headers().unwrap().clone();
    assert_eq!(a.headers().unwrap(), &headers);
    let mut a = a.records();
    for row in b.records() {
        let actual = a.next().unwrap().unwrap();
        let expected = row.unwrap();
        let left: Vec<_> = actual.iter().collect();
        let right: Vec<_> = expected.iter().collect();
        assert_eq!(left.len(), right.len());
        for ((field, got), expected) in headers.iter().zip(left).zip(right) {
            if got.starts_with('[') {
                json_value(
                    &serde_json::from_str::<Value>(got).unwrap(),
                    &serde_json::from_str::<Value>(expected).unwrap(),
                );
            } else {
                let float = if name.starts_with("vcftools/") {
                    matches!(
                        field,
                        "PI" | "TajimaD" | "WEIR_AND_COCKERHAM_FST" | "WEIGHTED_FST" | "MEAN_FST"
                    )
                } else if name.ends_with("sites.tsv.gz") {
                    field.ends_with("_mpd") || matches!(field, "hudson_num" | "hudson_den")
                } else if name.ends_with("windows.tsv.gz") {
                    !matches!(field, "start1" | "stop1")
                        && !field.ends_with("_records")
                        && !field.ends_with("_accessible_bp")
                } else {
                    false
                };
                if float {
                    numeric(
                        got.parse().unwrap(),
                        expected.parse().unwrap(),
                        name.starts_with("vcftools/"),
                    );
                } else {
                    assert_eq!(got, expected, "{name}/{field}");
                }
            }
        }
    }
    assert!(a.next().is_none());
}

#[test]
fn neutral_workflow_complete() {
    let scratch = tempfile::tempdir_in(std::env::temp_dir()).unwrap();
    let one = scratch.path().join("one");
    let four = scratch.path().join("four");
    for (threads, out) in [(1, &one), (4, &four)] {
        let run = invoke(&fixture().join("casebook.vcf"), out, threads);
        assert!(
            run.status.success(),
            "{}",
            String::from_utf8_lossy(&run.stderr)
        );
    }
    let names: Vec<String> =
        serde_json::from_slice(&fs::read(fixture().join("scientific-files.json")).unwrap())
            .unwrap();
    for name in names {
        let expected = text(&fixture().join("expected").join(&name));
        let first = text(&one.join(&name));
        let second = text(&four.join(&name));
        assert_eq!(first, second, "worker payload {name}");
        if name.ends_with(".json") {
            json_value(
                &serde_json::from_str::<Value>(&first).unwrap(),
                &serde_json::from_str::<Value>(&expected).unwrap(),
            );
        } else {
            table(&first, &expected, &name);
        }
    }
    let admission: Value =
        serde_json::from_slice(&fs::read(one.join("input-admission.json")).unwrap()).unwrap();
    let expected: Value =
        serde_json::from_slice(&fs::read(fixture().join("expected-admission.json")).unwrap())
            .unwrap();
    for key in [
        "gt_shape",
        "max_alt_number",
        "filter_counts",
        "population_indices",
        "selected_gt_stream_sha256",
        "record_stream_sha256",
    ] {
        assert_eq!(admission[key], expected[key], "admission {key}");
    }
    assert_eq!(admission["legacy"]["files"], expected["legacy"]["files"]);
    assert_eq!(admission["legacy"]["materialized"], false);
    for (name, value) in [("invalid", "99/0"), ("partial", "0/.")] {
        let source = fs::read_to_string(fixture().join("casebook.vcf")).unwrap();
        let altered = scratch.path().join(format!("{name}.vcf"));
        fs::write(&altered, source.replacen("0/0", value, 1)).unwrap();
        let out = scratch.path().join(name);
        let run = invoke(&altered, &out, 4);
        assert!(!run.status.success());
        assert!(!out.join("output-manifest.json").exists());
    }
}
