use crate::model::CoordinateWindow;

pub(crate) fn position_windows(stop: u64, size: u64) -> Vec<CoordinateWindow> {
    let mut result = Vec::new();
    let mut start = 1;
    while start < stop {
        let end = start.saturating_add(size);
        if end >= stop {
            result.push(CoordinateWindow { start, stop });
            break;
        }
        result.push(CoordinateWindow {
            start,
            stop: end - 1,
        });
        start += size;
    }
    result
}
