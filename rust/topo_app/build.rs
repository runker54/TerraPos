fn main() {
    // 构建时间戳(UTC+8), 注入 UI 顶栏用于版本辨识
    let now = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_secs() as i64
        + 8 * 3600;
    let days = now.div_euclid(86400);
    let rem = now.rem_euclid(86400);
    let z = days + 719468;
    let era = z.div_euclid(146097);
    let doe = z.rem_euclid(146097);
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let mut y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = (doy - (153 * mp + 2) / 5 + 1) as u32;
    let m = (if mp < 10 { mp + 3 } else { mp - 9 }) as u32;
    if m <= 2 {
        y += 1;
    }
    println!(
        "cargo:rustc-env=BUILD_STAMP={:04}-{:02}-{:02} {:02}:{:02}",
        y,
        m,
        d,
        rem / 3600,
        (rem % 3600) / 60
    );

    if std::env::var("CARGO_CFG_TARGET_OS").as_deref() == Ok("windows") {
        let mut res = winresource::WindowsResource::new();
        res.set_icon("assets/icon.ico");
        res.compile().expect("compile windows resources");
    }
}
