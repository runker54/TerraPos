//! 输入契约测试: 投影米制校验、NoData 解析、边界无效区保留与内部小孔填补。

mod common;

use common::meta_with_keys;
use topo_core::geotiff::{read_f32, write_f32, GeoMeta};
use topo_core::input::{validate_meta, InputConfig, prepare_values};

fn cfg() -> InputConfig {
    InputConfig {
        max_interior_hole_area_m2: 10_000.0,
        geomorph_smooth_radius_m: 5.0,
    }
}

fn flat_values(n: usize) -> Vec<f32> {
    vec![900.0; n]
}

#[test]
fn accepts_float32_projected_metre() {
    let meta = meta_with_keys(32, 32, 5.0, 1, 9001);
    let shape = validate_meta(&meta).unwrap();
    assert_eq!(shape.width, 32);
    assert_eq!(shape.height, 32);
    assert_eq!(shape.resolution_m, 5.0);
}

#[test]
fn rejects_geographic_crs_even_when_pixel_size_is_small() {
    let meta = meta_with_keys(32, 32, 0.00005, 2, 9102);
    let err = validate_meta(&meta).unwrap_err().to_string();
    assert!(err.contains("投影坐标系") && err.contains("米"), "{err}");
}

#[test]
fn rejects_projected_feet() {
    let meta = meta_with_keys(32, 32, 16.4, 1, 9003);
    let err = validate_meta(&meta).unwrap_err().to_string();
    assert!(err.contains("9003"), "错误应包含观测到的线性单位键: {err}");
}

#[test]
fn rejects_missing_geokeys() {
    let mut meta = meta_with_keys(32, 32, 5.0, 1, 9001);
    meta.geo_keys.clear();
    let err = validate_meta(&meta).unwrap_err().to_string();
    assert!(err.contains("投影坐标系"), "{err}");
}

#[test]
fn rejects_non_square_or_non_positive_pixels() {
    let mut meta = meta_with_keys(32, 32, 5.0, 1, 9001);
    meta.pixel_scale = [5.0, 6.0, 0.0];
    assert!(validate_meta(&meta).is_err(), "X/Y 分辨率差超 1% 应拒绝");

    meta.pixel_scale = [0.0, 0.0, 0.0];
    assert!(validate_meta(&meta).is_err(), "非正分辨率应拒绝");
}

#[test]
fn rejects_raster_smaller_than_three_cells() {
    let meta = meta_with_keys(2, 32, 5.0, 1, 9001);
    assert!(validate_meta(&meta).is_err());
}

#[test]
fn parses_gdal_nodata_tag_on_roundtrip() {
    let dir = std::env::temp_dir().join("topo_core_input_test");
    std::fs::create_dir_all(&dir).unwrap();
    let p = dir.join("nodata_tag.tif");
    let mut meta = meta_with_keys(16, 16, 5.0, 1, 9001);
    meta.nodata = Some(-9999.0);
    let data: Vec<f32> = (0..256).map(|i| 800.0 + i as f32).collect();
    write_f32(&p, &meta, &data).unwrap();
    let (_, back) = read_f32(&p).unwrap();
    assert_eq!(back.nodata, Some(-9999.0));
}

#[test]
fn boundary_connected_nodata_is_preserved() {
    let (dem, _) = common::with_border_nodata(64, 64, 5.0);
    let meta = meta_with_keys(64, 64, 5.0, 1, 9001);
    let prepared = prepare_values(dem, meta, &cfg()).unwrap();
    let border_invalid = prepared
        .valid
        .iter()
        .enumerate()
        .filter(|(i, _)| {
            let x = i % 64;
            let y = i / 64;
            x < 8 || y < 8 || x >= 56 || y >= 56
        })
        .all(|(_, &v)| !v);
    assert!(border_invalid, "边界连通无效区必须保持无效");
    let inner_valid = (8..56).all(|y| (8..56).all(|x| prepared.valid[y * 64 + x]));
    assert!(inner_valid, "内部有效区不得被侵蚀");
}

#[test]
fn small_interior_hole_is_filled() {
    let (mut dem, _) = common::conical_hill(64, 64, 5.0);
    // 中心 10x10=100 像元 = 2500 m² 孔(低于 10000 m² 上限)
    for y in 27..37 {
        for x in 27..37 {
            dem[y * 64 + x] = f32::NAN;
        }
    }
    let meta = meta_with_keys(64, 64, 5.0, 1, 9001);
    let prepared = prepare_values(dem, meta, &cfg()).unwrap();
    assert!(
        prepared.valid.iter().all(|&v| v),
        "不超过面积上限的内部小孔应被填补为有效"
    );
    let filled = (27..37).all(|y| (27..37).all(|x| prepared.raw[y * 64 + x].is_finite()));
    assert!(filled, "填补后 raw 必须为有限值");
}

#[test]
fn large_interior_hole_is_preserved_invalid() {
    let (mut dem, _) = common::conical_hill(64, 64, 5.0);
    // 中心 30x30=900 像元 = 22500 m² 孔(超上限, 保持无效)
    for y in 17..47 {
        for x in 17..47 {
            dem[y * 64 + x] = f32::NAN;
        }
    }
    let meta = meta_with_keys(64, 64, 5.0, 1, 9001);
    let prepared = prepare_values(dem, meta, &cfg()).unwrap();
    let hole_invalid = (17..47).all(|y| (17..47).all(|x| !prepared.valid[y * 64 + x]));
    let outside_valid = (2..15).all(|y| (2..15).all(|x| prepared.valid[y * 64 + x]));
    assert!(hole_invalid, "超上限内部孔必须保持无效");
    assert!(outside_valid, "孔外有效区不受影响");
}

#[test]
fn sentinel_nodata_values_are_excluded() {
    let mut values = flat_values(64 * 64);
    values[0] = -9999.0;
    let mut meta = meta_with_keys(64, 64, 5.0, 1, 9001);
    meta.nodata = Some(-9999.0);
    let prepared = prepare_values(values, meta, &cfg()).unwrap();
    assert!(!prepared.valid[0], "NoData 哨兵像元必须无效");
    assert!(prepared.valid[1], "正常像元必须有效");
}

#[test]
fn rejects_insufficient_valid_ratio() {
    let mut values = flat_values(32 * 32);
    for v in values.iter_mut().take(32 * 32 - 10) {
        *v = f32::NAN;
    }
    let meta = meta_with_keys(32, 32, 5.0, 1, 9001);
    let err = prepare_values(values, meta, &cfg()).unwrap_err().to_string();
    assert!(err.contains("有效"), "{err}");
}

#[test]
fn geomorph_surface_is_smoothed_valid_aware() {
    // 平面坡内部平滑不应改变高程(平面均值=自身); 边缘缺邻域的均值偏移
    // 是 valid-aware 滤波的已知边缘特性, 不在断言范围
    let (dem, _) = common::planar_slope(32, 32, 5.0);
    let meta = meta_with_keys(32, 32, 5.0, 1, 9001);
    let prepared = prepare_values(dem.clone(), meta, &cfg()).unwrap();
    for y in 1..31 {
        for x in 1..31 {
            let i = y * 32 + x;
            assert!(
                (prepared.geomorph[i] - dem[i]).abs() < 1e-3,
                "平面坡内部平滑前后应一致: {} vs {}",
                prepared.geomorph[i],
                dem[i]
            );
        }
    }
}
