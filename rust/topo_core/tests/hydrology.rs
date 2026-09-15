//! 水文契约测试: 有界修正、嵌套河网与流向关联 HAND。

mod common;

use topo_core::hydro::{build_hydro, hand_to_stream, HydroConfig};

fn cfg() -> HydroConfig {
    HydroConfig {
        coarse_res_m: 25.0,
        z_limit_m: 15.0,
        stream_areas_km2: [0.05, 0.20, 1.00, 5.00],
    }
}

/// HAND 相对沿 flow_to 链遇到的第一个河网像元(计划契约示例)
#[test]
fn hand_follows_flow_connected_stream() {
    let dem = vec![110., 109., 108., 107., 106., 105., 104., 103., 102.];
    let flow_to = vec![1, 2, 3, 4, 5, 6, 7, 8, u32::MAX];
    let stream = vec![false, false, false, false, false, false, true, false, true];
    let hand = hand_to_stream(&dem, &flow_to, &stream, &[true; 9]).unwrap();
    assert_eq!(hand[0], 6.0);
}

/// HAND 使用链上首个河网而非欧氏最近河网: 像元 3 的链上首河网是
/// 5(高程 90), 而欧氏最近河网是 2(高程 108, 欧氏 HAND 会为负)
#[test]
fn hand_prefers_flow_connection_over_euclidean() {
    let dem = vec![110., 109., 108., 103., 100., 90.];
    let flow_to = vec![1, 2, 3, 4, 5, 5];
    let stream = vec![false, false, true, false, false, true];
    let hand = hand_to_stream(&dem, &flow_to, &stream, &[true; 6]).unwrap();
    assert_eq!(hand[0], 2.0, "链上首河网为像元 2");
    assert_eq!(hand[3], 13.0, "必须沿链到河网 5, 而非欧氏最近的河网 2");
}

/// 河网像元 HAND 恒为 0
#[test]
fn hand_is_zero_on_stream_cells() {
    let dem = vec![110., 109., 108., 107., 106., 105., 104., 103., 102.];
    let flow_to = vec![1, 2, 3, 4, 5, 6, 7, 8, u32::MAX];
    let stream = vec![false, false, false, false, false, false, true, false, true];
    let hand = hand_to_stream(&dem, &flow_to, &stream, &[true; 9]).unwrap();
    for (i, &s) in stream.iter().enumerate() {
        if s {
            assert_eq!(hand[i], 0.0, "stream cell {i}");
        }
    }
}

/// 下游无任何河网(直排出口)时 HAND 为 Float32 NoData(NaN)
#[test]
fn hand_is_nodata_without_stream_downstream() {
    let dem = vec![10., 9., 8.];
    let flow_to = vec![1, 2, 2];
    let stream = vec![false, false, false];
    let hand = hand_to_stream(&dem, &flow_to, &stream, &[true; 3]).unwrap();
    assert!(hand[0].is_nan() && hand[1].is_nan() && hand[2].is_nan());
}

/// 爬升路由边(深洼伪影)必须切段而非跨段锚定:
/// 链 5 -> 10 -> 8(stream) 中 0->1 为爬升边, 若跨段锚定将产生
/// hand[0] = 5-8 = -3 的负值; 正确行为是 0 为 NoData、1 正常锚定。
/// (Task 3 原负值报错保护的演进: 切段后结构性负值不可达,
/// -0.05 报错保留为内部不变量兜底)
#[test]
fn ascending_route_edges_are_segmented_not_negative() {
    let dem = vec![5., 10., 8.];
    let flow_to = vec![1, 2, 2];
    let stream = vec![false, false, true];
    let hand = hand_to_stream(&dem, &flow_to, &stream, &[true; 3]).unwrap();
    assert!(hand[0].is_nan(), "爬升边外侧像元不得跨段锚定: {}", hand[0]);
    assert_eq!(hand[1], 2.0);
    assert_eq!(hand[2], 0.0);
}

/// 每个路由像元沿 flow_to 最多 n 步到达出口; 汇流沿下游非降
#[test]
fn routing_reaches_outlet_and_accumulation_monotone() {
    let (dem, shape) = common::v_valley(201, 201, 10.0);
    let valid = vec![true; dem.len()];
    let prepared = topo_core::input::prepare_values(
        dem,
        common::meta_with_keys(201, 201, 10.0, 1, 9001),
        &Default::default(),
    )
    .unwrap();
    let model = build_hydro(&prepared, &cfg()).unwrap();
    let n = model.flow_to.len();
    for s in 0..n {
        let mut cur = s;
        let mut steps = 0usize;
        while model.flow_to[cur] != cur as u32 {
            cur = model.flow_to[cur] as usize;
            steps += 1;
            assert!(steps <= n, "flow path from {s} exceeds n steps");
        }
    }
    for i in 0..n {
        let t = model.flow_to[i] as usize;
        if t != i {
            assert!(
                model.accumulation_cells[t] >= model.accumulation_cells[i],
                "accumulation decreased downstream at {i}"
            );
        }
    }
}

/// 四级河网按物理面积嵌套: 粗等级河网是细等级子集, stream_level 一致
#[test]
fn streams_are_nested_by_physical_area() {
    let (dem, _) = common::v_valley(401, 401, 10.0);
    let prepared = topo_core::input::prepare_values(
        dem,
        common::meta_with_keys(401, 401, 10.0, 1, 9001),
        &Default::default(),
    )
    .unwrap();
    let model = build_hydro(&prepared, &cfg()).unwrap();
    for k in 0..3 {
        for i in 0..model.streams[k + 1].len() {
            assert!(
                !model.streams[k + 1][i] || model.streams[k][i],
                "level {} stream not nested inside level {} at {}",
                k + 1,
                k,
                i
            );
        }
    }
    let cell_area = cfg().coarse_res_m * cfg().coarse_res_m;
    for (k, &area) in cfg().stream_areas_km2.iter().enumerate() {
        let th = (area * 1_000_000.0 / cell_area).ceil() as u32;
        for i in 0..model.stream_level.len() {
            let expect =
                (model.accumulation_cells[i] >= th) as u8;
            assert_eq!(
                (model.stream_level[i] >= k as u8 + 1) as u8, expect,
                "stream_level mismatch at {i}"
            );
        }
    }
}

/// 修正深度不超过 z_limit; 深洼保持原始高程(内流区合法出口)
#[test]
fn conditioning_depth_bounded_and_deep_sinks_kept() {
    let (dem, _) = common::closed_pit(201, 201, 10.0);
    let prepared = topo_core::input::prepare_values(
        dem,
        common::meta_with_keys(201, 201, 10.0, 1, 9001),
        &Default::default(),
    )
    .unwrap();
    let model = build_hydro(&prepared, &cfg()).unwrap();
    for i in 0..model.conditioning_depth.len() {
        assert!(
            model.conditioning_depth[i] <= cfg().z_limit_m + 1e-3,
            "conditioning depth {} exceeds z-limit at {}",
            model.conditioning_depth[i],
            i
        );
    }
    assert!(
        model.deep_sink.iter().any(|&b| b),
        "闭合洼地中心应被记录为深洼"
    );
    // 河网像元 HAND 为 0; 上坡像元 HAND 为正
    let hand = hand_to_stream(
        &model.conditioned,
        &model.flow_to,
        &model.streams[0],
        &vec![true; model.conditioned.len()],
    )
    .unwrap();
    for (i, (&s, &v)) in model.streams[0].iter().zip(hand.iter()).enumerate() {
        if s {
            assert_eq!(v, 0.0, "HAND not zero on stream {i}");
        } else if v.is_finite() {
            assert!(v >= 0.0, "negative HAND at {i}: {v}");
        }
    }
}

/// 无效区不参与: 边缘 NoData 场景的粗层路由全部有效像元可达出口
#[test]
fn routing_skips_nodata_border() {
    let (dem, shape) = common::with_border_nodata(201, 201, 10.0);
    let prepared = topo_core::input::prepare_values(
        dem,
        common::meta_with_keys(shape.width as u32, shape.height as u32, 10.0, 1, 9001),
        &Default::default(),
    )
    .unwrap();
    let model = build_hydro(&prepared, &cfg()).unwrap();
    // 粗层路由森林无环且到达出口
    let n = model.flow_to.len();
    for s in 0..n {
        let mut cur = s;
        let mut steps = 0usize;
        while model.flow_to[cur] != cur as u32 {
            cur = model.flow_to[cur] as usize;
            steps += 1;
            assert!(steps <= n);
        }
    }
}
