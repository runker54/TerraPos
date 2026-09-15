//! 端态契约锚(Task 1 冻结): 在 Task 12 暴露 run_arrays_for_test 前
//! 保持编译红, 用于证明测试连接的是新管线而非旧管线。

mod common;

use topo_core::pipeline::run_arrays_for_test;

/// 有效分类码域为 {1,3,4,5,6,7,8}: 代码 2 为保留码不得出现,
/// 全部有效像元必须非 0(0 仅用于 NoData)。
#[test]
fn final_codes_obey_contract() {
    let (dem, shape) = common::nested_ridges(257, 257, 10.0);
    let valid = vec![true; dem.len()];
    let out = run_arrays_for_test(&dem, &valid, shape, &Default::default()).unwrap();
    assert!(
        out.terrain
            .iter()
            .all(|c| matches!(c, 1 | 3 | 4 | 5 | 6 | 7 | 8)),
        "terrain codes outside contract"
    );
    assert!(!out.terrain.contains(&2), "reserved code 2 must never be emitted");
    assert!(
        out.terrain.iter().all(|&c| c != 0),
        "valid cells must all receive an effective class"
    );
}

/// 平面坡上每条脊→谷轨迹的坡位等级必须单调不增
/// (上(3/6)→中(4/7)→下(5/8), 允许缺级与同级重复, 禁止下→上反转)。
#[test]
fn planar_traces_never_invert() {
    let (dem, shape) = common::planar_slope(257, 257, 10.0);
    let valid = vec![true; dem.len()];
    let out = run_arrays_for_test(&dem, &valid, shape, &Default::default()).unwrap();
    let rank = |c: u8| match c {
        3 | 6 => 2u8,     // 上部
        4 | 7 => 1u8,     // 中部
        5 | 8 | 1 => 0u8, // 下部/盆地均视为谷端
        other => panic!("unexpected code {other} on planar trace"),
    };
    for col in (0..shape.width).step_by(8) {
        let mut prev = 2u8;
        for row in (0..shape.height).rev() {
            let i = row * shape.width + col;
            let r = rank(out.terrain[i]);
            assert!(
                r <= prev,
                "trace col {col} inverted at row {row}: {prev} -> {r}"
            );
            prev = r;
        }
    }
}
