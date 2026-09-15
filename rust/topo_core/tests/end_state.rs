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

/// 从脊像元沿距离到谷的梯度下降到谷锚的轨迹, 坡位等级必须单调不增
/// (上(3/6)→中(4/7)→下(5/8), 允许缺级与同级重复, 禁止下→上反转);
/// 与计划 14.2 同口径: >=95% 轨迹满足单调
#[test]
fn planar_traces_never_invert() {
    let (dem, shape) = common::nested_ridges(257, 257, 10.0);
    let valid = vec![true; dem.len()];
    let out = run_arrays_for_test(&dem, &valid, shape, &Default::default()).unwrap();
    let rank = |c: u8| match c {
        3 | 6 => 2i32,
        4 | 7 => 1,
        5 | 8 | 1 => 0,
        _ => -1,
    };
    let w = shape.width;
    // 确定性采样 500 条种子(脊顶区: 编码 3/6 像元), 沿编码"梯度"贪心下降
    let mut ok = 0usize;
    let mut total = 0usize;
    for s0 in (0..dem.len()).step_by(dem.len() / 500) {
        let start_rank = rank(out.terrain[s0]);
        if start_rank < 2 {
            continue;
        }
        total += 1;
        let x = (s0 % w) as i64;
        let y = (s0 / w) as i64;
        // 下降方向: 8 邻中高程最低(谷向), 用高程而非内部 dv(测试无诊断层)
        let mut cur = s0;
        let mut trace = vec![cur];
        let mut prev_rank = start_rank;
        let mut good = true;
        loop {
            let cx = (cur % w) as i64;
            let cy = (cur / w) as i64;
            let mut next = None;
            let mut best_z = f32::INFINITY;
            for (dx, dy) in [(-1i64, 0i64), (1, 0), (0, -1), (0, 1), (-1, -1), (1, -1), (-1, 1), (1, 1)] {
                let nx = cx + dx;
                let ny = cy + dy;
                if nx < 0 || ny < 0 || nx >= w as i64 || ny >= shape.height as i64 {
                    continue;
                }
                let j = (ny as usize) * w + nx as usize;
                if !trace.contains(&j) {
                    continue;
                }
                let _ = (x, y);
            }
            // 高程最低的未访问邻像元, 且必须严格低于当前像元
            // (单一坡面下降段; 到谷底即终止, 不横穿谷底爬对坡)
            let cur_z = dem[cur];
            for (dx, dy) in [(-1i64, 0i64), (1, 0), (0, -1), (0, 1), (-1, -1), (1, -1), (-1, 1), (1, 1)] {
                let nx = cx + dx;
                let ny = cy + dy;
                if nx < 0 || ny < 0 || nx >= w as i64 || ny >= shape.height as i64 {
                    continue;
                }
                let j = (ny as usize) * w + nx as usize;
                if trace.contains(&j) || dem[j] >= cur_z {
                    continue;
                }
                if dem[j] < best_z {
                    best_z = dem[j];
                    next = Some(j);
                }
            }
            match next {
                Some(j) => {
                    trace.push(j);
                    cur = j;
                    let r = rank(out.terrain[j]);
                    if r >= 0 {
                        if r > prev_rank {
                            good = false;
                            break;
                        }
                        prev_rank = r;
                    }
                }
                None => break,
            }
            if trace.len() >= 120 {
                break;
            }
        }
        if good {
            ok += 1;
        }
    }
    let rate = if total > 0 { ok as f64 / total as f64 } else { 0.0 };
    assert!(total >= 100, "脊顶种子不足: {total}");
    assert!(rate >= 0.95, "梯度轨迹单调率 {rate:.3} < 95%");
}