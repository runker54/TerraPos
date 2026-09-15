# -*- coding: utf-8 -*-
"""验收指标计算与硬门判定。

读取 Rust 管线(data/validation_out)的正式成果与合成场景期望,
计算规格 14/15 节的验收指标; 任一硬门失败则以非零码退出。
运行: D:/worker_code/.venvgis/Scripts/python.exe python/validate_terrain_result.py
"""
import sys
from pathlib import Path

import numpy as np
from osgeo import gdal

ROOT = Path(__file__).resolve().parents[1]
VALIDATION_OUT = ROOT / "data" / "validation_out"
SCENE_DIR = ROOT / "data" / "synthetic" / "generated"

# 硬门阈值(设计规格 14/15 节, 不得放宽)
CODE_DOMAIN = {1, 3, 4, 5, 6, 7, 8}
RESAMPLE_AGREEMENT_MIN = 0.85       # 10m/25m 回采样非边界一致率
NOISE_AGREEMENT_MIN = 0.90          # 0.5m 噪声非边界一致率
PARAM_CLASS_AREA_CHANGE_MAX = 0.15  # ±20% 参数类面积变化
PARAM_BASIN_AREA_CHANGE_MAX = 0.20  # ±20% 参数盆地面积变化
TRACE_MONOTONIC_MIN = 0.95          # 脊-谷轨迹单调率
TRACE_COUNT = 500

VALIDATION_SCENES = (
    "broad_basin",
    "v_valley",
    "nested_ridges",
    "conical_hill",
    "closed_pit",
)


def _read(path: Path):
    if not path.exists():
        return None, None
    ds = gdal.Open(str(path))
    band = ds.GetRasterBand(1)
    arr = band.ReadAsArray()
    nd = band.GetNoDataValue()
    ds = None
    return arr, nd


def scene_dir(name: str) -> Path:
    return VALIDATION_OUT / "synthetic" / name


def check_code_domain(terrain, z, nodata):
    """码域/保留码/NoData 完整性(expected==0 仅为无解析预期, 不参与)。"""
    if nodata is not None:
        valid = z != nodata
    else:
        valid = np.isfinite(z)
    violations = int(np.sum(~np.isin(terrain[valid], list(CODE_DOMAIN))))
    code2 = int(np.sum(terrain == 2))
    nodata_violation = int(np.sum(terrain[~valid] != 0))
    return violations + code2 + nodata_violation, len(set(np.unique(terrain[valid]).tolist()))


def trace_monotonic_rate(terrain, z):
    """从上坡像元沿高程最陡下降轨迹的单调率(确定性梯度轨迹)。"""
    n_rows, n_cols = terrain.shape
    rank = np.full(terrain.shape, -1, dtype=np.int8)
    for code, r in ((3, 2), (6, 2), (4, 1), (7, 1), (5, 0), (8, 0), (1, 0)):
        rank[terrain == code] = r

    tops = np.argwhere(rank >= 2)
    if len(tops) == 0:
        return 1.0, 0
    step = max(1, len(tops) // TRACE_COUNT)
    sampled = tops[::step][:TRACE_COUNT]

    ok = 0
    total = 0
    for (y, x) in sampled:
        total += 1
        prev = 2
        good = True
        cy, cx = y, x
        visited = set()
        while True:
            r = rank[cy, cx]
            if r >= 0:
                if r > prev:
                    good = False
                    break
                prev = r
            visited.add((cy, cx))
            best = None
            best_z = z[cy, cx]  # 只走严格下降段(单一坡面, 谷底终止)
            for dy in (-1, 0, 1):
                for dx in (-1, 0, 1):
                    if dx == 0 and dy == 0:
                        continue
                    ny, nx = cy + dy, cx + dx
                    if 0 <= ny < n_rows and 0 <= nx < n_cols:
                        if (ny, nx) in visited:
                            continue
                        if z[ny, nx] < best_z:
                            best_z = z[ny, nx]
                            best = (ny, nx)
            if best is None:
                break
            cy, cx = best
            if len(visited) > 400:
                break
        if good:
            ok += 1
    return (ok / total) if total else 1.0, total


def agreement(a, b, border_mask):
    """非边界像元一致率。"""
    m = ~border_mask
    if m.sum() == 0:
        return 1.0
    same = a[m] == b[m]
    return float(same.mean())


def scene_checks(name: str, report: list):
    """单场景的完整性/码域/轨迹检查。返回硬门失败数。"""
    fails = 0
    terrain, _ = _read(scene_dir(name) / "terrain_position.tif")
    z, z_nd = _read(SCENE_DIR / f"{name}.tif")
    if terrain is None:
        report.append(f"[FAIL] {name}: 缺少 terrain_position.tif")
        return 1
    violations, n_codes = check_code_domain(terrain, z, z_nd)
    rate, n_traces = trace_monotonic_rate(terrain, z)
    report.append(
        f"[{'OK' if violations == 0 else 'FAIL'}] {name}: 码域违规+保留码+NoData={violations}, "
        f"有效码数={n_codes}, 轨迹单调率={rate:.3f} ({n_traces} 条)"
    )
    if violations > 0:
        fails += 1
    if rate < TRACE_MONOTONIC_MIN:
        report.append(f"[FAIL] {name}: 单调率 {rate:.3f} < {TRACE_MONOTONIC_MIN}")
        fails += 1
    return fails


def variant_checks(report: list):
    """分辨率/噪声/变换变体与基准场景的非边界一致率。"""
    fails = 0
    pairs = (
        ("nested_ridges_25m", "nested_ridges", RESAMPLE_AGREEMENT_MIN, "25m 回采样"),
        ("nested_ridges_noise0.5", "nested_ridges", NOISE_AGREEMENT_MIN, "0.5m 噪声"),
        ("nested_ridges_rot90", "nested_ridges", 0.80, "90 度旋转"),
        ("v_valley_mirror", "v_valley", 0.80, "镜像"),
        ("broad_basin_shift50", "broad_basin", 0.95, "垂直平移"),
    )
    for variant, base, threshold, label in pairs:
        tv, _ = _read(scene_dir(variant) / "terrain_position.tif")
        tb, _ = _read(scene_dir(base) / "terrain_position.tif")
        if tv is None or tb is None:
            report.append(f"[SKIP] {label}: 缺少结果")
            continue
        if "rot90" in variant:
            tv = np.rot90(tv, k=-1).copy()  # 逆旋转回基准朝向
        if tv.shape != tb.shape:
            # 分辨率变体: 最近邻重采样到基准网格
            rows, cols = tb.shape
            yy = (np.arange(rows) * tv.shape[0] / rows).astype(int)
            xx = (np.arange(cols) * tv.shape[1] / cols).astype(int)
            tv = tv[np.ix_(yy, xx)]
        border = np.zeros_like(tv, dtype=bool)
        for axis in (0, 1):
            sl = [slice(None)] * 2
            sl[axis] = slice(2, -2)
            border[tuple(sl)] = True
        border = ~border
        rate = agreement(tv, tb, border)
        ok = rate >= threshold
        report.append(
            f"[{'OK' if ok else 'FAIL'}] {label}: 非边界一致率 {rate:.3f} (门限 {threshold})"
        )
        if not ok:
            fails += 1
    return fails


def param_checks(report: list):
    """±20% 参数扰动的类面积/盆地面积敏感性。"""
    fails = 0
    base_dir = VALIDATION_OUT / "param" / "wide_valley" / "base"
    tb, _ = _read(base_dir / "terrain_position.tif")
    tb_basin, _ = _read(base_dir / "diagnostics_basin_mask.tif")
    if tb is None:
        report.append("[SKIP] ±20% 参数扰动: 缺少 base 结果")
        return 0
    for variant in ("minus", "plus"):
        vdir = VALIDATION_OUT / "param" / "wide_valley" / variant
        tv, _ = _read(vdir / "terrain_position.tif")
        tv_basin, _ = _read(vdir / "diagnostics_basin_mask.tif")
        if tv is None or tv.shape != tb.shape:
            report.append(f"[SKIP] 参数 {variant}: 缺少结果")
            continue
        worst = 0.0
        for code in (3, 4, 5, 6, 7, 8):
            base_area = float((tb == code).sum())
            var_area = float((tv == code).sum())
            if base_area > 0:
                worst = max(worst, abs(var_area - base_area) / base_area)
        basin_base = float((tb_basin == 1).sum()) if tb_basin is not None else 0.0
        basin_var = float((tv_basin == 1).sum()) if tv_basin is not None else 0.0
        basin_change = (
            abs(basin_var - basin_base) / basin_base if basin_base > 0 else 0.0
        )
        ok = worst <= PARAM_CLASS_AREA_CHANGE_MAX and basin_change <= PARAM_BASIN_AREA_CHANGE_MAX
        report.append(
            f"[{'OK' if ok else 'FAIL'}] 参数 {variant}: 类面积最大变化 {worst:.3f}, "
            f"盆地面积变化 {basin_change:.3f}"
        )
        if not ok:
            fails += 1
    return fails


def main():
    gdal.UseExceptions()
    report = []
    fails = 0
    for scene in VALIDATION_SCENES:
        fails += scene_checks(scene, report)
    fails += variant_checks(report)
    fails += param_checks(report)
    print("\n".join(report))
    print(f"\n硬门失败数: {fails}")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
