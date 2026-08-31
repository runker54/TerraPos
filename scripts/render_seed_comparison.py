# -*- coding: utf-8 -*-
"""种子/个体分割策略四方案视觉对比渲染(固定路径, 无命令行参数)。

输入: rust/target/compare/{challenge_dem.tif, truth.txt, fixed|zoned|prom|scale/...}
输出: docs/seed_mode_comparison.png
"""
import os
import numpy as np
from osgeo import gdal
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei"]
plt.rcParams["axes.unicode_minus"] = False

BASE = os.path.join(os.path.dirname(__file__), "..", "rust", "target", "compare")
OUT = os.path.join(os.path.dirname(__file__), "..", "docs", "seed_mode_comparison.png")
# 丘陵链放大窗口 (x0, x1, y0, y1)
ZOOM = (420, 540, 380, 450)


def read_tif(path):
    ds = gdal.Open(path)
    band = ds.GetRasterBand(1)
    a = band.ReadAsArray()
    ds = None
    return a


def hillshade(dem, res=25.0, az=315.0, alt=45.0):
    gy, gx = np.gradient(dem, res)
    slope = np.arctan(np.hypot(gx, gy))
    aspect = np.arctan2(-gx, -gy)
    az, alt = np.radians(az), np.radians(alt)
    hs = np.sin(alt) * np.cos(slope) + np.cos(alt) * np.sin(slope) * np.cos(az - aspect)
    return np.clip(hs, 0, 1)


def unit_bounds(units):
    b = np.zeros_like(units, dtype=bool)
    b[:-1, :] |= units[:-1, :] != units[1:, :]
    b[1:, :] |= units[:-1, :] != units[1:, :]
    b[:, :-1] |= units[:, :-1] != units[:, 1:]
    b[:, 1:] |= units[:, :-1] != units[:, 1:]
    return b & (units > 0)


TERRAIN_COLORS = {
    1: ((51, 178, 229), "山间盆地"),
    3: ((250, 217, 89), "丘陵上"),
    4: ((217, 237, 166), "丘陵中"),
    5: ((153, 199, 102), "丘陵下"),
    6: ((250, 165, 60), "山地坡上"),
    7: ((222, 100, 50), "山地坡中"),
    8: ((107, 68, 35), "山地坡下"),
}


def terrain_rgba(tp):
    out = np.ones(tp.shape + (3,), dtype=np.float32)
    for code, (rgb, _) in TERRAIN_COLORS.items():
        out[tp == code] = np.array(rgb, dtype=np.float32) / 255.0
    return out


def main():
    dem = read_tif(os.path.join(BASE, "challenge_dem.tif"))
    hs = hillshade(dem)
    truth = []
    with open(os.path.join(BASE, "truth.txt"), encoding="utf-8") as f:
        for line in f:
            x, y, name = line.split()
            truth.append((int(x), int(y), name))

    modes = [
        ("fixed", "A. 单一固定窗口 500m (现状默认)"),
        ("zoned", "B. 分亚类动态窗口 250~900m"),
        ("prom", "C. 地形突起度 prominence ≥50m"),
        ("scale", "D. 多尺度 TPI 特征尺度投票"),
        ("hybrid", "E. 混合: 突出度∪距离语义 (推荐)"),
    ]

    fig, axes = plt.subplots(5, 3, figsize=(19, 31),
                             gridspec_kw={"width_ratios": [1.6, 1.6, 1.0]})
    x0, x1, y0, y1 = ZOOM
    n_units_all = {}
    for r, (mdir, title) in enumerate(modes):
        d = os.path.join(BASE, mdir)
        units = read_tif(os.path.join(d, "units_coarse.tif")).astype(int)
        tp = read_tif(os.path.join(d, "terrain_position.tif")).astype(int)
        n_units_all[mdir] = int(units.max())
        eb = unit_bounds(units)

        # 列1: 全域山体单元分割
        ax = axes[r][0]
        ax.imshow(hs, cmap="gray", vmin=0.15, vmax=0.95)
        ax.imshow(np.ma.masked_where(~eb, eb), cmap="autumn", vmin=0, vmax=1)
        for (x, y, name) in truth:
            ax.plot(x, y, "r+", ms=11, mew=2)
            ax.annotate(name, (x, y), color="yellow", fontsize=8,
                        textcoords="offset points", xytext=(4, 3))
        ax.add_patch(plt.Rectangle((x0, y0), x1 - x0, y1 - y0,
                                   fill=False, ec="lime", lw=1.5))
        ax.set_title(f"{title} — 山体单元 {units.max()} 个", fontsize=11, loc="left")
        ax.set_xlim(0, 600)
        ax.set_ylim(600, 0)

        # 列2: 8 类地形部位
        ax = axes[r][1]
        ax.imshow(hs, cmap="gray", vmin=0.15, vmax=0.55)
        ax.imshow(terrain_rgba(tp), alpha=0.82)
        ax.set_title("8 类地形部位划分", fontsize=11, loc="left")
        ax.set_xlim(0, 600)
        ax.set_ylim(600, 0)

        # 列3: 丘陵链放大(单元边界)
        ax = axes[r][2]
        ax.imshow(hs[y0:y1, x0:x1], cmap="gray", vmin=0.15, vmax=0.95)
        ax.imshow(np.ma.masked_where(~eb[y0:y1, x0:x1], eb[y0:y1, x0:x1]),
                  cmap="cool", vmin=0, vmax=1)
        for (x, y, name) in truth:
            if x0 <= x < x1 and y0 <= y < y1:
                ax.plot(x - x0, y - y0, "r+", ms=10, mew=2)
        ax.set_title("丘陵链放大", fontsize=11, loc="left")
        ax.set_xlim(0, x1 - x0)
        ax.set_ylim(y1 - y0, 0)

    handles = [Patch(facecolor=np.array(rgb) / 255, label=nm)
               for (_, (rgb, nm)) in TERRAIN_COLORS.items()]
    handles.append(Line2D([0], [0], marker="+", color="r", lw=0, ms=10, label="合成真值峰位"))
    fig.legend(handles=handles, loc="lower center", ncol=8, fontsize=10, frameon=False)
    fig.suptitle("山体个体分割策略四方案对比(合成挑战地形 15km×15km @25m)", fontsize=15)
    fig.tight_layout(rect=(0, 0.025, 1, 0.975))
    fig.savefig(OUT, dpi=105)
    print("已输出:", os.path.abspath(OUT))
    print("单元数:", n_units_all)


if __name__ == "__main__":
    main()
