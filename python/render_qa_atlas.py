# -*- coding: utf-8 -*-
"""固定版式 QA 图集: 4x4 面板对照 DEM 与全部诊断层。

读取 data/validation_out 下指定场景的正式成果与诊断层,
输出 data/qa/<场景>_atlas.png。常量在文件顶部; 无命令行参数。
运行: D:/worker_code/.venvgis/Scripts/python.exe python/render_qa_atlas.py
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from osgeo import gdal

ROOT = Path(__file__).resolve().parents[1]
VALIDATION_OUT = ROOT / "data" / "validation_out" / "synthetic"
QA_DIR = ROOT / "data" / "qa"
SCENE = "broad_basin"  # 图集目标场景
CLASS_CMAP = {
    0: (1.0, 1.0, 1.0),
    1: (0.47, 0.75, 0.47),
    3: (0.84, 0.62, 0.62),
    4: (0.90, 0.78, 0.55),
    5: (0.67, 0.86, 0.63),
    6: (0.66, 0.44, 0.28),
    7: (0.87, 0.77, 0.47),
    8: (0.52, 0.66, 0.38),
}
PANELS = (
    ("DEM 高程", "dem.tif", "terrain"),
    ("水文修正深度", "diagnostics/hydro_conditioning_depth.tif", "continuous"),
    ("河网等级", "diagnostics/stream_level.tif", "stream"),
    ("脊线", "diagnostics/ridge_mask.tif", "mask"),
    ("坡面单元", "diagnostics/slope_unit.tif", "unit"),
    ("自适应尺度 R*", "diagnostics/adaptive_scale_m.tif", "continuous"),
    ("HAND", "diagnostics/hand_m.tif", "continuous"),
    ("相对位置 q", "diagnostics/relative_position.tif", "continuous"),
    ("原始坡位", "diagnostics/slope_position_raw.tif", "class"),
    ("盆地候选", "diagnostics/basin_candidate.tif", "mask"),
    ("盆地核心", "diagnostics/basin_core.tif", "mask"),
    ("盆地重建", "diagnostics/basin_mask.tif", "mask"),
    ("最终坡位", "terrain_position.tif", "class"),
    ("置信度", "terrain_confidence.tif", "continuous"),
    ("类别直方图", "@histogram", "hist"),
    ("报告摘要", "@report", "text"),
)


def _read(scene_dir: Path, rel: str):
    path = scene_dir / rel
    if not path.exists():
        return None
    ds = gdal.Open(str(path))
    arr = ds.GetRasterBand(1).ReadAsArray()
    ds = None
    return arr


def _show(ax, arr, kind, title):
    ax.set_title(title, fontsize=9)
    ax.set_xticks([])
    ax.set_yticks([])
    if arr is None:
        ax.text(0.5, 0.5, "缺失", ha="center", va="center", transform=ax.transAxes)
        return
    if kind == "terrain" or kind == "continuous":
        im = ax.imshow(arr, cmap="terrain" if kind == "terrain" else "viridis")
        plt.colorbar(im, ax=ax, fraction=0.04)
    elif kind == "class":
        cmap = matplotlib.colors.ListedColormap(
            [CLASS_CMAP.get(c, (0.8, 0.8, 0.8)) for c in range(9)]
        )
        ax.imshow(arr, cmap=cmap, vmin=0, vmax=8, interpolation="nearest")
    elif kind in ("mask", "stream", "unit"):
        ax.imshow(arr, cmap="viridis", interpolation="nearest")
    else:
        ax.axis("off")


def main():
    gdal.UseExceptions()
    scene_dir = VALIDATION_OUT / SCENE
    if not scene_dir.exists():
        print(f"缺少场景结果: {scene_dir}")
        sys_exit = 1
        raise SystemExit(sys_exit)
    QA_DIR.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(4, 4, figsize=(18, 14))
    for ax, (title, rel, kind) in zip(axes.ravel(), PANELS):
        if rel == "@histogram":
            terrain = _read(scene_dir, "terrain_position.tif")
            if terrain is not None:
                vals, counts = np.unique(terrain, return_counts=True)
                ax.bar([str(int(v)) for v in vals], counts)
            ax.set_title("类别直方图", fontsize=9)
            continue
        if rel == "@report":
            report_path = scene_dir / "class_report.txt"
            text = report_path.read_text(encoding="utf-8") if report_path.exists() else "缺失"
            ax.text(0.0, 1.0, text[:1200], va="top", fontsize=5, family="monospace")
            ax.set_title("class_report 摘要", fontsize=9)
            ax.set_xticks([])
            ax.set_yticks([])
            continue
        _show(ax, _read(scene_dir, rel), kind, title)
    fig.suptitle(f"QA 图集: {SCENE}", fontsize=13)
    fig.tight_layout()
    out = QA_DIR / f"{SCENE}_atlas.png"
    fig.savefig(out, dpi=110)
    print(f"QA 图集 -> {out}")


if __name__ == "__main__":
    main()
