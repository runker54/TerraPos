# -*- coding: utf-8 -*-
"""高倍放大验证：大山体小地形内部的上中下区分（DEM阴影 与 坡位图 同区对照）"""
import os
import numpy as np
import rasterio
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

for f in ['Microsoft YaHei', 'SimHei']:
    if any(f.lower() == x.name.lower() for x in font_manager.fontManager.ttflist):
        plt.rcParams['font.sans-serif'] = [f]
        break
plt.rcParams['axes.unicode_minus'] = False

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORK = os.path.join(ROOT, 'work')
OUT = os.path.join(ROOT, 'output')

NAMES = {1: '山间盆地', 6: '山地坡上', 7: '山地坡中', 8: '山地坡下'}
COLORS = {1: '#33B2E5', 6: '#FAA53C', 7: '#DE6432', 8: '#6B4423'}

with rasterio.open(os.path.join(WORK, 'dem25.tif')) as src:
    dem25 = src.read(1).astype(np.float32)
result = np.load(os.path.join(WORK, 'result_5m.npy'))[::5, ::5][:dem25.shape[0], :dem25.shape[1]]

def hillshade(a, res, azi=315, alt=45):
    gy, gx = np.gradient(a.astype(np.float64), res)
    az = np.radians(450 - azi); al = np.radians(alt)
    s = np.arctan(np.hypot(gx, gy)); asp = np.arctan2(-gx, gy)
    return np.clip(np.sin(al)*np.cos(s) + np.cos(al)*np.sin(s)*np.cos(az - asp), 0, 1)

hs = hillshade(dem25, 25)

# 高倍放大区：中部山体区 600x600 像元(15km栅格/25m) -> 每幅显示约15km,细看丘包
regions = [
    ('A区 (y1500-1900, x1100-1700)', np.s_[1500:1900, 1100:1700]),
    ('B区 (y300-700, x1500-2100)', np.s_[300:700, 1500:2100]),
]
fig, axes = plt.subplots(2, 2, figsize=(19, 13), dpi=105)
for row, (title, sl) in enumerate(regions):
    ax0 = axes[row, 0]
    ax0.imshow(hs[sl], cmap='gray')
    ax0.set_title(f'DEM山体阴影 {title}', fontsize=12)
    ax0.axis('off')
    ax1 = axes[row, 1]
    ax1.imshow(hs[sl], cmap='gray', alpha=0.35)
    rgb = np.zeros(result[sl].shape + (4,))
    for v, c in COLORS.items():
        rgb[result[sl] == v] = matplotlib.colors.to_rgba(c, 0.85)
    ax1.imshow(rgb)
    ax1.set_title(f'地形部位 {title}', fontsize=12)
    ax1.axis('off')
handles = [plt.Rectangle((0, 0), 1, 1, fc=c) for c in COLORS.values()]
fig.legend(handles, [NAMES[v] for v in COLORS], loc='lower center', ncol=4, fontsize=12)
plt.tight_layout(rect=[0, 0.035, 1, 1])
p = os.path.join(WORK, 'check_zoom_small.png')
plt.savefig(p, bbox_inches='tight')
plt.close()
print('已保存:', p)
