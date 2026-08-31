# -*- coding: utf-8 -*-
"""
结果验证 v3：空间检查图(全区+3放大) + 地形剖面验证图(高程-坡位对应关系)
"""
import os
import numpy as np
import rasterio
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import ListedColormap, BoundaryNorm

for f in ['Microsoft YaHei', 'SimHei']:
    if any(f.lower() == x.name.lower() for x in font_manager.fontManager.ttflist):
        plt.rcParams['font.sans-serif'] = [f]
        break
plt.rcParams['axes.unicode_minus'] = False

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORK = os.path.join(ROOT, 'work')
OUT = os.path.join(ROOT, 'output')
STEP = 5  # 5m->25m 渲染步长

NAMES = {1: '山间盆地', 3: '丘陵上部', 4: '丘陵中部', 5: '丘陵下部',
         6: '山地坡上', 7: '山地坡中', 8: '山地坡下'}
COLORS = {1: '#33B2E5', 3: '#F9D959', 4: '#D9EDA6', 5: '#99C766',
          6: '#FAA53C', 7: '#DE6432', 8: '#6B4423'}
SUB_NAMES = {0: '无', 1: '低丘', 2: '高丘', 3: '低山', 4: '中山', 5: '高山', 6: '极高山', 7: '平坝'}
SUB_COLORS = {0: '#000000', 1: '#B4E696', 2: '#6EC36E', 3: '#FAE182',
              4: '#EBAA5A', 5: '#CD6E4B', 6: '#96413C', 7: '#55B9EB'}

with rasterio.open(os.path.join(WORK, 'dem25.tif')) as src:
    dem25 = src.read(1).astype(np.float32)
    _nod = src.nodata
if _nod is not None:
    dem25[dem25 == _nod] = np.nan
_bad = ~np.isfinite(dem25)
if _bad.any():
    from scipy import ndimage as _ndi
    _, (_iy, _ix) = _ndi.distance_transform_edt(_bad, return_indices=True)
    dem25[_bad] = dem25[_iy[_bad], _ix[_bad]]
result = np.load(os.path.join(WORK, 'result_5m.npy'))[::STEP, ::STEP]
with rasterio.open(os.path.join(OUT, 'geomorph_subclass_5m.tif')) as src:
    sub = src.read(1)[::STEP, ::STEP]
# 与25m底图对齐(5m降采样可能多1行/列)
result = result[:dem25.shape[0], :dem25.shape[1]]
sub = sub[:dem25.shape[0], :dem25.shape[1]]

def hillshade(a, res, azi=315, alt=45):
    gy, gx = np.gradient(a.astype(np.float64), res)
    az = np.radians(450 - azi); al = np.radians(alt)
    s = np.arctan(np.hypot(gx, gy)); asp = np.arctan2(-gx, gy)
    return np.clip(np.sin(al)*np.cos(s) + np.cos(al)*np.sin(s)*np.cos(az - asp), 0, 1)

hs = hillshade(dem25, 25)

def class_rgba(cls_arr, colors):
    rgb = np.zeros(cls_arr.shape + (4,))
    for v, c in colors.items():
        rgb[cls_arr == v] = matplotlib.colors.to_rgba(c, 0.82)
    return rgb

# ---------- 图1: 空间检查 ----------
fig, axes = plt.subplots(2, 2, figsize=(20, 13), dpi=100)
views = [(np.s_[:, :], '全区地形部位概览'),
         (np.s_[600:1300, 900:1700], '放大① 河谷坝区'),
         (np.s_[100:700, 2300:2966], '放大② 山地坡位分带'),
         (np.s_[1300:1812, 1200:2000], '放大③ 丘陵/低山丘包')]
for ax, (sl, title) in zip(axes.flat, views):
    ax.imshow(hs[sl], cmap='gray')
    ax.imshow(class_rgba(result[sl], COLORS), interpolation='nearest')
    ax.set_title(title, fontsize=13); ax.axis('off')
handles = [plt.Rectangle((0, 0), 1, 1, fc=c) for c in COLORS.values()]
fig.legend(handles, [NAMES[v] for v in COLORS], loc='lower center', ncol=7, fontsize=11)
plt.tight_layout(rect=[0, 0.04, 1, 1])
p1 = os.path.join(WORK, 'check_full.png')
plt.savefig(p1, bbox_inches='tight'); plt.close()
print('已保存:', p1)

# ---------- 图2: 亚类检查 ----------
fig, axes = plt.subplots(1, 2, figsize=(20, 7), dpi=100)
for ax, (sl, title) in zip(axes, [(np.s_[:, :], '地貌亚类 全区'),
                                  (np.s_[100:900, 1200:2400], '地貌亚类 放大(西北部)')]):
    ax.imshow(hs[sl], cmap='gray')
    ax.imshow(class_rgba(sub[sl], SUB_COLORS), interpolation='nearest')
    ax.set_title(title, fontsize=13); ax.axis('off')
handles = [plt.Rectangle((0, 0), 1, 1, fc=c) for k, c in SUB_COLORS.items() if k > 0]
fig.legend(handles, [SUB_NAMES[k] for k in SUB_COLORS if k > 0], loc='lower center', ncol=7, fontsize=11)
plt.tight_layout(rect=[0, 0.05, 1, 1])
p2 = os.path.join(WORK, 'check_subclass.png')
plt.savefig(p2, bbox_inches='tight'); plt.close()
print('已保存:', p2)

# ---------- 图3: 地形剖面验证 ----------
# 剖面A: 东西向 y=820 穿中部河谷; 剖面B: 南北向 x=2400 穿东部山地
profiles = [('剖面A  东西向 (y=820)', dem25[820, :], result[820, :], sub[820, :]),
            ('剖面B  南北向 (x=2400)', dem25[:, 2400], result[:, 2400], sub[:, 2400])]
fig, axes = plt.subplots(len(profiles), 1, figsize=(18, 9), dpi=110, sharex=False)
for ax, (title, z, cls, subp) in zip(axes, profiles):
    x = np.arange(len(z)) * 25 / 1000.0  # km
    ymin = z.min() - 60
    # 分段着色填充
    for v in COLORS:
        m = cls == v
        if m.any():
            ax.fill_between(x[m], ymin, z[m], color=COLORS[v], alpha=0.75, linewidth=0)
    ax.plot(x, z, 'k-', lw=0.9)
    # 亚类分段条带
    for k in SUB_COLORS:
        if k > 0:
            m = subp == k
            if m.any():
                ax.fill_between(x[m], ymin - 90, ymin - 15, color=SUB_COLORS[k], linewidth=0)
    ax.set_ylim(ymin - 95, z.max() + 60)
    ax.set_title(title + '    (上部:高程-坡位着色 / 下部色条:地貌亚类)', fontsize=12)
    ax.set_ylabel('高程 (m)')
    ax.set_xlim(x[0], x[-1])
axes[-1].set_xlabel('距离 (km)')
legend_items = [plt.Rectangle((0, 0), 1, 1, fc=c) for c in COLORS.values()]
legend_labels = [NAMES[v] for v in COLORS]
for k in [3, 4, 7]:
    legend_items.append(plt.Rectangle((0, 0), 1, 1, fc=SUB_COLORS[k]))
    legend_labels.append('亚类:' + SUB_NAMES[k])
axes[0].legend(legend_items, legend_labels, loc='upper right', ncol=4, fontsize=9)
plt.tight_layout()
p3 = os.path.join(WORK, 'check_profile.png')
plt.savefig(p3, bbox_inches='tight'); plt.close()
print('已保存:', p3)

# 同步成品tif色表
import rasterio as _rio
def hex2rgb(h):
    h = h.lstrip('#')
    return tuple(int(h[i:i+2], 16) for i in (0, 2, 4))
with _rio.open(os.path.join(OUT, 'terrain_position_5m.tif'), 'r+') as dst:
    dst.write_colormap(1, {v: hex2rgb(c) + (255,) for v, c in COLORS.items()} | {0: (0, 0, 0, 0)})
print('已更新 tif 色表')
