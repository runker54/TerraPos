# -*- coding: utf-8 -*-
"""
矢量化输出：地形部位栅格 -> shapefile
同时更新输出tif的色表(与最新配色一致)
"""
import os
import numpy as np
import rasterio
from rasterio.features import shapes
from shapely.geometry import shape
import geopandas as gpd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, 'output')
TIF = os.path.join(OUT, 'terrain_position_5m.tif')
SHP = os.path.join(OUT, 'terrain_position_5m.shp')
SIMPLIFY_TOL = 15.0  # 拓扑简化容差(米, 3像元)

CLASS_NAMES = {
    1: '山间盆地', 2: '宽谷盆地',
    3: '丘陵上部', 4: '丘陵中部', 5: '丘陵下部',
    6: '山地坡上', 7: '山地坡中', 8: '山地坡下',
}
CLASS_COLORS = {
    1: (51, 178, 229), 2: (102, 217, 242),
    3: (250, 217, 89), 4: (217, 237, 166), 5: (153, 199, 102),
    6: (250, 165, 60), 7: (222, 100, 50), 8: (140, 60, 40),
}

print('读取栅格...', flush=True)
with rasterio.open(TIF) as src:
    arr = src.read(1)
    transform = src.transform
    crs = src.crs

# 更新色表
with rasterio.open(TIF, 'r+') as dst:
    dst.write_colormap(1, {c: CLASS_COLORS[c] + (255,) for c in CLASS_COLORS} | {0: (0, 0, 0, 0)})

print('栅格转面...', flush=True)
recs = []
for geom, val in shapes(arr, mask=arr > 0, transform=transform, connectivity=8):
    recs.append((int(val), shape(geom)))
print(f'图斑总数: {len(recs)}', flush=True)

print(f'拓扑简化(容差{SIMPLIFY_TOL}m)...', flush=True)
geoms = [g.simplify(SIMPLIFY_TOL, preserve_topology=True) for _, g in recs]
vals = [v for v, _ in recs]

gdf = gpd.GeoDataFrame({
    'DLBW': vals,
    'MC': [CLASS_NAMES[v] for v in vals],
}, geometry=geoms, crs=crs)

print('写出shapefile...', flush=True)
gdf.to_file(SHP, driver='ESRI Shapefile', encoding='utf-8')

# 图斑面积统计
gdf['AREA_KM2'] = gdf.geometry.area / 1e6
stat = gdf.groupby(['DLBW', 'MC'])['AREA_KM2'].agg(['count', 'sum'])
print(stat.to_string(), flush=True)
print('完成:', SHP, flush=True)
