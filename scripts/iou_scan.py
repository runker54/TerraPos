# -*- coding: utf-8 -*-
"""融合版(多河网阈值) vs Manba 参数版原始坝区 IoU 扫描(固定路径)"""
import numpy as np
import geopandas as gpd
import rasterio
from rasterio.features import rasterize
from rasterio.transform import Affine

DEM = r"G:\tif_features\county_feature\hhgq\dem.tif"
HIS_SHP = r"G:\soil_shp_data\hhgq\river\alluvial_areas_final_500_5_5_5_5000_50.shp"
FUSED = r"E:\zcode_worker\Topographic\rust\target\basin_fused\fused_{a}.tif"

ds = rasterio.open(DEM)
W, H = ds.width, ds.height
g0 = ds.transform
step = 4
w, h = W // step, H // step
transform = Affine(g0.a * step, g0.b, g0.c, g0.d, g0.e * step, g0.f)
gdf = gpd.read_file(HIS_SHP)
his = rasterize([(g, 1) for g in gdf.geometry], out_shape=(h, w), transform=transform, fill=0, dtype="uint8") == 1
del gdf
cell = 20.0 * 20.0 / 1e6
a_his = his.sum() * cell
print(f"Manba 参数版: {a_his:.1f} km2")

for a in ["0.3", "1", "2", "5", "10", "20"]:
    mine = rasterio.open(FUSED.format(a=a)).read(1, out_shape=(h, w)) == 1
    inter = (his & mine).sum() * cell
    union = (his | mine).sum() * cell
    miss = 100 * (his & ~mine).sum() / max(1, his.sum())
    extra = (mine & ~his).sum() * cell
    print(f"acc={a:>4}km2: mine {mine.sum()*cell:5.1f} | IoU {100*inter/union:5.1f}% | 漏检 {miss:4.0f}% | 多检 {extra:5.1f}km2")
