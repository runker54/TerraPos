# -*- coding: utf-8 -*-
"""检查 Manba result.shp 各面位置与我的坝子覆盖的关系"""
import numpy as np
from osgeo import gdal

DEM = r"G:\tif_features\county_feature\hhgq\dem.tif"
MINE = r"E:\zcode_worker\Topographic\rust\target\hhgq_out\basin_mask.tif"

# 用 gdal 读 shp 的 envvelope(简化, 不用 geopandas 避免 segfault)
ds_vec = gdal.OpenEx(r"G:\soil_shp_data\hhgq\river\alluvial_areas_final_result.shp", gdal.OF_VECTOR)
lyr = ds_vec.GetLayer()
print(f"面要素 {lyr.GetFeatureCount()} 个")
gt_dem = gdal.Open(DEM).GetGeoTransform()
W = gdal.Open(DEM).RasterXSize
H = gdal.Open(DEM).RasterYSize
mine_ds = gdal.Open(MINE)
mine_b = mine_ds.GetRasterBand(1)
mine_gt = mine_ds.GetGeoTransform()

for feat in lyr:
    geom = feat.GetGeometryRef()
    env = geom.GetEnvelope()
    cx = (env[0] + env[1]) / 2
    cy = (env[2] + env[3]) / 2
    area_km2 = geom.GetArea() / 1e6
    col = int((cx - gt_dem[0]) / gt_dem[1])
    row = int((gt_dem[3] - cy) / -gt_dem[5])
    in_dem = 0 <= col < W and 0 <= row < H
    # 读该位置的坝子值
    has = "N/A"
    if in_dem:
        val = mine_b.ReadAsArray(row, col, 1, 1)[0][0]
        has = f"坝子={val}"
    print(f"  面积={area_km2:.2f}km² 中心=({col},{row}) 在DEM={in_dem} {has}")
ds_vec = None
