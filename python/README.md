# Python 基准实现

地形部位划分的 Python 参考实现（ArcGIS Pro arcpy + scipy/rasterio）。
用于验证 Rust 应用的结果一致性（对拍一致率 97.17%），以及需要 arcpy 高级
水文分析能力的场景。**生产使用请用 `rust/` 桌面应用。**

## 依赖

| 环境 | 用途 |
|---|---|
| ArcGIS Pro (arcpy, Advanced) | 填洼 / 流向 / 汇流累积 / 负地形脊线 |
| Python: rasterio, scipy, numpy, matplotlib, PIL | 其余因子与可视化 |

## 运行顺序

```bash
# 0) 生成 25m 中间层
gdalwarp -tr 25 25 -r average data/dem.tif work/dem25.tif

# 1) 水文分析(arcpy 环境): 正地形河网
python extract_river_arcpy.py    # -> work/acc25_arcpy.tif
#    (extract_ridge_arcpy.py 为可选的负地形脊线提取, 当前流程未使用)

# 2) 主流程(GIS python 环境): 输出 output/terrain_position_5m.tif 等
python terrain_position_main.py

# 3) 矢量化与验证图
python raster_to_shp.py
python check_result.py           # 空间/亚类/剖面三张检查图
```

脚本内路径均为写死的绝对路径（data/dem.tif、work/、output/），
换区域时同步修改脚本头部配置区。

## 脚本说明

| 脚本 | 说明 |
|---|---|
| `terrain_position_main.py` | 主流程：地貌亚类 + HAND/RHAND 坡位 + 坝子 + 后处理 |
| `extract_river_arcpy.py` | 正地形填洼(15m z-limit)→D8→汇流累积 |
| `extract_ridge_arcpy.py` | 负地形汇流累积(山脊线对偶) |
| `raster_to_shp.py` | 分类栅格 → Shapefile |
| `check_result.py` | 空间检查图 / 亚类图 / 地形剖面图 |
| `check_zoom.py` | 小地形高倍 DEM 阴影与坡位对照 |

## 与 Rust 版的差异

- 本实现的河网/脊线由 arcpy 生成；Rust 版已内置 Priority-Flood 水文，无需 arcpy；
- 两版在 watershed 平局顺序、重采样实现上存在固有差异，对拍总体一致率 97.17%。
