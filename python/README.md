# Python 验证工具(科学验证层)

本目录脚本是与 Rust 生产实现分离的**验证与诊断工具**。
全部脚本把常量放在文件顶部, **不使用命令行参数**。
统一使用 `D:/worker_code/.venvgis/Scripts/python.exe` 运行(需 GDAL/NumPy/SciPy/matplotlib)。

## 验证脚本(自适应地形部位划分管线)

| 顺序 | 脚本 | 输入 | 输出 | 用途 |
|---|---|---|---|---|
| 1 | `build_synthetic_dem.py` | (无) | `data/synthetic/generated/*.tif` | 生成 13 个解析合成场景(宽盆/V 谷/嵌套脊/锥丘/闭合洼地/参考小样区)与平移、镜像、旋转、噪声、分辨率变体; 每场景附期望语义区掩膜(`*_expected.tif`, 1 盆地 2 上部 3 中部 4 下部 0 无解析预期) |
| 2 | `rust/topo_core/examples/run_validation.rs` | 场景 tif | `data/validation_out/synthetic/<场景>/` 与 `param/` | Rust 管线批量运行(含 ±20% 参数扰动组), 供验收判定 |
| 3 | `adaptive_terrain_reference.py` | `ref_small` 小样区 | `data/synthetic/reference/*.npz` | NumPy/SciPy 参考公式(稳健起伏/尺度选择/D8-HAND/qd-qz-q/模糊隶属度)——独立公式 oracle, 不是第二套生产分类器 |
| 4 | `validate_terrain_result.py` | `data/validation_out/` | 终端指标 + 退出码 | 验收硬门: 码域({1,3-8})/保留码 2 计数 0/NoData 保留/25m 回采样一致率≥85%/0.5m 噪声≥90%/90 度旋转≥0.8/镜像/垂直平移/±20% 参数类面积≤15% 盆地≤20%/脊-谷轨迹单调率≥95%。任一失败退出码非零 |
| 5 | `render_qa_atlas.py` | `data/validation_out/synthetic/broad_basin` | `data/qa/broad_basin_atlas.png` | 4×4 固定版式图集(高程/修正深度/河网/脊线/单元/尺度/HAND/q/坡位/盆候选/核/重建/最终类别/置信度/直方图/报告) |

### 用户可编辑常量

- `build_synthetic_dem.py`: `EPSG`(4545)、`NOISE_SEED`(42)、`NOISE_SIGMAS`、`RESOLUTIONS_M`
- `validate_terrain_result.py`: 全部 `*_MIN/*_MAX` 阈值(来自设计规格 14/15 节, **不得为通过验收而放宽**)、`VALIDATION_SCENES`
- `adaptive_terrain_reference.py`: `SCALES_M`、`GROWTH_THRESHOLD`、`MEMBERS`(隶属函数参数)
- `render_qa_atlas.py`: `SCENE`(图集场景)、`CLASS_CMAP`(类别配色)

### 结果语义

- 分类编码与生产完全一致: `0=NoData, 1=山间/宽谷盆地, 2=保留码(永不生成), 3/4/5=丘陵上/中/下, 6/7/8=山地坡上/中/下`
- 合成场景无人工真值, 验收依赖**解析预期+拓扑一致性+尺度稳定性**三类证据;
  通过验收不等于传统意义的"分类精度"达标。
- 参考实现的朴素滤波内存受限, 仅在小样区(`ref_small`)上运行。

---

# 旧版 Python 基准(已被自适应管线替代, 仅留档对照)

`terrain_position_main.py` 及 arcpy 脚本为旧固定窗口方案(对拍一致率 97.17%),
新算法验证通过后已不再用于生产; 依赖 ArcGIS Pro(arcpy, Advanced)、
rasterio/scipy/numpy, 运行顺序与说明见该目录下各脚本头部配置区。
