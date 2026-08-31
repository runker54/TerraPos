# TerraPos 桌面应用（Rust）

地形部位划分工具的产品实现：eframe/egui 桌面应用 + topo_core 纯 Rust 算法库。
无 GDAL / ArcGIS / Python 依赖，静态链接 CRT，单文件免安装。

## 构建

```bash
cargo build --release
# 产物: target/release/topo_app.exe
```

要求 Rust 1.75+（MSVC 工具链，已配置静态 CRT，无 VC++ Redistributable 依赖）。

## 命令行批处理（无 UI）

```bash
cargo run --release -p topo_core --example run_full -- <dem.tif> <out_dir>
```

## 界面

- **划分工作台**：左侧参数分组卡片（全部判断指标可调，含恢复默认），
  中央结果预览（地形部位/地貌亚类图层切换、鼠标锚点缩放、拖拽平移、浮动图例），
  顶部统计卡带（各类面积与占比条形）；
- **参数详解**：全部判断指标的地貌含义、默认值与调整方向说明；
- 运行在后台线程，分阶段进度条 + 日志 + 取消。

## 算法模块（topo_core）

| 模块 | 内容 |
|---|---|
| `geotiff` | 轻量 GeoTIFF 读写（IFD 解析 + GeoTags 注入，无 GDAL） |
| `filter` | 分离均值/极值滤波（rayon 并行，nearest 外推与 scipy 一致） |
| `terrain` | 坡度(numpy.gradient 对齐)/TPI/起伏度/box 下采样 |
| `hydro` | Priority-Flood 填洼(z-limit) + D8 流向 + 汇流累积 |
| `distance` | 精确欧氏 EDT（Felzenszwalb 两遍 + 最近源索引） |
| `segment` | 峰顶 NMS + marker watershed + 连通域 |
| `pipeline` | 参数模型 + 编排 + 进度/取消 |

验收：17 项单元测试、与 Python 基准对拍一致率 97.17%、3358km²@5m 全流程约 22 秒。
