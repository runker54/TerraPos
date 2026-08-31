# TerraPos 地形部位划分工具

[简体中文](README.md) | [English](README_EN.md)

从 DEM 自动划分西南区 8 类地形部位：山间盆地、丘陵上/中/下、山地坡上/中/下。

## 地形部位编码

| 编码 | 部位 | 判据摘要 |
|---|---|---|
| 1 | 山间盆地 | 坡度<6° + TPI<-25m(2km) + 面积≥0.5km² + 宽度>250m |
| 2 | 宽谷盆地 | 编码保留 |
| 3/4/5 | 丘陵上/中/下 | 海拔<500m 且山体单元内位置三分 |
| 6/7/8 | 山地坡上/中/下 | 海拔≥500m 且山体单元内位置三分 |

地貌亚类栅格：低丘 / 高丘 / 低山 / 中山 / 高山 / 极高山 / 平坝。

## 使用

```bash
cd rust
cargo build --release
./target/release/topo_app.exe
```

命令行批处理：

```bash
cargo run --release -p topo_core --example run_full -- <dem.tif> <out_dir>
```

界面：导入 DEM → 调整参数（可默认）→ 运行 → 预览导出。
参数含义见应用内「参数详解」页。

## 算法要点

1. Priority-Flood 填洼（z-limit 保留喀斯特深洼）→ D8 汇流累积；
2. 精确欧氏 EDT → HAND 阶地/坡麓修正；
3. 峰顶 watershed 山体单元，单元内位置三等分坡位；
4. 坝子逐像元宽度核心化（窄谷带归入坡下）；
5. 丘陵/山地按海拔 500m 分界。

## 发布

发布包（TerraPos-v0.0.1-win64.zip）见 GitHub [Releases](https://github.com/runker54/TerraPos/releases)，解压即用；
发布说明见 docs/RELEASE_NOTES-v0.0.1.md。
