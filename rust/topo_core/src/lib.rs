//! topo_core: 地形部位划分核心算法库（纯 Rust，无 GDAL 依赖）
//!
//! 模块:
//! - [`geotiff`]：轻量 GeoTIFF 读写
//! - [`input`]：输入检查、有效区管理与地貌表面构建
//! - [`filter`]：并行分离滤波(均值/min/max)
//! - [`terrain`]：坡度/起伏度/TPI
//! - [`hydro`]：水文修正 + D8 流向/汇流/嵌套河网/HAND
//! - [`distance`]：精确欧氏距离变换(含最近源索引)
//! - [`segment`]：峰顶提取 + 分水岭分割 + 连通域
//! - [`geomorphon`]：geomorphon 地貌形态模式
//! - [`pipeline`]：参数模型 + 全流程编排

pub mod basin;
pub mod distance;
pub mod geomorphon;
pub mod error;
pub mod filter;
pub mod geotiff;
pub mod hydro;
pub mod input;
pub mod pipeline;
pub mod postprocess;
pub mod ridge;
pub mod scale;
pub mod slope_position;
pub mod slope_unit;
pub mod segment;
pub mod terrain;
