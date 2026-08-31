# TerraPos — Terrain Position Classification Tool

[简体中文](README.md) | [English](README_EN.md)

Automatically classifies DEM data into 8 terrain positions for southwest China:
mountain-basin, hills (upper/middle/lower), and mountain slopes (upper/middle/lower).

## Terrain Position Codes

| Code | Position | Criteria Summary |
|---|---|---|
| 1 | Mountain basin | Slope <6° + TPI <−25 m (2 km) + Area ≥0.5 km² + Width >250 m |
| 2 | Broad-valley basin | Reserved code |
| 3/4/5 | Hills upper/middle/lower | Elevation <500 m, tertile position within landform unit |
| 6/7/8 | Mountain slopes upper/middle/lower | Elevation ≥500 m, tertile position within landform unit |

Geomorphological subclass raster: low hill / high hill / low mountain /
middle mountain / high mountain / extremely high mountain / flat basin.

## Usage

```bash
cd rust
cargo build --release
./target/release/topo_app.exe
```

Batch processing via CLI:

```bash
cargo run --release -p topo_core --example run_full -- <dem.tif> <out_dir>
```

UI workflow: import DEM (elevation / hillshade preview shown immediately) →
adjust parameters (defaults work out of the box) → run → preview & export.
Parameter meanings are documented in the in-app "参数详解" (Parameter Guide) page.

Five landform-unit segmentation strategies are available (`SeedMode`, hybrid
prominence ∪ distance by default); see the in-app guide and
`docs/seed_mode_comparison.png` for a strategy comparison.

## Algorithm Highlights

1. Priority-Flood fill (z-limit preserves karst deep sinks) → D8 flow accumulation;
2. Exact Euclidean distance transform → HAND terrace/foot-slope correction;
3. Peak-based watershed landform units; tertile position within each unit;
4. Per-pixel width core-normalization for basins (narrow valley bands → lower slope);
5. Hill/mountain split at 500 m elevation.

## Release

The portable package (`TerraPos-v0.0.1-win64.zip`) is available under
[Releases](https://github.com/runker54/TerraPos/releases) — unzip and run,
no installation or runtime dependencies required.
Release notes: `docs/RELEASE_NOTES-v0.0.1.md`.
