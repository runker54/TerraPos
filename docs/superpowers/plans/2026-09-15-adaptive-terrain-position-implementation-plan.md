# Adaptive Terrain Position Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the current fixed-window and global-ratio workflow with a DEM-only, multiscale, hydrologically connected terrain-position classifier that emits the approved seven effective classes (code 2 reserved and never emitted), preserves NoData, and produces auditable diagnostic rasters.

**Architecture:** Keep `topo_core` as the production engine and the desktop app as a thin caller. Build a dual-surface DEM model, derive a coarse multiscale hydrology–ridge–slope-unit graph, calculate adaptive relative slope position, classify fuzzy upper/middle/lower memberships, detect basin objects independently, then compose and clean the full-resolution result under terrain barriers. Python scripts are independent reference/validation tools with constants at file top and no command-line parameters.

**Tech Stack:** Rust 1.94.1, Cargo, `tiff`, `rayon`, `serde`, existing `topo_core`; Python 3.10.11 at `D:/worker_code/.venvgis/Scripts/python.exe` with GDAL, NumPy, SciPy, matplotlib, and Pillow already present; GeoTIFF in projected metre CRS.

## Global Constraints

- The approved design in `docs/superpowers/specs/2026-09-15-adaptive-terrain-position-design.md` is the source of truth.
- Required input is one Float32 DEM. Do not introduce river, catchment, land-cover, or training-label inputs.
- Do not retain a legacy fixed-window execution branch. Replace the current production pipeline directly.
- Effective output codes are fixed: `0=NoData`, `1=山间/宽谷盆地`, `2=reserved`, `3=丘陵坡下`, `4=丘陵坡中`, `5=丘陵坡上`, `6=山地坡下`, `7=山地坡中`, `8=山地坡上`.
- Every distance, window, and area parameter is represented in metres or square metres and converted to pixels exactly once at the consuming module boundary.
- Full-resolution arrays are tiled or released stage-by-stage; persistent global analytical arrays live on the coarse grid. Peak memory on `data/dem.tif` must remain at or below 6 GiB.
- NoData connected to the DEM boundary is never interpolated and is restored as code 0 in every categorical output.
- Python scripts use top-of-file constants and must not add command-line arguments.
- Do not tune toward fixed output-area percentages. Distribution checks are invariants and warning signals, not optimization targets.
- Every user-visible parameter needs a sensitivity test proving that changing it changes the intended intermediate layer or result.
- After any source or documentation change, run the complete project quality, packaging, Git, and GitHub Release sequence in Task 16. Do not publish an intermediate failing state.

---

## Target File Map

### Rust production code

| Path | Responsibility |
|---|---|
| `rust/topo_core/src/input.rs` | DEM metadata validation, valid mask, small interior-hole fill, geomorphic surface |
| `rust/topo_core/src/scale.rs` | metre-scale family, robust multiscale deviation/relief, characteristic-scale selection |
| `rust/topo_core/src/hydro.rs` | bounded conditioning, D8 route, accumulation, nested streams, flow-connected HAND |
| `rust/topo_core/src/ridge.rs` | subcatchment divides plus high-confidence topographic ridge evidence |
| `rust/topo_core/src/slope_unit.rs` | slope-unit labels, constrained valley/ridge distances, adaptive analysis scale and relative position |
| `rust/topo_core/src/geomorphon.rs` | ten-form geomorphon lookup and adaptive-scale morphology evidence |
| `rust/topo_core/src/slope_position.rs` | fuzzy upper/middle/lower memberships and confidence |
| `rust/topo_core/src/basin.rs` | basin candidate, width core, object tests, constrained reconstruction |
| `rust/topo_core/src/postprocess.rs` | final code composition, terrain-constrained patch merge, monotonic correction |
| `rust/topo_core/src/pipeline.rs` | configuration, orchestration, cancellation/progress, GeoTIFF/report output only |
| `rust/topo_core/src/geotiff.rs` | NoData tag and projected-metre CRS inspection/writing |
| `rust/topo_core/src/filter.rs` | reusable percentile/robust neighbourhood primitives |
| `rust/topo_core/src/lib.rs` | module exports |
| `rust/topo_app/src/main.rs` | compact basic controls, advanced panel, diagnostics preview |

### Rust tests and fixtures

| Path | Responsibility |
|---|---|
| `rust/topo_core/tests/common/mod.rs` | deterministic synthetic DEM builders and assertion helpers |
| `rust/topo_core/tests/input_contract.rs` | CRS, NoData, interior-hole behavior |
| `rust/topo_core/tests/hydrology.rs` | routing, nested streams, true HAND |
| `rust/topo_core/tests/adaptive_scale.rs` | characteristic scale and resolution stability |
| `rust/topo_core/tests/geomorphon_forms.rs` | all ten geomorphon forms and metre/pixel conversion regression |
| `rust/topo_core/tests/slope_geometry.rs` | ridges, slope units, constrained distances, relative position |
| `rust/topo_core/tests/basin_objects.rs` | valley rejection and complete basin reconstruction |
| `rust/topo_core/tests/classification.rs` | code domain, class ordering, transformations and robustness |
| `rust/topo_core/tests/parameter_effects.rs` | every exposed parameter has an observable intended effect |

### Python and documentation

| Path | Responsibility |
|---|---|
| `python/build_synthetic_dem.py` | write reproducible scene GeoTIFFs and expected-zone masks |
| `python/adaptive_terrain_reference.py` | readable reference calculations for selected small scenes |
| `python/validate_terrain_result.py` | acceptance metrics for real and synthetic runs |
| `python/render_qa_atlas.py` | fixed-layout PNG diagnostic atlas |
| `python/README.md` | script inputs, constants, outputs, interpretation |
| `README.md` | product workflow and seven-class semantics |
| `docs/RELEASE_NOTES-v0.0.1.md` | implemented algorithm, validation evidence, limitations |

---

## Shared Data Contracts

Implement these names and meanings before module work diverges:

```rust
// input.rs
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RasterShape {
    pub width: usize,
    pub height: usize,
    pub resolution_m: f64,
}

#[derive(Debug)]
pub struct PreparedDem {
    pub raw: Vec<f32>,
    pub geomorph: Vec<f32>,
    pub valid: Vec<bool>,
    pub shape: RasterShape,
    pub meta: GeoMeta,
}

// pipeline.rs
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PrecisionPreset { Fast, Standard, Detailed }

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BasinTendency { Strict, Standard, Loose }

#[derive(Debug, Clone)]
pub struct Params {
    pub dem_path: String,
    pub out_dir: String,
    pub precision: PrecisionPreset,
    pub basin_tendency: BasinTendency,
    pub basin_min_area_m2: f64,
    pub postprocess_strength: f64,
    pub write_diagnostics: bool,
    pub advanced: AdvancedParams,
}

#[derive(Debug, Clone)]
pub struct AdvancedParams {
    pub coarse_res_m: f64,
    pub hydro_z_limit_m: f32,
    pub scale_growth_threshold: f32,
    pub low_relief_m: f32,
    pub hill_elevation_max_m: f32,
    pub stream_areas_km2: [f64; 4],
}
```

Default values are `Standard`, `Standard`, `66_666.67 m²`, `1.0`, diagnostics enabled, and advanced values `25 m`, `15 m`, `0.15`, `20 m`, `500 m`, `[0.05, 0.20, 1.00, 5.00] km²`.

The production pipeline consumes and produces flat row-major buffers. Use `u32::MAX` as the sole “no downstream/source cell” sentinel. Categorical diagnostics use 0 for invalid/unassigned; continuous diagnostics use Float32 NoData.

---

## Task 1: Freeze Synthetic Fixtures and Classification Invariants

**Files:**

- Create: `rust/topo_core/tests/common/mod.rs`
- Create: `rust/topo_core/tests/classification.rs`
- Modify: `rust/topo_core/examples/run_sample.rs`

- [x] **Step 1: Add deterministic in-memory terrain builders**

Write builders with no random component: `planar_slope`, `v_valley`, `broad_basin`, `closed_pit`, `conical_hill`, `nested_ridges`, and `with_border_nodata`. Each returns `(Vec<f32>, RasterShape)` and uses centre-of-cell metre coordinates.

```rust
pub fn broad_basin(width: usize, height: usize, res: f64) -> (Vec<f32>, RasterShape) {
    let cx = width as f64 * res * 0.5;
    let cy = height as f64 * res * 0.5;
    let z = (0..height).flat_map(|row| {
        (0..width).map(move |col| {
            let x = (col as f64 + 0.5) * res - cx;
            let y = (row as f64 + 0.5) * res - cy;
            let r = x.hypot(y);
            let floor = 800.0 + 0.002 * x;
            let rim = ((r - 900.0) / 250.0).clamp(0.0, 1.0);
            (floor + 120.0 * rim * rim) as f32
        })
    }).collect();
    (z, RasterShape { width, height, resolution_m: res })
}
```

- [x] **Step 2: Add failing end-state invariant tests**

The first test calls the future array-level entry point `run_arrays_for_test`; it must fail to compile until Task 12 exposes that test API. Assert valid codes are only `{1,3,4,5,6,7,8}`, code 2 count is zero, all valid cells are nonzero, invalid cells are zero, and every planar ridge-to-valley trace is non-increasing in the order upper→middle→lower.

```rust
#[test]
fn final_codes_obey_contract() {
    let (dem, shape) = common::nested_ridges(257, 257, 10.0);
    let valid = vec![true; dem.len()];
    let out = topo_core::pipeline::run_arrays_for_test(&dem, &valid, shape, &Default::default()).unwrap();
    assert!(out.terrain.iter().all(|c| matches!(c, 1 | 3 | 4 | 5 | 6 | 7 | 8)));
    assert!(!out.terrain.contains(&2));
}
```

- [x] **Step 3: Run the test and record the expected failure**

Run: `cd rust; cargo test --release --test classification`

Expected: compile failure stating `run_arrays_for_test` does not exist. This proves the end-state test is connected to the new API rather than the legacy pipeline.

- [x] **Step 4: Turn the sample example into a smoke gate without target ratios**

Keep the printed statistics, remove “8类” wording, and add assertions: code 2 absent; at least one upper, middle, and lower class exists; no non-basin class exceeds 90% of valid pixels; NoData is excluded from the denominator. Do not require a basin because a sample tile may contain no complete basin object.

- [x] **Step 5: Commit the fixture contract**

```powershell
git add rust/topo_core/tests/common/mod.rs rust/topo_core/tests/classification.rs rust/topo_core/examples/run_sample.rs
git commit -m "test: define adaptive terrain classification invariants"
```

---

## Task 2: Enforce DEM Metadata and NoData Contracts

**Files:**

- Create: `rust/topo_core/src/input.rs`
- Create: `rust/topo_core/tests/input_contract.rs`
- Modify: `rust/topo_core/src/geotiff.rs`
- Modify: `rust/topo_core/src/lib.rs`
- Modify: `rust/topo_core/src/error.rs`

- [x] **Step 1: Write failing GeoTIFF metadata tests**

Add tests for: Float32 projected-metre raster accepted; geographic EPSG:4490 rejected; projected feet rejected; non-square or non-positive pixels rejected; GDAL NoData tag parsed; boundary-connected NoData preserved; a finite interior hole up to `10_000 m²` filled; a larger interior hole preserved invalid.

```rust
#[test]
fn rejects_geographic_crs_even_when_pixel_size_is_small() {
    let meta = common::meta_with_keys(32, 32, 0.00005, 2, 9102);
    let err = validate_meta(&meta).unwrap_err().to_string();
    assert!(err.contains("投影坐标系") && err.contains("米"));
}
```

- [x] **Step 2: Run the tests and verify failure**

Run: `cd rust; cargo test --release --test input_contract`

Expected: unresolved imports for `input::validate_meta` and missing `GeoMeta::nodata`.

- [x] **Step 3: Extend `GeoMeta` and GeoTIFF parsing**

Add `pub nodata: Option<f32>` to `GeoMeta`. Parse TIFF tag 42113 as trimmed ASCII Float32. Add GeoKey lookup helpers for `GTModelTypeGeoKey=1024` and `ProjLinearUnitsGeoKey=3076`; accept model type 1 with unit EPSG 9001 only. Preserve raw keys when writing and write tag 42113 for Float32 diagnostics and categorical outputs.

```rust
impl GeoMeta {
    pub fn geo_key_u16(&self, key_id: u16) -> Option<u16> {
        let entry_count = *self.geo_keys.get(3)? as usize;
        (0..entry_count).find_map(|i| {
            let base = 4 + i * 4;
            let entry = self.geo_keys.get(base..base + 4)?;
            (entry[0] == key_id && entry[1] == 0 && entry[2] == 1)
                .then_some(entry[3])
        })
    }
    pub fn is_projected_metre(&self) -> bool {
        self.geo_key_u16(1024) == Some(1) && self.geo_key_u16(3076) == Some(9001)
    }
}
```

The implementation of `geo_key_u16` must honor the GeoKeyDirectory header count and only return inline values where `TIFFTagLocation == 0`; unsupported indirection yields `None` and therefore an explicit validation error.

- [x] **Step 4: Implement input preparation**

Expose:

```rust
pub struct InputConfig {
    pub max_interior_hole_area_m2: f64,
    pub geomorph_smooth_radius_m: f64,
}

pub fn validate_meta(meta: &GeoMeta) -> Result<RasterShape>;
pub fn prepare_values(values: Vec<f32>, meta: GeoMeta, cfg: &InputConfig) -> Result<PreparedDem>;
pub fn prepare_input(path: &Path, cfg: &InputConfig) -> Result<PreparedDem>;
```

Build `valid` from finite values not equal to metadata NoData. Flood-fill invalid cells from all raster edges to mark external NoData. Connected internal invalid objects are filled by nearest valid cell only when `cell_count * resolution_m² <= max_interior_hole_area_m2`; larger objects stay invalid. Produce `geomorph` by one valid-aware circular mean filter with the requested metre radius. Never mutate `raw` outside accepted small interior holes.

- [x] **Step 5: Add actionable errors and export the module**

Use `CoreError::Invalid` messages containing the observed pixel size, model type, and unit key. Add `pub mod input;` to `lib.rs` and correct the stale module list in its crate docs.

- [x] **Step 6: Run tests**

Run: `cd rust; cargo test --release --test input_contract`

Expected: all input contract tests pass.

- [x] **Step 7: Commit**

```powershell
git add rust/topo_core/src/input.rs rust/topo_core/src/geotiff.rs rust/topo_core/src/lib.rs rust/topo_core/src/error.rs rust/topo_core/tests/input_contract.rs
git commit -m "feat: validate projected DEM input and preserve nodata"
```

---

## Task 3: Build Bounded Hydrological Conditioning, Nested Streams, and True HAND

**Files:**

- Modify: `rust/topo_core/src/hydro.rs`
- Create: `rust/topo_core/tests/hydrology.rs`

- [x] **Step 1: Write failing route and HAND tests**

Cover: every routed cell reaches an outlet or recorded deep sink in at most `n` steps; accumulation is non-decreasing downstream; stream masks are nested from 0.05 to 5.00 km²; HAND is exactly zero on stream cells; an upslope cell receives elevation relative to the first stream encountered along `flow_to`, not the Euclidean-nearest stream; conditioning depth never exceeds the z-limit.

```rust
#[test]
fn hand_follows_flow_connected_stream() {
    let dem = vec![110., 109., 108., 107., 106., 105., 104., 103., 102.];
    let flow_to = vec![1, 2, 3, 4, 5, 6, 7, 8, u32::MAX];
    let stream = vec![false, false, false, false, false, false, true, false, true];
    let hand = hand_to_stream(&dem, &flow_to, &stream, &[true; 9]).unwrap();
    assert_eq!(hand[0], 6.0);
}
```

- [x] **Step 2: Run and confirm failure**

Run: `cd rust; cargo test --release --test hydrology`

Expected: missing `HydroConfig`, `HydroModel`, `build_hydro`, and `hand_to_stream`.

- [x] **Step 3: Replace the monolithic hydrology result with explicit products**

```rust
#[derive(Debug, Clone)]
pub struct HydroConfig {
    pub coarse_res_m: f64,
    pub z_limit_m: f32,
    pub stream_areas_km2: [f64; 4],
}

pub struct HydroModel {
    pub conditioned: Vec<f32>,
    pub conditioning_depth: Vec<f32>,
    pub flow_to: Vec<u32>,
    pub pop_order: Vec<u32>,
    pub accumulation_cells: Vec<u32>,
    pub deep_sink: Vec<bool>,
    pub streams: [Vec<bool>; 4],
    pub stream_level: Vec<u8>,
    pub shape: RasterShape,
}

pub fn build_hydro(prepared: &PreparedDem, cfg: &HydroConfig) -> Result<HydroModel>;
pub fn hand_to_stream(dem: &[f32], flow_to: &[u32], stream: &[bool], valid: &[bool]) -> Result<Vec<f32>>;
```

Downsample the hydrological surface to `coarse_res_m` with a valid-aware block minimum/mean hybrid: use the block minimum as drainage support and clamp it no lower than `block_mean - z_limit_m`. Apply bounded Priority-Flood. Cells requiring more than `z_limit_m` of raising stay `deep_sink=true`; they are allowed outlets, not flattened basins.

- [x] **Step 4: Derive stream thresholds from physical area**

For each threshold compute `ceil(area_km2 * 1_000_000 / coarse_cell_area_m2)` once. Sort and validate strictly increasing configuration. `stream_level` is 0 off-stream and 1–4 for the finest-to-coarsest nested level a cell reaches. Reject a threshold that maps below one cell or exceeds valid area.

- [x] **Step 5: Implement memoized downstream HAND**

Trace downstream until reaching the first stream cell, outlet, deep sink, or already-solved cell; unwind the stack once. On stream, HAND is 0. At an outlet/deep sink without a stream encounter, HAND is Float32 NoData. Clamp small negative numerical values to zero but return an error when a negative value is below `-0.05 m`, because that indicates inconsistent surfaces or routing.

- [x] **Step 6: Run tests and existing hydro users**

Run: `cd rust; cargo test --release --test hydrology; cargo test --release hydro::`

Expected: all hydrology tests pass; compile errors in legacy pipeline are acceptable until Task 12, but library unit tests must compile by temporarily retaining only the existing `fill_and_route` wrapper if another untouched example still imports it. This wrapper must be deleted in Task 12, not shipped as a compatibility path.

- [x] **Step 7: Commit**

```powershell
git add rust/topo_core/src/hydro.rs rust/topo_core/tests/hydrology.rs
git commit -m "feat: derive nested streams and flow-connected hand"
```

---

## Task 4: Implement the Robust Multiscale Terrain Pyramid

**Files:**

- Create: `rust/topo_core/src/scale.rs`
- Create: `rust/topo_core/tests/adaptive_scale.rs`
- Modify: `rust/topo_core/src/filter.rs`
- Modify: `rust/topo_core/src/lib.rs`

- [x] **Step 1: Write failing characteristic-scale tests**

Use geometrically identical Gaussian hills sampled at 5 m, 10 m, and 25 m. Assert median selected `R0` over the central landform differs by at most one adjacent scale step. Add a narrow hill and broad mountain test and assert the broad feature selects a larger median scale. Add a constant vertical offset and assert identical selected scales.

- [x] **Step 2: Run and confirm failure**

Run: `cd rust; cargo test --release --test adaptive_scale`

Expected: missing `scale` module.

- [x] **Step 3: Add valid-aware robust filters**

Implement `focal_quantile_valid`, `focal_median_valid`, and `focal_mad_valid`. Use an exact sorted neighbourhood only in unit tests and a separable/histogram or tiled selection implementation in production so a 4000 m radius never creates an `O(n*r²)` loop. The functions accept metre radius and convert using `RasterShape::resolution_m` internally once.

- [x] **Step 4: Add scale contracts and bounds**

```rust
pub const BASE_SCALES_M: [f64; 6] = [125.0, 250.0, 500.0, 1000.0, 2000.0, 4000.0];

pub struct ScaleLayer {
    pub radius_m: f64,
    pub deviation_m: Vec<f32>,
    pub relief_m: Vec<f32>,
    pub normalized_deviation: Vec<f32>,
}

pub struct ScalePyramid {
    pub layers: Vec<ScaleLayer>,
    pub characteristic_scale_m: Vec<f32>,
    pub characteristic_relief_m: Vec<f32>,
    pub shape: RasterShape,
}

pub fn usable_scales(shape: RasterShape) -> Vec<f64>;
pub fn build_scale_pyramid(dem: &[f32], valid: &[bool], shape: RasterShape, growth_threshold: f32) -> Result<ScalePyramid>;
```

Drop scales smaller than five coarse pixels and larger than one quarter of the shorter valid-extent dimension. Require at least three usable scales; otherwise return a clear “DEM extent too small for multiscale classification” error.

- [x] **Step 5: Implement scale metrics and selection**

At each scale compute `DEV=z-local_median`, `H=P95-P05`, `NDEV=DEV/(1.4826*MAD+0.1)`, and `g=(H_next-H)/max(H_next,0.1)`. Select the first scale with `g < threshold` for two consecutive transitions and `H > 5 m`; if none qualifies, select the scale minimizing `abs(g-threshold)` among layers with `H > 5 m`, else the smallest usable scale. Store only coarse-grid layers.

- [x] **Step 6: Add landform scale clamping helper**

Expose `scale_bounds_for(elevation_m, relief_m, data_max_m) -> (f32, f32)` with exact approved ranges: low hill 125–750, high hill 250–1500, low mountain 250–2000, middle mountain 500–4000, high mountain 750–6000, extreme mountain 1000–min(8000,data maximum). Apply the project’s current elevation/relief subclass thresholds only here, so classification code does not duplicate them.

- [x] **Step 7: Run tests and commit**

Run: `cd rust; cargo test --release --test adaptive_scale`

Expected: all adaptive-scale tests pass.

```powershell
git add rust/topo_core/src/scale.rs rust/topo_core/src/filter.rs rust/topo_core/src/lib.rs rust/topo_core/tests/adaptive_scale.rs
git commit -m "feat: select terrain scale from robust multiscale relief"
```

---

## Task 5: Replace Count-Based Geomorphons with the Ten Standard Forms

**Files:**

- Modify: `rust/topo_core/src/geomorphon.rs`
- Modify: `rust/topo_core/src/terrain.rs`
- Create: `rust/topo_core/tests/geomorphon_forms.rs`

- [x] **Step 1: Write failing form-table and distance tests**

Construct representative `Pattern8` values for flat, peak, ridge, shoulder, spur, slope, hollow, footslope, valley, and pit. Assert each maps to its own form. Add a regression where a 1000 m search at 5 m resolution reaches 200 cells, not 40 cells.

```rust
#[repr(u8)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Landform {
    Flat = 1, Peak = 2, Ridge = 3, Shoulder = 4, Spur = 5,
    Slope = 6, Hollow = 7, Footslope = 8, Valley = 9, Pit = 10,
}
```

- [x] **Step 2: Run and confirm the legacy mapping fails**

Run: `cd rust; cargo test --release --test geomorphon_forms`

Expected: missing `Landform` or multiple forms collapsing into the same old level.

- [x] **Step 3: Implement canonical rotational lookup**

Encode each pattern as trits `{lower=0, flat=1, higher=2}`. Generate its eight rotations and eight reflected rotations, choose the minimum canonical code, and map the canonical code to the standard ten-form decision table. Keep `Pattern8` public for diagnostics. Remove `pattern_to_level` and `geomorphon_levels` after all callers move in Task 12.

- [x] **Step 4: Correct the unit boundary**

`geomorphon_pattern` continues to accept `search_m` and `skip_m`. It alone calculates `ceil(distance_m / resolution_m)`. All callers pass metre values unchanged. Add debug assertions that `skip_m >= resolution_m`, `search_m > skip_m`, and generated sample offsets do not exceed `ceil(search_m/resolution_m)`.

- [x] **Step 5: Add terrain derivatives**

Expose valid-aware `profile_curvature` and `plan_curvature` in `terrain.rs` using 3×3 quadratic finite differences. The edge and invalid-neighbour policy is Float32 NoData, not replicated values.

- [x] **Step 6: Run and commit**

Run: `cd rust; cargo test --release --test geomorphon_forms; cargo test --release terrain::`

Expected: ten forms and unit regression pass.

```powershell
git add rust/topo_core/src/geomorphon.rs rust/topo_core/src/terrain.rs rust/topo_core/tests/geomorphon_forms.rs
git commit -m "feat: implement ten-form geomorphon evidence"
```

---

## Task 6: Extract Ridges and Label Slope Units

**Files:**

- Create: `rust/topo_core/src/ridge.rs`
- Create: `rust/topo_core/src/slope_unit.rs`
- Create: `rust/topo_core/tests/slope_geometry.rs`
- Modify: `rust/topo_core/src/lib.rs`
- Reuse: `rust/topo_core/src/segment.rs`

- [x] **Step 1: Write failing ridge and unit tests**

On a symmetric V-valley with two flanking ridges, assert the central valley is not marked ridge, each flank contains a continuous ridge, no slope unit crosses a ridge or stream barrier, and every non-flat valid cell has one reachable valley and ridge anchor. On a mirrored DEM, unit areas and distance distributions must match within one coarse cell.

- [x] **Step 2: Run and confirm failure**

Run: `cd rust; cargo test --release --test slope_geometry`

Expected: unresolved `ridge` and `slope_unit` modules.

- [x] **Step 3: Implement ridge evidence**

```rust
pub struct RidgeModel {
    pub mask: Vec<bool>,
    pub strength: Vec<f32>,
    pub subcatchment_id: Vec<u32>,
    pub shape: RasterShape,
}

pub fn build_ridges(
    dem: &[f32],
    valid: &[bool],
    hydro: &HydroModel,
    pyramid: &ScalePyramid,
    landform: &[Landform],
) -> Result<RidgeModel>;
```

Label each fine-stream subcatchment by reverse traversal of `flow_to`. Mark contacts between different labels as primary divide candidates. Score supplements from inverse-terrain flow accumulation, positive normalized deviation, and `Peak|Ridge|Shoulder|Spur`. Retain supplements only when strength ≥0.7 and connected to a primary divide within `0.1 * characteristic_scale`, capped at 250 m. Thin the combined mask to one coarse cell without breaking connectivity.

- [x] **Step 4: Label slope units**

```rust
pub struct SlopeUnits {
    pub unit_id: Vec<u32>,
    pub aspect_sector: Vec<u8>,
    pub valley_mask: Vec<bool>,
    pub ridge_mask: Vec<bool>,
    pub shape: RasterShape,
}

pub fn build_slope_units(
    dem: &[f32],
    valid: &[bool],
    hydro: &HydroModel,
    ridges: &RidgeModel,
) -> Result<SlopeUnits>;
```

Use the finest valid stream mask as valley barrier, ridge mask as crest barrier, and eight aspect sectors as secondary partitions. Seed each connected region between barriers and flood only to neighbours whose circular aspect difference is ≤90°. Merge regions smaller than four coarse cells into the neighbour with longest non-barrier boundary.

- [x] **Step 5: Run and commit**

Run: `cd rust; cargo test --release --test slope_geometry ridge_continuity slope_units_do_not_cross_barriers`

Expected: named tests pass.

```powershell
git add rust/topo_core/src/ridge.rs rust/topo_core/src/slope_unit.rs rust/topo_core/src/lib.rs rust/topo_core/tests/slope_geometry.rs
git commit -m "feat: derive ridge-bounded slope units"
```

---

## Task 7: Compute Constrained Distances, Adaptive Scale, and Relative Position

**Files:**

- Modify: `rust/topo_core/src/slope_unit.rs`
- Modify: `rust/topo_core/tests/slope_geometry.rs`

- [ ] **Step 1: Add failing constrained-distance tests**

Create two adjacent valleys separated by a high ridge. Assert the Euclidean-nearest valley across the ridge is never selected. Assert `dv=0` on valley, `dr=0` on ridge, `L=dv+dr`, all finite `q` values lie in `[0,1]`, and low-relief areas weight horizontal position more strongly than vertical position.

- [ ] **Step 2: Run and verify failure**

Run: `cd rust; cargo test --release --test slope_geometry constrained_distance`

Expected: missing `SlopeGeometry` or assertion failure from unconstrained distance.

- [ ] **Step 3: Add exact geometry output**

```rust
pub struct SlopeGeometry {
    pub valley_anchor: Vec<u32>,
    pub ridge_anchor: Vec<u32>,
    pub distance_to_valley_m: Vec<f32>,
    pub distance_to_ridge_m: Vec<f32>,
    pub local_width_m: Vec<f32>,
    pub adaptive_scale_m: Vec<f32>,
    pub hand_m: Vec<f32>,
    pub q_distance: Vec<f32>,
    pub q_elevation: Vec<f32>,
    pub relative_position: Vec<f32>,
    pub shape: RasterShape,
}

pub fn build_slope_geometry(
    dem: &[f32], valid: &[bool], units: &SlopeUnits, hydro: &HydroModel,
    pyramid: &ScalePyramid, low_relief_m: f32,
) -> Result<SlopeGeometry>;
```

- [ ] **Step 4: Implement barrier-aware multi-source Dijkstra**

Use eight-neighbour step length and edge cost:

```text
cost = delta_s * (1 + 2*abs(sin(delta_aspect)) + 4*cross_barrier)
```

Run one pass from valley seeds and one from ridge seeds. Do not enqueue a neighbour with a different nonzero `unit_id`; therefore `cross_barrier` is normally zero and remains four only for a one-cell gap repair inside the same unit. Store the originating anchor index with each shortest path.

- [ ] **Step 5: Select the local stream/analysis scale**

For each stream level `k`, compute its constrained valley/ridge width `L_k`. Choose the level minimizing:

```text
abs(ln(L_k / (2*R0))) + 0.5*abs(ln(L_(k+1) / L_k))
```

Use the last term as zero at the coarsest level. Then set `R*=clip(0.5*L, Rmin, Rmax)` with `Rmin/Rmax` from `scale_bounds_for`. Record the chosen level for diagnostics.

- [ ] **Step 6: Compute true relative position**

Sample valley and ridge elevations at stored anchors. Compute:

```text
qd = dv / (dv + dr + 0.001)
qz = clamp((z - zv) / (zr - zv + 0.001), 0, 1)
w_distance = lerp(0.45, 0.75, clamp((20 - H*) / 20, 0, 1))
q = clamp(w_distance*qd + (1-w_distance)*qz, 0, 1)
```

If either anchor is absent, use the available component only and set a low-confidence flag consumed by Task 9; do not fabricate an anchor with Euclidean nearest-neighbour lookup.

- [ ] **Step 7: Run tests and commit**

Run: `cd rust; cargo test --release --test slope_geometry`

Expected: all ridge, unit, distance, scale, and relative-position tests pass.

```powershell
git add rust/topo_core/src/slope_unit.rs rust/topo_core/tests/slope_geometry.rs
git commit -m "feat: compute adaptive relative slope geometry"
```

---

## Task 8: Calculate Adaptive Morphology Evidence

**Files:**

- Modify: `rust/topo_core/src/geomorphon.rs`
- Modify: `rust/topo_core/tests/geomorphon_forms.rs`

- [x] **Step 1: Write failing adaptive-evidence tests**

For a hill-to-valley transect, assert peak/ridge/shoulder evidence favours upper position, hollow/footslope/valley favours lower position, and planar slope favours middle position. Double the feature width and assert the sampled geomorphon radius follows `adaptive_scale_m` rather than staying fixed.

- [x] **Step 2: Run and verify failure**

Run: `cd rust; cargo test --release --test geomorphon_forms adaptive_evidence`

Expected: missing `MorphEvidence` and `adaptive_morphology_evidence`.

- [x] **Step 3: Implement scale-indexed morphology layers**

```rust
pub struct MorphEvidence {
    pub form: Vec<Landform>,
    pub upper: Vec<f32>,
    pub middle: Vec<f32>,
    pub lower: Vec<f32>,
    pub slope_deg: Vec<f32>,
    pub profile_curvature: Vec<f32>,
    pub plan_curvature: Vec<f32>,
    pub low_confidence: Vec<bool>,
}

pub fn adaptive_morphology_evidence(
    dem: &[f32], valid: &[bool], shape: RasterShape,
    adaptive_scale_m: &[f32], usable_scales_m: &[f64],
) -> Result<MorphEvidence>;
```

Compute one geomorphon raster per usable scale using `search_m=scale` and `skip_m=max(2*resolution, 0.05*scale)`. For each cell take the form from its nearest scale layer. Use these exact base contributions:

| Form | upper | middle | lower |
|---|---:|---:|---:|
| Peak | 1.00 | 0.00 | 0.00 |
| Ridge | 0.90 | 0.10 | 0.00 |
| Shoulder | 0.70 | 0.30 | 0.00 |
| Spur | 0.55 | 0.45 | 0.00 |
| Slope | 0.10 | 0.80 | 0.10 |
| Hollow | 0.00 | 0.45 | 0.55 |
| Footslope | 0.00 | 0.25 | 0.75 |
| Valley | 0.00 | 0.10 | 0.90 |
| Pit | 0.00 | 0.00 | 1.00 |
| Flat | 0.05 | 0.25 | 0.70 |

Refine each score by at most 0.10 from the sign and robust magnitude of profile curvature and plan curvature. Renormalize the three scores to sum to one. Mark a cell low-confidence when the selected scale is at a scale-family endpoint or derivative support is incomplete.

- [x] **Step 4: Limit memory by releasing scale layers**

Process scale layers from small to large. Copy a cell’s form into the final buffer when that layer is its nearest selected scale, then release the layer. Keep only the final form/evidence and current working buffers.

- [x] **Step 5: Run and commit**

Run: `cd rust; cargo test --release --test geomorphon_forms`

Expected: all ten-form and adaptive-evidence tests pass.

```powershell
git add rust/topo_core/src/geomorphon.rs rust/topo_core/tests/geomorphon_forms.rs
git commit -m "feat: sample morphology at local terrain scale"
```

---

## Task 9: Classify Fuzzy Upper, Middle, and Lower Slope Positions

**Files:**

- Create: `rust/topo_core/src/slope_position.rs`
- Modify: `rust/topo_core/src/lib.rs`
- Modify: `rust/topo_core/tests/classification.rs`

- [x] **Step 1: Write failing membership tests**

Test exact ordering at `q=0.1`, `0.5`, and `0.9`; continuity either side of 0.38 and 0.62; tie behavior at the centre; morphology evidence resolving a near tie; confidence equal to the largest membership minus the second largest; low-confidence geometry reducing confidence without changing membership order.

- [x] **Step 2: Run and confirm failure**

Run: `cd rust; cargo test --release --test classification fuzzy_membership`

Expected: missing `slope_position` module.

- [x] **Step 3: Implement the fuzzy model**

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SlopePosition { Lower = 1, Middle = 2, Upper = 3 }

pub struct Memberships {
    pub upper: Vec<f32>,
    pub middle: Vec<f32>,
    pub lower: Vec<f32>,
    pub raw: Vec<SlopePosition>,
    pub confidence: Vec<f32>,
}

pub fn classify_slope_positions(
    relative_position: &[f32], morph: &MorphEvidence, valid: &[bool],
) -> Result<Memberships>;
```

Use:

```text
U0 = sigmoid((q-0.62)/0.08)
L0 = sigmoid((0.38-q)/0.08)
M0 = exp(-((q-0.50)/0.22)^2)
U = U0 + 0.20*morph.upper
M = M0 + 0.15*morph.middle
L = L0 + 0.20*morph.lower
```

Normalize `U/M/L`. Resolve exact ties in favour of middle, then lower, then upper; this prevents small floating-point changes from creating false crest pixels. Confidence is `max-second_max`; multiply it by 0.6 for `morph.low_confidence` or missing slope anchors.

- [x] **Step 4: Run and commit**

Run: `cd rust; cargo test --release --test classification fuzzy_membership`

Expected: all fuzzy membership tests pass.

```powershell
git add rust/topo_core/src/slope_position.rs rust/topo_core/src/lib.rs rust/topo_core/tests/classification.rs
git commit -m "feat: classify fuzzy upper middle and lower positions"
```

---

## Task 10: Detect Basin Objects and Reconstruct Their Full Boundary

**Files:**

- Create: `rust/topo_core/src/basin.rs`
- Create: `rust/topo_core/tests/basin_objects.rs`
- Modify: `rust/topo_core/src/lib.rs`

- [x] **Step 1: Write failing basin discrimination tests**

Use four scenes: broad enclosed basin, broad stream-connected valley basin, narrow V-valley, and flat upland. Assert the first two are accepted, the latter two rejected. For the accepted basin, assert the final mask includes at least 90% of the known flat-floor candidate while the width core occupies a smaller area. Assert increasing `basin_min_area_m2` removes only undersized objects.

- [x] **Step 2: Run and confirm failure**

Run: `cd rust; cargo test --release --test basin_objects`

Expected: unresolved `basin` module.

- [x] **Step 3: Define basin products and policy**

```rust
pub struct BasinConfig {
    pub tendency: BasinTendency,
    pub min_area_m2: f64,
}

pub struct BasinObject {
    pub id: u32,
    pub area_m2: f64,
    pub max_width_m: f32,
    pub median_width_m: f32,
    pub inner_relief_m: f32,
    pub surround_rise_m: f32,
    pub stream_connected: bool,
    pub closed_depression: bool,
    pub accepted: bool,
}

pub struct BasinResult {
    pub candidate: Vec<bool>,
    pub core: Vec<bool>,
    pub mask: Vec<bool>,
    pub objects: Vec<BasinObject>,
}

pub fn detect_basins(
    dem: &[f32], valid: &[bool], shape: RasterShape,
    hydro: &HydroModel, geometry: &SlopeGeometry,
    pyramid: &ScalePyramid, morph: &MorphEvidence,
    cfg: &BasinConfig,
) -> Result<BasinResult>;
```

- [x] **Step 4: Build the pixel candidate with adaptive thresholds**

Candidate requires lower geomorphic evidence; stream connectivity or recorded closed depression; `q<0.45`; low local slope, relief, and HAND; and low-flat persistence in at least three adjacent scale layers. For each cell compute:

```text
r_micro = clip(0.05*L, 25, 100) m
Hmicro = clip(0.03*H*, 3, 12) m
HANDmax = clip(0.02*H* + 2, 5, 20) m
```

Slope limits by relief/elevation subclass are exactly `6, 5.5, 5, 5, 4, 4` degrees from low hill through extreme mountain. Map `Strict/Standard/Loose` to multiplicative factors `0.85/1.00/1.15` on `Hmicro`, `HANDmax`, and slope limit only; do not change object area or width rules.

- [x] **Step 5: Evaluate connected objects**

For every candidate object derive:

```text
Wmin = clip(0.12*median(L), 150, 500) m
Amin = max(configured_min_area, clip(0.015*median(L)^2, 50_000, 500_000)) m²
Hinner = clip(0.08*median(H*), 10, 30) m
ring_radius = clip(0.15*median(L), 150, 500) m
required_surround = max(10, 0.05*median(H*)) m
```

Measure width with Euclidean distance to the candidate boundary: `width=2*distance`. Require area ≥Amin, maximum width ≥Wmin, median width ≥0.5Wmin, `P95-P05 ≤ Hinner`, ring median elevation rise ≥required surround, hydrological connection or closed-depression flag, and three-scale persistence.

- [x] **Step 6: Reconstruct the complete basin mask**

Define core cells as candidate width ≥Wmin. Reject objects with no core. Starting from accepted cores, perform geodesic reconstruction constrained to the original candidate component and stop at ridge barriers. The output `mask` is the reconstructed candidate, not the eroded width core. Fill only interior holes smaller than `min(0.1*object_area, 20_000 m²)`.

- [x] **Step 7: Run and commit**

Run: `cd rust; cargo test --release --test basin_objects`

Expected: accepted broad basins, rejected narrow valleys/upland flats, and full-boundary reconstruction assertions pass.

```powershell
git add rust/topo_core/src/basin.rs rust/topo_core/src/lib.rs rust/topo_core/tests/basin_objects.rs
git commit -m "feat: classify and reconstruct geomorphic basin objects"
```

---

## Task 11: Compose Codes and Apply Terrain-Constrained Cleanup

**Files:**

- Create: `rust/topo_core/src/postprocess.rs`
- Modify: `rust/topo_core/src/lib.rs`
- Modify: `rust/topo_core/tests/classification.rs`

- [x] **Step 1: Write failing composition and cleanup tests**

Assert basin overrides slope position; hill/mountain is selected solely by the configured 500 m elevation boundary; code 2 is never emitted; invalid cells remain zero. Add a one-cell speckle adjacent to a ridge and assert cleanup cannot move it across the ridge. Add 500 deterministic ridge-to-valley traces and require at least 95% to be monotonic after correction.

- [x] **Step 2: Run and verify failure**

Run: `cd rust; cargo test --release --test classification compose_codes`

Expected: missing `postprocess` module.

- [x] **Step 3: Implement direct composition**

```rust
pub struct FinalClassification {
    pub terrain: Vec<u8>,
    pub geomorph_subclass: Vec<u8>,
    pub confidence: Vec<f32>,
}

pub fn compose_codes(
    elevation_m: &[f32], valid: &[bool], basin: &BasinResult,
    positions: &Memberships, hill_elevation_max_m: f32,
) -> Result<FinalClassification>;

pub fn constrained_cleanup(
    result: &mut FinalClassification, positions: &Memberships,
    units: &SlopeUnits, geometry: &SlopeGeometry, basin: &BasinResult,
    shape: RasterShape, strength: f64,
) -> Result<()>;
```

Map lower/middle/upper to codes `3/4/5` below 500 m and `6/7/8` at or above 500 m. `geomorph_subclass` stores the existing relief/elevation subclass code independently of slope-position output.

- [x] **Step 4: Merge small patches without crossing terrain barriers**

Set smoothing distance `clip(5*strength, 5, 50) m`, further capped by `0.02*median(L)` per unit. Label same-code patches inside each slope unit and hill/mountain zone. A patch below `smoothing_distance²` merges into the adjacent class with longest shared boundary, then smallest total membership loss. Never merge across ridge, valley, basin edge, NoData, or hill/mountain boundary.

- [x] **Step 5: Enforce monotonic traces with dynamic programming**

Trace from every retained ridge seed to its paired valley anchor along decreasing `distance_to_valley_m`. Solve the three-state sequence with allowed transitions `upper→upper|middle|lower`, `middle→middle|lower`, `lower→lower`; the per-cell cost is negative log membership. Basin cells are fixed and excluded. Apply a correction only when total sequence cost improves and affected cells have confidence below 0.35. Record corrected-cell count in the report.

- [x] **Step 6: Run and commit**

Run: `cd rust; cargo test --release --test classification`

Expected: composition, barrier, code-domain, and ≥95% monotonic trace tests pass.

```powershell
git add rust/topo_core/src/postprocess.rs rust/topo_core/src/lib.rs rust/topo_core/tests/classification.rs
git commit -m "feat: compose and clean terrain positions within slope units"
```

---

## Task 12: Replace the Production Pipeline and Write Formal Outputs

**Files:**

- Rewrite: `rust/topo_core/src/pipeline.rs`
- Modify: `rust/topo_core/src/geotiff.rs`
- Modify: `rust/topo_core/examples/run_sample.rs`
- Delete or migrate obsolete calls in: `rust/topo_core/examples/*.rs`
- Modify: `rust/topo_core/tests/classification.rs`

- [ ] **Step 1: Add a failing array-level end-to-end test**

Complete the Task 1 test by specifying this API:

```rust
pub struct ArrayOutputs {
    pub terrain: Vec<u8>,
    pub geomorph_subclass: Vec<u8>,
    pub confidence: Vec<f32>,
    pub diagnostics: DiagnosticLayers,
    pub stats: Vec<(u8, f64)>,
    pub report: String,
}

pub fn run_arrays_for_test(
    dem: &[f32], valid: &[bool], shape: RasterShape, params: &Params,
) -> Result<ArrayOutputs>;
```

Assert `nested_ridges` produces all three slope positions and no class dominates more than 90% of valid cells.

- [ ] **Step 2: Run and confirm the old pipeline fails the contract**

Run: `cd rust; cargo test --release --test classification final_codes_obey_contract`

Expected: compile failure or distribution assertion failure before the rewrite.

- [ ] **Step 3: Reduce `pipeline.rs` to orchestration**

Keep public `run(params, progress, cancelled) -> Result<Outputs>`. Its stage order is fixed:

1. validate and prepare input;
2. create dual surfaces;
3. build hydrology and nested streams;
4. build multiscale terrain pyramid;
5. calculate scale-indexed geomorphons;
6. extract ridges and slope units;
7. calculate constrained geometry, HAND, and relative position;
8. calculate morphology evidence and fuzzy slope positions;
9. detect/reconstruct basins;
10. compose and constrained-clean outputs;
11. upsample coarse products with nearest neighbour for categories and bilinear interpolation for continuous layers, masked by the full-resolution valid mask;
12. write formal outputs, diagnostics, statistics, and report.

At every stage invoke cancellation checks before allocating another global buffer. Drop stage-local vectors once their last consumer finishes.

- [ ] **Step 4: Define diagnostic output exactly**

```rust
pub struct DiagnosticLayers {
    pub hydro_conditioning_depth: Vec<f32>,
    pub stream_level: Vec<u8>,
    pub ridge_mask: Vec<u8>,
    pub slope_unit: Vec<u32>,
    pub adaptive_scale_m: Vec<f32>,
    pub hand_m: Vec<f32>,
    pub relative_position: Vec<f32>,
    pub slope_position_raw: Vec<u8>,
    pub basin_candidate: Vec<u8>,
    pub basin_core: Vec<u8>,
    pub basin_mask: Vec<u8>,
}
```

Extend `geotiff.rs` with `write_u32` for `slope_unit.tif`. When `write_diagnostics=false`, retain diagnostics only until formal output creation and skip the directory entirely.

- [ ] **Step 5: Write filenames and report schema**

Always write:

- `terrain_position.tif`
- `geomorph_subclass.tif`
- `terrain_confidence.tif`
- `class_report.txt`

When enabled, write the eleven fields above under `diagnostics/` with matching snake-case filenames. `class_report.txt` contains input metadata, resolved parameter values, usable scales, class pixel count/area/percentage, basin object table, invalid count, low-confidence percentage, monotonic correction count, elapsed time per stage, and warnings. It explicitly states “code 2 reserved; emitted count must be 0”.

- [ ] **Step 6: Remove legacy production behavior**

Delete `SeedMode`, fixed `slope_search_m` production use, global histogram-matching code, Euclidean nearest-river logic, elevation-versus-focal-mean basin percentile, and unused basin bridge/merge fields. Migrate diagnostic examples to the new module APIs or delete an example when it exists only to compare retired schemes. This is a direct replacement; do not preserve a hidden old-mode flag.

- [ ] **Step 7: Run end-to-end tests and sample smoke**

Run:

```powershell
cd rust
cargo test --release --test classification
cargo run -p topo_core --example run_sample --release
```

Expected: all invariant tests pass; sample run has code 2 count zero, upper/middle/lower all present, and no non-basin class above 90% of valid area.

- [ ] **Step 8: Commit**

```powershell
git add rust/topo_core/src/pipeline.rs rust/topo_core/src/geotiff.rs rust/topo_core/examples rust/topo_core/tests/classification.rs
git commit -m "feat: integrate adaptive dem-only terrain pipeline"
```

---

## Task 13: Prove Every Exposed Parameter Has an Intended Effect

**Files:**

- Create: `rust/topo_core/tests/parameter_effects.rs`
- Modify: `rust/topo_core/src/pipeline.rs`

- [ ] **Step 1: Add one sensitivity test per visible setting**

Use the same deterministic mixed terrain. Required assertions:

- `Fast→Detailed` increases coarse analytical resolution and/or usable scale count while preserving dimensions.
- `Strict→Loose` monotonically increases or preserves basin candidate area.
- increasing `basin_min_area_m2` monotonically decreases or preserves accepted basin area.
- increasing `postprocess_strength` changes only low-confidence small patches and never barriers.
- `write_diagnostics` changes file emission only, not terrain codes.
- changing `hydro_z_limit_m`, growth threshold, low-relief threshold, hill elevation threshold, or stream thresholds changes its named intermediate layer or expected final zoning.

- [ ] **Step 2: Run and confirm any dead parameter fails**

Run: `cd rust; cargo test --release --test parameter_effects`

Expected before adjustment: at least one assertion exposes any still-unused field.

- [ ] **Step 3: Centralize preset resolution**

Add `Params::resolved() -> Result<ResolvedParams>`. Map presets exactly:

| Preset | coarse resolution | scale family | tile edge |
|---|---:|---|---:|
| Fast | max(DEM resolution, 50 m) | 250–4000 m | 1024 px |
| Standard | max(DEM resolution, 25 m) | 125–4000 m | 1024 px |
| Detailed | max(DEM resolution, 10 m) | 125–min(8000,data max) m | 768 px |

Advanced non-default values override only their matching resolved field. Validate ordered stream thresholds, `0.05≤growth_threshold≤0.40`, `0≤postprocess_strength≤2`, and positive physical distances/areas.

- [ ] **Step 4: Run and commit**

Run: `cd rust; cargo test --release --test parameter_effects`

Expected: all parameter-effect tests pass and no field is read only for reporting.

```powershell
git add rust/topo_core/src/pipeline.rs rust/topo_core/tests/parameter_effects.rs
git commit -m "test: verify adaptive terrain parameter effects"
```

---

## Task 14: Build Python Reference, Validation, and QA Atlas Scripts

**Files:**

- Create: `python/build_synthetic_dem.py`
- Create: `python/adaptive_terrain_reference.py`
- Create: `python/validate_terrain_result.py`
- Create: `python/render_qa_atlas.py`
- Modify: `python/README.md`

- [ ] **Step 1: Create synthetic raster generator**

Put these constants at file top:

```python
OUTPUT_DIR = Path(__file__).resolve().parents[1] / "data" / "synthetic"
CRS_WKT_SOURCE = Path(__file__).resolve().parents[1] / "data" / "dem.tif"
RESOLUTIONS_M = (5.0, 10.0, 25.0)
NODATA = -9999.0
```

Write broad basin, V-valley, nested ridge, conical hill, closed pit, vertical-shift, mirror, rotation, and noise variants. Write each expected zone as a byte GeoTIFF beside the DEM. Delete and recreate only `data/synthetic/generated/`, after resolving and asserting it is inside `data/synthetic`; never delete `data/synthetic` itself.

- [ ] **Step 2: Implement the readable reference calculations**

Implement small-raster NumPy/SciPy versions of robust relief, scale selection, flow-connected HAND, `qd/qz/q`, and fuzzy memberships. It reads only synthetic scenes and writes `.npz` arrays. This is an independent oracle for formulas, not a second production classifier.

- [ ] **Step 3: Implement acceptance metrics**

Use fixed top-of-file paths. Calculate:

- valid-pixel coverage and code-domain violations;
- code 2 count;
- per-class count/area/percentage;
- non-boundary agreement after 10 m and 25 m results are resampled to 5 m;
- agreement after vertical shift, mirror, and rotation inversion;
- agreement after 0.25 m and 0.5 m deterministic noise;
- class-area changes for ±20% parameter runs;
- basin-area changes for ±20% parameter runs;
- 500 ridge-to-valley trace monotonicity rate.

Exit code is nonzero unless all hard gates pass: all valid pixels classified; NoData preserved; code 2 zero; 10 m non-boundary agreement ≥85%; 0.5 m noise non-boundary agreement ≥90%; parameter class-area change ≤15%; basin-area change ≤20%; trace monotonicity ≥95%.

- [ ] **Step 4: Render a fixed QA atlas**

Create a 4×4 PNG containing DEM/hillshade, conditioning depth, stream level, ridges, slope units, adaptive scale, HAND, relative position, raw slope position, basin candidate/core/final, terrain class, confidence, class histogram, and warning text. Use discrete legends for categorical rasters and identical map extent for every panel.

- [ ] **Step 5: Run the scripts**

Run from repository root:

```powershell
D:/worker_code/.venvgis/Scripts/python.exe python/build_synthetic_dem.py
D:/worker_code/.venvgis/Scripts/python.exe python/adaptive_terrain_reference.py
D:/worker_code/.venvgis/Scripts/python.exe python/validate_terrain_result.py
D:/worker_code/.venvgis/Scripts/python.exe python/render_qa_atlas.py
```

Expected: all four exit 0; validation prints every hard metric and `PASS`; atlas PNG is created in the configured validation output directory.

- [ ] **Step 6: Document and commit**

Document the constants users edit, exact files read/written, code semantics, thresholds, and the distinction between reference and production results.

```powershell
git add python/build_synthetic_dem.py python/adaptive_terrain_reference.py python/validate_terrain_result.py python/render_qa_atlas.py python/README.md
git commit -m "test: add terrain reference and validation toolkit"
```

---

## Task 15: Simplify the Desktop UI and Surface Diagnostics

**Files:**

- Modify: `rust/topo_app/src/main.rs`
- Modify: `README.md`
- Modify: `docs/RELEASE_NOTES-v0.0.1.md`

- [ ] **Step 1: Add UI-to-core mapping tests**

Move pure mapping into testable functions. Assert every UI selection maps to the corresponding new core field and that retired fixed-window/seed controls are absent from serialized application state.

- [ ] **Step 2: Run the focused app test and observe failure**

Run: `cd rust; cargo test -p topo_app --release`

Expected before UI migration: struct field mismatches with the rewritten core `Params`.

- [ ] **Step 3: Replace basic controls**

The always-visible form contains only DEM path, output directory, precision (`快速/标准/精细`), basin tendency (`严格/标准/宽松`), minimum basin area in mu and m², postprocess strength (`弱/标准/强` mapped to `0.6/1.0/1.5`), and diagnostics checkbox. The advanced collapsible panel exposes the six `AdvancedParams` fields with units and valid ranges.

- [ ] **Step 4: Update progress and preview**

Show the twelve pipeline stages from Task 12. After success preview `terrain_position.tif`; add a layer selector for `terrain_confidence.tif`, `adaptive_scale_m.tif`, `hand_m.tif`, `relative_position.tif`, and `basin_mask.tif` when present. The legend labels code 1 as `山间/宽谷盆地` and omits code 2.

- [ ] **Step 5: Update user documentation**

Document DEM-only input, projected-metre requirement, seven effective classes, automatic hydrology extraction, diagnostic meanings, default strategy, limitations without ground truth, and the acceptance report. Remove claims about fixed slope windows, target proportions, external river extraction, and separate wide-valley code 2.

- [ ] **Step 6: Run tests and commit**

Run: `cd rust; cargo test -p topo_app --release`

Expected: all app mapping tests pass.

```powershell
git add rust/topo_app/src/main.rs README.md docs/RELEASE_NOTES-v0.0.1.md
git commit -m "feat: expose adaptive terrain controls and diagnostics"
```

---

## Task 16: Validate the Real DEM, Complete Quality Gates, Package, and Publish

**Files:**

- Modify only if evidence requires: files owned by Tasks 2–15
- Update with measured results: `docs/RELEASE_NOTES-v0.0.1.md`
- Create during packaging, then remove: `TerraPos-v0.0.1-win64/`
- Replace: `dist/TerraPos-v0.0.1-win64.zip`

- [ ] **Step 1: Run the full classifier on `data/dem.tif`**

Use Standard precision and default parameters. Write to a dedicated validation directory outside committed source files. Record elapsed time, peak memory, class distribution, low-confidence area, basin objects, and warnings. Confirm input CRS is CGCS2000 projected metres and note that all elevations are above 500 m, so hill classes 3–5 are legitimately absent under the approved business rule.

- [ ] **Step 2: Inspect the QA atlas and numeric diagnostics**

Reject the run if streams climb ridges, ridges occupy valley centres, slope units cross primary barriers, HAND is discontinuous along one flow path, basin masks are only eroded cores, or broad rectangular seams follow processing tiles. Numerical acceptance is necessary but not sufficient; record visual findings in release notes without claiming external accuracy.

- [ ] **Step 3: Run transformation and robustness validation**

Run the four Python scripts from Task 14 and preserve their metric summary in the release notes. If a hard gate fails, fix the owning module, add a regression test, rerun that module’s focused test, then repeat this step. Do not relax a threshold solely to make the current DEM pass.

- [ ] **Step 4: Run the mandatory Rust quality gate**

```powershell
cd rust
cargo test --release
cargo clippy --release
cargo build --release
cargo run -p topo_core --example run_sample --release
```

Expected: every test passes, clippy emits zero warnings, release build succeeds, and sample smoke assertions pass with a plausible non-degenerate distribution.

- [ ] **Step 5: Rebuild the portable package**

From repository root, resolve each source and destination path first and assert the temporary directory is exactly `E:\zcode_worker\Topographic\TerraPos-v0.0.1-win64`. Then create:

```text
TerraPos-v0.0.1-win64/
├── TerraPos.exe
├── README.txt
└── RELEASE_NOTES.txt
```

Copy `rust/target/release/topo_app.exe` to `TerraPos.exe`; create `README.txt` from the final product instructions without mentioning or including a `sample/` directory; copy `docs/RELEASE_NOTES-v0.0.1.md` to `RELEASE_NOTES.txt`. Package with:

```powershell
powershell -NoProfile -Command "Compress-Archive -Path 'TerraPos-v0.0.1-win64' -DestinationPath 'dist\TerraPos-v0.0.1-win64.zip' -Force"
```

Inspect the archive, verify it contains exactly the three files, and verify the executable byte size and SHA-256 differ from the previous archive when code changed. Remove only the verified temporary package directory with PowerShell `Remove-Item -LiteralPath 'E:\zcode_worker\Topographic\TerraPos-v0.0.1-win64' -Recurse -Force`.

- [ ] **Step 6: Review the final diff before publishing**

Run:

```powershell
git status --short
git diff --check
git diff --stat
rg -n "兼容|SeedMode|target.*percent|宽谷盆地.*2" rust python README.md docs
```

Expected: no whitespace errors; no accidental compatibility branch, dead legacy strategy, target-ratio logic, or documentation that says code 2 is emitted. Review every remaining search match in context and remove only matches that violate the approved design.

- [ ] **Step 7: Commit and push the complete verified state**

```powershell
git add -A
git commit -m "release: validate adaptive terrain position v0.0.1"
git push origin main
```

If earlier task commits already contain every tracked change, do not create an empty commit; verify `git status --short` is clean and push the existing commits.

- [ ] **Step 8: Replace the existing GitHub Release asset**

```powershell
gh release upload TerraPos-v0.0.1 dist/TerraPos-v0.0.1-win64.zip --clobber
gh release edit TerraPos-v0.0.1 --draft=false
gh release view TerraPos-v0.0.1
git log origin/main -1 --oneline
```

Expected: release tag remains `TerraPos-v0.0.1`, `draft=false`, asset timestamp/size match the new archive, and `origin/main` points to the verified release commit.

---

## Final Acceptance Checklist

- [ ] Only `dem.tif` is required at runtime.
- [ ] The DEM is rejected unless it is Float32, projected, metre-based, square-pixel, and georeferenced.
- [ ] Boundary NoData is preserved through every output.
- [ ] Streams are nested physical-area thresholds derived from the DEM.
- [ ] HAND follows D8 drainage to the first connected stream.
- [ ] Ridges and streams jointly bound slope units.
- [ ] Characteristic scale changes with terrain width and remains stable across 5/10/25 m sampling.
- [ ] Geomorphons expose ten distinct standard forms and consume metre distances exactly once.
- [ ] Upper/middle/lower memberships follow relative valley-to-ridge geometry rather than global target ratios.
- [ ] Basin acceptance uses width, area, internal relief, surrounding rise, hydrology, and three-scale persistence.
- [ ] Basin output reconstructs the full accepted low-flat candidate rather than retaining only its core.
- [ ] Code 2 is reserved and never emitted; code 1 merges mountain/intermontane and broad-valley basin semantics.
- [ ] All valid cells receive one effective class and invalid cells receive code 0.
- [ ] At least 95% of sampled ridge-to-valley traces are monotonic after constrained correction.
- [ ] 10 m non-boundary agreement is at least 85%; 0.5 m-noise agreement is at least 90%.
- [ ] ±20% parameter perturbations stay within the approved class-area and basin-area sensitivity limits.
- [ ] `data/dem.tif` completes within 6 GiB peak memory and has no tile seams.
- [ ] Rust tests, clippy, build, sample smoke, Python validation, package inspection, Git push, and release verification all pass.

## Execution Notes for the Implementing Agent

- Execute tasks in numerical order; Tasks 3–5 may be developed independently only after Task 2 fixes the shared input types.
- Use focused tests after each red/green step. Run the full suite at Task 12 and Task 16.
- If a design decision conflicts with the approved specification, stop and ask Manba rather than inventing a third behavior.
- Keep commits task-scoped. Do not amend or discard unrelated user changes in a dirty worktree.
- A failed quality gate blocks packaging and publication. Report the measured failure and continue fixing the owning task until the gate passes.
