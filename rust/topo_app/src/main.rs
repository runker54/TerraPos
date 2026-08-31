//! TerraPos 地形部位划分工具 — 现代深色 GIS 工作站 UI
//!
//! 设计语言: 分层深色背景 / 青蓝主色 / 卡片化参数(彩色竖条+徽章) /
//! 状态胶囊 / 浮动图例 / 渐变进度条 / 统计卡带 / 空状态引导

#![windows_subsystem = "windows"]
use eframe::egui;
use std::sync::mpsc::{channel, Receiver, Sender};
use std::sync::Arc;
use std::thread::JoinHandle;

// ============================== 设计令牌 ==============================
mod theme {
    use egui::Color32;
    pub const BG_WIN: Color32 = Color32::from_rgb(243, 245, 248);
    pub const BG_PANEL: Color32 = Color32::from_rgb(255, 255, 255);
    pub const BG_CARD: Color32 = Color32::from_rgb(255, 255, 255);
    pub const BG_INPUT: Color32 = Color32::from_rgb(255, 255, 255);
    pub const BG_HOVER: Color32 = Color32::from_rgb(239, 246, 255);
    pub const STROKE: Color32 = Color32::from_rgb(222, 228, 236);
    pub const ACCENT: Color32 = Color32::from_rgb(37, 99, 235);
    pub const ACCENT_DIM: Color32 = Color32::from_rgb(29, 78, 216);
    pub const TEXT: Color32 = Color32::from_rgb(31, 41, 55);
    pub const TEXT_SUB: Color32 = Color32::from_rgb(90, 103, 120);
    pub const TEXT_DIM: Color32 = Color32::from_rgb(150, 160, 175);
    pub const WARN: Color32 = Color32::from_rgb(180, 83, 9);
    pub const OK: Color32 = Color32::from_rgb(22, 130, 93);
    pub const ERR: Color32 = Color32::from_rgb(200, 30, 40);

    pub const SEC: [Color32; 6] = [
        Color32::from_rgb(8, 145, 178),
        Color32::from_rgb(37, 99, 235),
        Color32::from_rgb(22, 163, 74),
        Color32::from_rgb(234, 88, 12),
        Color32::from_rgb(124, 58, 237),
        Color32::from_rgb(100, 116, 139),
    ];

    pub const CLASS_COLORS: [(u8, &str, Color32); 8] = [
        (1, "山间盆地", Color32::from_rgb(51, 178, 229)),
        (3, "丘陵上部", Color32::from_rgb(250, 217, 89)),
        (4, "丘陵中部", Color32::from_rgb(217, 237, 166)),
        (5, "丘陵下部", Color32::from_rgb(153, 199, 102)),
        (6, "山地坡上", Color32::from_rgb(250, 165, 60)),
        (7, "山地坡中", Color32::from_rgb(222, 100, 50)),
        (8, "山地坡下", Color32::from_rgb(107, 68, 35)),
        (2, "宽谷盆地", Color32::from_rgb(102, 217, 242)),
    ];
}

// ============================== 参数模型 ==============================

#[derive(Clone)]
struct Params {
    dem_path: String,
    out_dir: String,
    coarse_res: f64,
    basin_river_acc_km2: f64,
    basin_buffer_m: f64,
    basin_elev_diff_m: f64,
    basin_slope_th: f64,
    basin_relief_m: f64,
    basin_min_area_mu: f64,
    basin_inner_relief_m: f64,
    basin_bridge_m: f64,
    basin_merge_m: f64,
    basin_merge_max_mu: f64,
    basin_smooth_m: f64,
    slope_tpi_focus_m: f64,
    slope_flat_deg: f64,
    slope_min_patch_m2: f64,
    hill_z_max: f64,
    relief_subclass_win: f64,
    relief_low_hill: f64,
    mode_filter_iter: usize,
    min_patch_m2: f64,
}

impl Params {
    fn defaults() -> Self {
        let d = topo_core::pipeline::Params::default();
        Params {
            dem_path: String::new(),
            out_dir: String::new(),
            coarse_res: d.coarse_res,
            basin_river_acc_km2: d.basin_river_acc_km2,
            basin_buffer_m: d.basin_buffer_m,
            basin_elev_diff_m: d.basin_elev_diff_m,
            basin_slope_th: d.basin_slope_th,
            basin_relief_m: d.basin_relief_m,
            basin_min_area_mu: d.basin_min_area_m2 / 666.6667,
            basin_inner_relief_m: d.basin_inner_relief_m,
            basin_bridge_m: d.basin_bridge_m,
            basin_merge_m: d.basin_merge_m,
            basin_merge_max_mu: d.basin_merge_max_m2 / 666.6667,
            basin_smooth_m: d.basin_smooth_m,
            slope_tpi_focus_m: d.slope_tpi_focus_m,
            slope_flat_deg: d.slope_flat_deg,
            slope_min_patch_m2: d.slope_min_patch_m2,
            hill_z_max: d.hill_z_max,
            relief_subclass_win: d.relief_subclass_win,
            relief_low_hill: d.relief_low_hill,
            mode_filter_iter: d.mode_filter_iter,
            min_patch_m2: d.min_patch_m2,
        }
    }
    fn to_core(&self) -> topo_core::pipeline::Params {
        topo_core::pipeline::Params {
            dem_path: self.dem_path.clone(),
            out_dir: self.out_dir.clone(),
            coarse_res: self.coarse_res,
            basin_river_acc_km2: self.basin_river_acc_km2,
            basin_buffer_m: self.basin_buffer_m,
            basin_elev_diff_m: self.basin_elev_diff_m,
            basin_slope_th: self.basin_slope_th,
            basin_relief_m: self.basin_relief_m,
            basin_min_area_m2: self.basin_min_area_mu * 666.6667,
            basin_inner_relief_m: self.basin_inner_relief_m,
            basin_bridge_m: self.basin_bridge_m,
            basin_merge_m: self.basin_merge_m,
            basin_merge_max_m2: self.basin_merge_max_mu * 666.6667,
            basin_smooth_m: self.basin_smooth_m,
            slope_tpi_focus_m: self.slope_tpi_focus_m,
            slope_flat_deg: self.slope_flat_deg,
            slope_min_patch_m2: self.slope_min_patch_m2,
            hill_z_max: self.hill_z_max,
            relief_subclass_win: self.relief_subclass_win,
            relief_low_hill: self.relief_low_hill,
            mode_filter_iter: self.mode_filter_iter,
            min_patch_m2: self.min_patch_m2,
        }
    }
}

enum WorkerMsg {
    Progress(String, f32, String),
    Done(String, Vec<(u8, f64)>),
    Failed(String),
    /// (代号, 高程rgba, 阴影rgba, w, h): DEM 双预览
    DemPreview(u64, Vec<u8>, Vec<u8>, usize, usize),
}

#[derive(Clone, Copy, PartialEq)]
enum Layer {
    Dem,
    DemShade,
    Terrain,
    Subclass,
}

#[derive(Clone, Copy, PartialEq)]
enum Page {
    Workbench,
    ParamDocs,
}

struct StatRow {
    name: &'static str,
    color: egui::Color32,
    area_km2: f64,
}


// ============================== 主题/字体 ==============================
fn apply_theme(ctx: &egui::Context) {
    let mut st = egui::Style::default();
    let v = &mut st.visuals;
    v.dark_mode = true;
    v.panel_fill = theme::BG_PANEL;
    v.extreme_bg_color = theme::BG_WIN;
    v.window_fill = theme::BG_CARD;
    v.widgets.noninteractive.bg_fill = theme::BG_CARD;
    v.widgets.noninteractive.fg_stroke = egui::Stroke::new(1.0, theme::TEXT_SUB);
    v.widgets.noninteractive.bg_stroke = egui::Stroke::new(1.0, theme::STROKE);
    v.widgets.inactive.bg_fill = theme::BG_INPUT;
    v.widgets.inactive.fg_stroke = egui::Stroke::new(1.0, theme::TEXT);
    v.widgets.inactive.bg_stroke = egui::Stroke::new(1.0, theme::STROKE);
    v.widgets.hovered.bg_fill = theme::BG_HOVER;
    v.widgets.hovered.fg_stroke = egui::Stroke::new(1.2, theme::ACCENT);
    v.widgets.hovered.bg_stroke = egui::Stroke::new(1.0, theme::ACCENT);
    v.widgets.active.bg_fill = theme::BG_HOVER;
    v.widgets.active.fg_stroke = egui::Stroke::new(1.0, theme::ACCENT_DIM);
    v.selection.stroke = egui::Stroke::new(1.6, theme::ACCENT);
    v.selection.bg_fill = theme::ACCENT.gamma_multiply(0.14);
    st.spacing.item_spacing = egui::vec2(8.0, 7.0);
    st.spacing.button_padding = egui::vec2(10.0, 5.0);
    st.text_styles = [
        (egui::TextStyle::Heading, egui::FontId::proportional(20.0)),
        (egui::TextStyle::Body, egui::FontId::proportional(13.5)),
        (egui::TextStyle::Small, egui::FontId::proportional(11.0)),
        (egui::TextStyle::Button, egui::FontId::proportional(13.5)),
        (egui::TextStyle::Monospace, egui::FontId::monospace(12.0)),
    ]
    .into();
    ctx.set_style(st);
    setup_chinese_fonts(ctx);
}

fn setup_chinese_fonts(ctx: &egui::Context) {
    let candidates = [
        r"C:\Windows\Fonts\msyh.ttc",
        r"C:\Windows\Fonts\simhei.ttf",
        r"C:\Windows\Fonts\Deng.ttf",
        r"C:\Windows\Fonts\simsun.ttc",
    ];
    let mut fonts = egui::FontDefinitions::default();
    for path in candidates {
        if let Ok(bytes) = std::fs::read(path) {
            fonts
                .font_data
                .insert("cjk".to_owned(), std::sync::Arc::new(egui::FontData::from_owned(bytes)));
            for fam in [egui::FontFamily::Proportional, egui::FontFamily::Monospace] {
                fonts.families.entry(fam).or_default().insert(0, "cjk".to_owned());
            }
            ctx.set_fonts(fonts);
            eprintln!("已加载中文字体: {path}");
            return;
        }
    }
    eprintln!("警告: 未找到系统中文字体, 界面中文将无法显示");
}

// ============================== 通用组件 ==============================
fn card<R>(
    ui: &mut egui::Ui,
    accent: egui::Color32,
    badge: &str,
    title: &str,
    body: impl FnOnce(&mut egui::Ui) -> R,
) -> R {
    egui::Frame::new()
        .fill(theme::BG_CARD)
        .stroke(egui::Stroke::new(1.0, theme::STROKE))
        .corner_radius(8.0)
        .inner_margin(egui::Margin::same(10))
        .outer_margin(egui::Margin { left: 4, right: 0, top: 0, bottom: 0 })
        .show(ui, |ui| {
            ui.horizontal(|ui| {
                let (r, _) = ui.allocate_exact_size(egui::vec2(3.5, 16.0), egui::Sense::hover());
                ui.painter().rect_filled(r, 2.0, accent);
                ui.label(egui::RichText::new(badge).monospace().color(accent).strong());
                ui.label(egui::RichText::new(title).strong());
            });
            ui.add_space(2.0);
            body(ui)
        })
    .inner
}

fn row<R>(ui: &mut egui::Ui, label: &str, add: impl FnOnce(&mut egui::Ui) -> R) -> R {
    ui.horizontal(|ui| {
        ui.add_sized(
            [148.0, 18.0],
            egui::Label::new(egui::RichText::new(label).color(theme::TEXT_SUB)),
        );
        add(ui)
    })
    .inner
}

fn num(ui: &mut egui::Ui, label: &str, v: &mut f64, speed: f64) {
    row(ui, label, |ui| {
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            ui.set_min_width(ui.available_width());
            ui.add_sized([ui.available_width(), 20.0], egui::DragValue::new(v).speed(speed).max_decimals(2))
        });
    });
}

fn stepper(ui: &mut egui::Ui, label: &str, v: &mut usize, step: usize, min: usize) {
    row(ui, label, |ui| {
        if ui.add(egui::Button::new("−").small()).clicked() && *v > min {
            *v -= step;
        }
        ui.monospace(format!("{v}"));
        if ui.add(egui::Button::new("+").small()).clicked() {
            *v += step;
        }
    });
}

fn status_pill(ui: &mut egui::Ui, running: bool, stage: &str, pct: f32) {
    let (color, text) = if running {
        (theme::WARN, format!("⟳ {stage} {pct:.0}%"))
    } else if !stage.is_empty() {
        (theme::OK, "✓ 就绪".to_string())
    } else {
        (theme::TEXT_DIM, "待机".to_string())
    };
    egui::Frame::new()
        .fill(color.gamma_multiply(0.16))
        .stroke(egui::Stroke::new(1.0, color.gamma_multiply(0.55)))
        .corner_radius(10.0)
        .inner_margin(egui::Margin::symmetric(10, 4))
        .show(ui, |ui| {
            ui.label(egui::RichText::new(text).color(color).small().strong());
        });
}

fn primary_button(ui: &mut egui::Ui, text: &str, enabled: bool) -> bool {
    let size = egui::vec2(ui.available_width(), 34.0);
    ui.add_enabled(
        enabled,
        egui::Button::new(egui::RichText::new(text).strong().size(14.5))
            .fill(if enabled { theme::ACCENT_DIM } else { theme::BG_INPUT })
            .corner_radius(8.0)
            .min_size(size),
    )
    .clicked()
}

fn ghost_button(ui: &mut egui::Ui, text: &str) -> bool {
    ui.add_sized(
        [ui.available_width(), 26.0],
        egui::Button::new(egui::RichText::new(text).color(theme::TEXT_SUB))
            .fill(egui::Color32::TRANSPARENT)
            .stroke(egui::Stroke::new(1.0, theme::STROKE))
            .corner_radius(6.0),
    )
    .clicked()
}

fn chip(ui: &mut egui::Ui, text: &str, selected: bool) -> bool {
    let (fill, stroke, txt) = if selected {
        (
            theme::ACCENT.gamma_multiply(0.25),
            egui::Stroke::new(1.2, theme::ACCENT),
            egui::RichText::new(text).color(theme::ACCENT).strong(),
        )
    } else {
        (
            theme::BG_CARD,
            egui::Stroke::new(1.0, theme::STROKE),
            egui::RichText::new(text).color(theme::TEXT_SUB),
        )
    };
    ui.add(egui::Button::new(txt).fill(fill).stroke(stroke).corner_radius(12.0))
        .clicked()
}

// ============================== App ==============================
struct App {
    params: Params,
    page: Page,
    running: bool,
    request_start: bool,
    cancel: Arc<std::sync::atomic::AtomicBool>,
    worker: Option<JoinHandle<()>>,
    rx: Option<Receiver<WorkerMsg>>,
    log: Vec<(String, bool)>,
    last_stage: String,
    last_pct: f32,
    result_summary: Option<String>,
    stats: Vec<StatRow>,
    tex_terrain: Option<egui::TextureHandle>,
    tex_sub: Option<egui::TextureHandle>,
    tex_dem: Option<egui::TextureHandle>,
    tex_dem_shade: Option<egui::TextureHandle>,
    dem_rx: Option<Receiver<WorkerMsg>>,
    dem_gen: u64,
    dem_loading: bool,
    layer: Layer,
    view_scale: f32,
    view_ox: f32,
    view_oy: f32,
    data_wh: (f32, f32),
}

impl App {
    fn new(cc: &eframe::CreationContext<'_>) -> Self {
        apply_theme(&cc.egui_ctx);
        let mut app = App {
            params: Params::defaults(),
            page: Page::Workbench,
            running: false,
            request_start: false,
            cancel: Arc::new(std::sync::atomic::AtomicBool::new(false)),
            worker: None,
            rx: None,
            log: Vec::new(),
            last_stage: String::new(),
            last_pct: 0.0,
            result_summary: None,
            stats: Vec::new(),
            tex_terrain: None,
            tex_sub: None,
            tex_dem: None,
            tex_dem_shade: None,
            dem_rx: None,
            dem_gen: 0,
            dem_loading: false,
            layer: Layer::Terrain,
            view_scale: 0.12,
            view_ox: 0.0,
            view_oy: 0.0,
            data_wh: (14831.0, 9058.0),
        };
        if std::env::var("TOPOSUITE_PAGE").as_deref() == Ok("docs") {
            app.page = Page::ParamDocs;
        }
        if let Ok(dp) = std::env::var("TOPOSUITE_DEMPATH") {
            app.params.dem_path = dp;
            if let Ok(od) = std::env::var("TOPOSUITE_OUT") {
                app.params.out_dir = od;
            }
        }
        if let Ok(demo) = std::env::var("TOPOSUITE_DEMO") {
            app.params.dem_path = demo.clone();
            if let Ok(outd) = std::env::var("TOPOSUITE_OUT") {
                app.params.out_dir = outd;
            } else {
                app.params.out_dir = "demo_out".into();
            }
            app.request_start = true;
        }
        app
    }

    /// 后台加载 DEM 并生成山体阴影预览纹理
    fn spawn_dem_preview(&mut self, path: String, ctx: &egui::Context) {
        self.dem_gen += 1;
        self.dem_loading = true;
        self.tex_dem = None;
        let gen = self.dem_gen;
        let (tx, rx): (Sender<WorkerMsg>, Receiver<WorkerMsg>) = channel();
        self.dem_rx = Some(rx);
        let ctx2 = ctx.clone();
        std::thread::spawn(move || {
            let res = (|| -> Result<(Vec<u8>, Vec<u8>, usize, usize), String> {
                let (dem, meta) =
                    topo_core::geotiff::read_f32(&path).map_err(|e| e.to_string())?;
                let (w, h) = (meta.width as usize, meta.height as usize);
                let res_m = meta.resolution();
                // 预览目标 ~1600px 宽; 降采样后立即释放原始大数组
                let (small, sw, sh) = if w > 1600 {
                    let dst = res_m * (w as f64 / 1600.0);
                    let (s, sw2, sh2) =
                        topo_core::terrain::box_downsample(&dem, w, h, res_m, dst);
                    (s, sw2, sh2)
                } else {
                    (dem, w, h)
                };
                let hypso = topo_core::terrain::hypsometric_rgba(&small, sw, sh);
                let shade = topo_core::terrain::hillshade_rgba(&small, sw, sh, res_m);
                Ok((hypso, shade, sw, sh))
            })();
            match res {
                Ok((hypso, shade, w, h)) => {
                    let _ = tx.send(WorkerMsg::DemPreview(gen, hypso, shade, w, h));
                }
                Err(e) => {
                    let _ = tx.send(WorkerMsg::Failed(e));
                }
            }
            ctx2.request_repaint();
        });
    }

    fn poll_dem(&mut self, ctx: &egui::Context) {
        if let Some(rx) = &self.dem_rx {
            while let Ok(msg) = rx.try_recv() {
                match msg {
                    WorkerMsg::DemPreview(gen, hypso, shade, w, h) if gen == self.dem_gen => {
                        self.tex_dem = Some(ctx.load_texture(
                            "tex_dem",
                            egui::ColorImage::from_rgba_unmultiplied([w, h], &hypso),
                            egui::TextureOptions::LINEAR,
                        ));
                        self.tex_dem_shade = Some(ctx.load_texture(
                            "tex_dem_shade",
                            egui::ColorImage::from_rgba_unmultiplied([w, h], &shade),
                            egui::TextureOptions::LINEAR,
                        ));
                        self.data_wh = (w as f32, h as f32);
                        self.dem_loading = false;
                        self.layer = Layer::Dem;
                        self.view_scale = ((1180.0 / w as f32).min(600.0 / h as f32)).max(0.02);
                        self.view_ox = 0.0;
                        self.view_oy = 0.0;
                        self.log.push(("✓ DEM 预览已加载".into(), false));
                    }
                    WorkerMsg::Failed(e) if self.dem_loading => {
                        self.dem_loading = false;
                        self.log.push((format!("✗ DEM 预览失败: {e}"), true));
                    }
                    _ => {}
                }
            }
        }
    }

    fn start_run(&mut self, ctx: &egui::Context) {
        let params = self.params.to_core();
        self.running = true;
        self.cancel
            .store(false, std::sync::atomic::Ordering::Relaxed);
        self.log.clear();
        self.result_summary = None;
        self.stats.clear();
        let (tx, rx): (Sender<WorkerMsg>, Receiver<WorkerMsg>) = channel();
        self.rx = Some(rx);
        let cancel = Arc::clone(&self.cancel);
        let ctx2 = ctx.clone();
        self.worker = Some(std::thread::spawn(move || {
            let res = topo_core::pipeline::run(
                &params,
                &|p| {
                    let _ = tx.send(WorkerMsg::Progress(p.stage.clone(), p.pct, p.msg.clone()));
                    !cancel.load(std::sync::atomic::Ordering::Relaxed)
                },
                &cancel,
            );
            match res {
                Ok(out) => {
                    let _ = tx.send(WorkerMsg::Done(out.report, out.stats));
                }
                Err(topo_core::error::CoreError::Cancelled) => {
                    let _ = tx.send(WorkerMsg::Failed("已取消".into()));
                }
                Err(e) => {
                    let _ = tx.send(WorkerMsg::Failed(format!("{e}")));
                }
            }
            // DemPreview 走独立通道(dem_rx), 该通道不会收到
            ctx2.request_repaint();
        }));
    }

    fn poll_worker(&mut self, ctx: &egui::Context) {
        if let Some(rx) = &self.rx {
            while let Ok(msg) = rx.try_recv() {
                match msg {
                    WorkerMsg::Progress(stage, pct, m) => {
                        self.last_stage = stage;
                        self.last_pct = pct;
                        self.log.push((format!("[{pct:5.1}%] {m}"), false));
                        if self.log.len() > 500 {
                            self.log.remove(0);
                        }
                    }
                    WorkerMsg::DemPreview(..) => {}
                    WorkerMsg::Done(report, stats) => {
                        self.result_summary = Some(report);
                        self.log.push(("✓ 划分完成, 成果已写入输出目录".into(), false));
                        let base = std::path::Path::new(&self.params.out_dir);
                        if let Ok((data, w, h)) = topo_core::geotiff::read_u8_preview(
                            base.join("terrain_position.tif"),
                            2400,
                        ) {
                            let rgba = colorize(&data);
                            self.tex_terrain = Some(ctx.load_texture(
                                "tex_terrain",
                                egui::ColorImage::from_rgba_unmultiplied([w, h], &rgba),
                                egui::TextureOptions::LINEAR,
                            ));
                            self.data_wh = (w as f32, h as f32);
                            let names = [
                                (1u8, "山间盆地", theme::CLASS_COLORS[0].2),
                                (6, "山地坡上", theme::CLASS_COLORS[4].2),
                                (7, "山地坡中", theme::CLASS_COLORS[5].2),
                                (8, "山地坡下", theme::CLASS_COLORS[6].2),
                                (3, "丘陵上部", theme::CLASS_COLORS[1].2),
                                (4, "丘陵中部", theme::CLASS_COLORS[2].2),
                                (5, "丘陵下部", theme::CLASS_COLORS[3].2),
                            ];
                            self.stats = stats
                                .iter()
                                .filter_map(|(c, a)| {
                                    names.iter().find(|(k, _, _)| k == c).map(|(_, n, col)| StatRow {
                                        name: n,
                                        color: *col,
                                        area_km2: *a,
                                    })
                                })
                                .collect();
                            self.layer = Layer::Terrain;
                            // 自适应画布(预览2400宽)
                            self.view_scale = ((1180.0 / w as f32).min(600.0 / h as f32)).max(0.02);
                            self.view_ox = 0.0;
                            self.view_oy = 0.0;
                        }
                        if let Ok((data, w, h)) = topo_core::geotiff::read_u8_preview(
                            base.join("geomorph_subclass.tif"),
                            2400,
                        ) {
                            let rgba = colorize_sub(&data);
                            self.tex_sub = Some(ctx.load_texture(
                                "tex_sub",
                                egui::ColorImage::from_rgba_unmultiplied([w, h], &rgba),
                                egui::TextureOptions::LINEAR,
                            ));
                        }
                        self.running = false;
                    }
                    WorkerMsg::Failed(e) => {
                        self.log.push((format!("✗ 失败: {e}"), true));
                        self.running = false;
                    }
                }
            }
        }
    }

    // ---------- 左侧面板: 导航 + 页面 ----------
    fn params_panel(&mut self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            let w = (ui.available_width() - 8.0) / 2.0;
            let nav = |ui: &mut egui::Ui, text: &str, sel: bool| {
                ui.add_sized(
                    [w, 30.0],
                    egui::Button::new(egui::RichText::new(text).strong())
                        .fill(if sel { theme::ACCENT.gamma_multiply(0.16) } else { theme::BG_CARD })
                        .stroke(egui::Stroke::new(
                            1.0,
                            if sel { theme::ACCENT } else { theme::STROKE },
                        ))
                        .corner_radius(8.0),
                )
                .clicked()
            };
            if nav(ui, "⚙  划分", self.page == Page::Workbench) {
                self.page = Page::Workbench;
            }
            if nav(ui, "📖  参数详解", self.page == Page::ParamDocs) {
                self.page = Page::ParamDocs;
            }
        });
        ui.add_space(6.0);
        match self.page {
            Page::Workbench => self.workbench_panel(ui),
            Page::ParamDocs => self.docs_panel(ui),
        }
    }

    // ---------- 参数详解页面 ----------
    fn docs_panel(&mut self, ui: &mut egui::Ui) {
        egui::ScrollArea::vertical().show(ui, |ui| {
            ui.label(
                egui::RichText::new(
                    "以下为全部判断指标的详细说明：每个参数控制划分结果的一个方面。                     先阅读分组说明理解其地貌含义，再按研究区实际情况调整；                     拿不准时保持默认值——默认值已与本区验证成果对齐。",
                )
                .color(theme::TEXT_SUB)
                .small(),
            );
            ui.add_space(8.0);
            for (group, accent, items) in PARAM_DOCS {
                egui::Frame::new()
                    .fill(theme::BG_CARD)
                    .stroke(egui::Stroke::new(1.0, theme::STROKE))
                    .corner_radius(8.0)
                    .inner_margin(egui::Margin::same(10))
                    .outer_margin(egui::Margin { bottom: 8, ..Default::default() })
                    .show(ui, |ui| {
                        ui.horizontal(|ui| {
                            let (r, _) = ui
                                .allocate_exact_size(egui::vec2(3.5, 16.0), egui::Sense::hover());
                            ui.painter().rect_filled(r, 2.0, *accent);
                            ui.label(egui::RichText::new(*group).strong().color(*accent));
                        });
                        ui.add_space(4.0);
                        for item in items.iter() {
                            egui::CollapsingHeader::new(
                                egui::RichText::new(item.name).strong().size(13.0),
                            )
                            .id_salt((*group, item.name))
                            .default_open(true)
                            .show(ui, |ui| {
                                ui.horizontal(|ui| {
                                    ui.label(
                                        egui::RichText::new("默认").small().color(theme::TEXT_DIM),
                                    );
                                    ui.label(
                                        egui::RichText::new(item.default)
                                            .monospace()
                                            .small()
                                            .color(theme::ACCENT_DIM),
                                    );
                                });
                                ui.label(egui::RichText::new(item.what).small());
                                ui.label(
                                    egui::RichText::new(item.effect).small().color(theme::TEXT_SUB),
                                );
                                ui.add_space(2.0);
                            });
                        }
                    });
            }
            ui.add_space(10.0);
            ui.label(
                egui::RichText::new(
                    "提示：坝子最低宽度、坡上/坡下阈值、山体单元窗口是对结果影响最大的三个旋钮；                     其余参数多数区域保持默认即可。",
                )
                .small()
                .color(theme::TEXT_DIM),
            );
        });
    }

    fn workbench_panel(&mut self, ui: &mut egui::Ui) {
        egui::ScrollArea::vertical().show(ui, |ui| {
            // DEM 预览未完成时禁用运行: 避免预览线程与管线并发各持一份全量 DEM
            let can_run =
                !self.params.dem_path.is_empty() && !self.params.out_dir.is_empty() && !self.dem_loading;

            card(ui, theme::SEC[0], "01", "数据导入", |ui| {
                if ui
                    .add_sized(
                        [ui.available_width(), 30.0],
                        egui::Button::new("🗂  选择 DEM GeoTIFF..."),
                    )
                    .clicked()
                {
                    if let Some(p) = rfd::FileDialog::new()
                        .add_filter("GeoTIFF", &["tif", "tiff"])
                        .pick_file()
                    {
                        self.params.dem_path = p.display().to_string();
                        if self.params.out_dir.is_empty() {
                            // 默认输出到程序所在目录(而非栅格所在位置)
                            let base = std::env::current_exe()
                                .ok()
                                .and_then(|e| e.parent().map(|d| d.to_path_buf()))
                                .or_else(|| p.parent().map(|d| d.to_path_buf()))
                                .unwrap_or_else(|| std::path::PathBuf::from("."));
                            self.params.out_dir = base.join("topo_out").display().to_string();
                        }
                        self.spawn_dem_preview(self.params.dem_path.clone(), ui.ctx());
                    }
                }
                ui.label(
                    egui::RichText::new(if self.params.dem_path.is_empty() {
                        "未选择 — 支持 float32 投影坐标系 GeoTIFF".to_string()
                    } else {
                        self.params
                            .dem_path
                            .rsplit(['\\', '/'])
                            .next()
                            .unwrap_or("")
                            .to_string()
                    })
                    .small()
                    .color(if self.params.dem_path.is_empty() {
                        theme::TEXT_DIM
                    } else {
                        theme::OK
                    }),
                );
                ui.add_space(2.0);
                ui.horizontal(|ui| {
                    ui.label(egui::RichText::new("输出:").small().color(theme::TEXT_DIM));
                    ui.label(
                        egui::RichText::new(if self.params.out_dir.is_empty() {
                            "(随 DEM 自动设定)".to_string()
                        } else {
                            self.params.out_dir.clone()
                        })
                        .small()
                        .color(theme::TEXT_DIM),
                    );
                    if ui.small_button("更改").clicked() {
                        if let Some(p) = rfd::FileDialog::new().pick_folder() {
                            self.params.out_dir = p.display().to_string();
                        }
                    }
                });
            });
            ui.add_space(8.0);

            card(ui, theme::SEC[1], "02", "山间盆地判别", |ui| {
                ui.label(
                    egui::RichText::new("河谷低平带(legacy 同参数) + 对象级内部起伏检验")
                        .small()
                        .color(theme::TEXT_DIM),
                );
                num(ui, "河网阈值 (km²)", &mut self.params.basin_river_acc_km2, 0.5);
                num(ui, "河流缓冲 (m)", &mut self.params.basin_buffer_m, 50.0);
                num(ui, "河流高程差上限 (m)", &mut self.params.basin_elev_diff_m, 0.5);
                num(ui, "坡度上限 (°)", &mut self.params.basin_slope_th, 0.5);
                num(ui, "局部起伏上限 (m)", &mut self.params.basin_relief_m, 0.5);
                ui.separator();
                num(ui, "最低面积 (亩, 不足转坡下)", &mut self.params.basin_min_area_mu, 5.0);
                num(ui, "内部起伏上限 (m)", &mut self.params.basin_inner_relief_m, 1.0);
                num(ui, "碎片桥接 (m)", &mut self.params.basin_bridge_m, 10.0);
                num(ui, "碎斑归并 (m, 0=关)", &mut self.params.basin_merge_m, 10.0);
                num(ui, "碎斑面积上限 (亩)", &mut self.params.basin_merge_max_mu, 1.0);
                num(ui, "平滑距离 (m)", &mut self.params.basin_smooth_m, 10.0);
            });
            ui.add_space(8.0);

            card(ui, theme::SEC[2], "03", "坡位判别", |ui| {
                ui.label(
                    egui::RichText::new("TPI 分类: 山谷/坡下/平坡/坡中/坡上/山脊 (脚本 focus=101)")
                        .small()
                        .color(theme::TEXT_SUB),
                );
                num(ui, "TPI 焦点窗 (m)", &mut self.params.slope_tpi_focus_m, 25.0);
                num(ui, "平坡坡度分界 (°)", &mut self.params.slope_flat_deg, 0.5);
                num(ui, "坡位小斑蚕食 (m²)", &mut self.params.slope_min_patch_m2, 1000.0);
            });
            ui.add_space(8.0);

            card(ui, theme::SEC[3], "04", "丘陵 / 山地", |ui| {
                num(ui, "丘陵海拔上限 (m)", &mut self.params.hill_z_max, 50.0);
                num(ui, "亚类起伏度窗口 (m)", &mut self.params.relief_subclass_win, 100.0);
                num(ui, "低丘起伏度上限 (m)", &mut self.params.relief_low_hill, 25.0);
            });
            ui.add_space(8.0);

            card(ui, theme::SEC[4], "05", "后处理", |ui| {
                stepper(ui, "众数滤波轮数", &mut self.params.mode_filter_iter, 1, 0);
                num(ui, "最小图斑 (m²)", &mut self.params.min_patch_m2, 1000.0);
            });
            ui.add_space(8.0);

            egui::Frame::new()
                .fill(theme::BG_CARD)
                .stroke(egui::Stroke::new(1.0, theme::STROKE))
                .corner_radius(8.0)
                .inner_margin(egui::Margin::same(10))
                .outer_margin(egui::Margin { left: 4, right: 0, top: 0, bottom: 0 })
                .show(ui, |ui| {
                    ui.horizontal(|ui| {
                        let (r, _) =
                            ui.allocate_exact_size(egui::vec2(3.5, 16.0), egui::Sense::hover());
                        ui.painter().rect_filled(r, 2.0, theme::SEC[5]);
                        ui.label(egui::RichText::new("06").monospace().color(theme::SEC[5]));
                        ui.label(egui::RichText::new("高级参数").strong());
                    });
                    egui::CollapsingHeader::new(
                        egui::RichText::new("展开 ▾").small().color(theme::TEXT_DIM),
                    )
                    .id_salt("adv")
                    .show_unindented(ui, |ui| {
                        num(ui, "中间层分辨率 (m)", &mut self.params.coarse_res, 5.0);
                    });
                });
            ui.add_space(10.0);

            if primary_button(
                ui,
                if self.running { "计算中…" } else { "▶  运行划分" },
                can_run && !self.running,
            ) {
                self.request_start = true;
            }
            ui.add_space(6.0);
            if ghost_button(ui, "恢复默认参数") {
                let d = Params::defaults();
                self.params.coarse_res = d.coarse_res;

                self.params.basin_river_acc_km2 = d.basin_river_acc_km2;
                self.params.basin_buffer_m = d.basin_buffer_m;
                self.params.basin_elev_diff_m = d.basin_elev_diff_m;
                self.params.basin_slope_th = d.basin_slope_th;
                self.params.basin_relief_m = d.basin_relief_m;
                self.params.basin_inner_relief_m = d.basin_inner_relief_m;
                self.params.basin_bridge_m = d.basin_bridge_m;
                self.params.basin_merge_m = d.basin_merge_m;
                self.params.basin_merge_max_mu = d.basin_merge_max_mu;
                self.params.basin_min_area_mu = d.basin_min_area_mu;
                self.params.basin_smooth_m = d.basin_smooth_m;
                self.params.slope_tpi_focus_m = d.slope_tpi_focus_m;
                self.params.slope_flat_deg = d.slope_flat_deg;
                self.params.slope_min_patch_m2 = d.slope_min_patch_m2;
                self.params.hill_z_max = d.hill_z_max;
                self.params.relief_subclass_win = d.relief_subclass_win;
                self.params.relief_low_hill = d.relief_low_hill;
                self.params.mode_filter_iter = d.mode_filter_iter;
                self.params.min_patch_m2 = d.min_patch_m2;
            }
            ui.add_space(4.0);
        });
    }

    // ---------- 中央预览 ----------
    fn preview_panel(&mut self, ui: &mut egui::Ui) {
        let has_result = self.tex_terrain.is_some();
        ui.add_space(6.0);
        ui.horizontal(|ui| {
            ui.add_space(8.0);
            if self.tex_dem.is_some() || self.dem_loading {
                ui.label(egui::RichText::new("预览").strong());
                ui.separator();
                if chip(ui, "DEM 高程", self.layer == Layer::Dem) {
                    self.layer = Layer::Dem;
                }
                if self.tex_dem_shade.is_some()
                    && chip(ui, "DEM 阴影", self.layer == Layer::DemShade)
                {
                    self.layer = Layer::DemShade;
                }
            }
            if has_result {
                if self.tex_dem.is_some() {
                    ui.separator();
                } else {
                    ui.label(egui::RichText::new("结果预览").strong());
                    ui.separator();
                }
                if chip(ui, "地形部位", self.layer == Layer::Terrain) {
                    self.layer = Layer::Terrain;
                }
                if chip(ui, "地貌亚类", self.layer == Layer::Subclass) {
                    self.layer = Layer::Subclass;
                }
                ui.separator();
                if ui.small_button("−").clicked() {
                    self.view_scale = (self.view_scale / 1.25).max(0.008);
                }
                ui.monospace(format!("{:.1}x", self.view_scale / 0.12));
                if ui.small_button("+").clicked() {
                    self.view_scale = (self.view_scale * 1.25).min(1.0);
                }
                if ui.small_button("复位").clicked() {
                    self.view_scale = 0.12;
                    self.view_ox = 0.0;
                    self.view_oy = 0.0;
                }
            }
        });
        ui.add_space(4.0);

        let (rect, resp) = ui.allocate_exact_size(ui.available_size(), egui::Sense::click_and_drag());
        let painter = ui.painter_at(rect);
        painter.rect_filled(rect, 0.0, theme::BG_WIN);

        let tex = match self.layer {
            Layer::Dem => &self.tex_dem,
            Layer::DemShade => &self.tex_dem_shade,
            Layer::Terrain => &self.tex_terrain,
            Layer::Subclass => &self.tex_sub,
        };
        if let Some(tex) = tex {
            if resp.hovered() {
                let scroll = ui.input(|i| i.raw_scroll_delta.y);
                if scroll != 0.0 {
                    let f = if scroll > 0.0 { 1.2 } else { 1.0 / 1.2 };
                    let mouse = resp.hover_pos().unwrap() - rect.left_top();
                    let mx = (mouse.x - self.view_ox * self.view_scale) / self.view_scale;
                    let my = (mouse.y - self.view_oy * self.view_scale) / self.view_scale;
                    self.view_scale = (self.view_scale * f).clamp(0.008, 1.0);
                    self.view_ox = mouse.x - mx * self.view_scale;
                    self.view_oy = mouse.y - my * self.view_scale;
                }
            }
            if resp.dragged() {
                self.view_ox += resp.drag_delta().x / self.view_scale;
                self.view_oy += resp.drag_delta().y / self.view_scale;
            }
            let (tw, th) = self.data_wh;
            let draw_w = tw * self.view_scale;
            let draw_h = th * self.view_scale;
            let x0 = rect.left() + self.view_ox * self.view_scale;
            let y0 = rect.top() + self.view_oy * self.view_scale;
            painter.image(
                tex.id(),
                egui::Rect::from_min_size(egui::pos2(x0, y0), egui::vec2(draw_w, draw_h)),
                egui::Rect::from_min_max(egui::pos2(0.0, 0.0), egui::pos2(1.0, 1.0)),
                egui::Color32::WHITE,
            );

            // 浮动图例(左下); DEM 阴影为灰度无图例
            let items: Vec<(&str, egui::Color32)> = if matches!(self.layer, Layer::Dem | Layer::DemShade)
            {
                Vec::new()
            } else if self.layer == Layer::Terrain {
                theme::CLASS_COLORS
                    .iter()
                    .filter(|(c, _, _)| *c != 2)
                    .map(|(_, n, c)| (*n, *c))
                    .collect()
            } else {
                vec![
                    ("低丘", egui::Color32::from_rgb(180, 230, 150)),
                    ("高丘", egui::Color32::from_rgb(110, 195, 110)),
                    ("低山", egui::Color32::from_rgb(250, 225, 130)),
                    ("中山", egui::Color32::from_rgb(235, 170, 90)),
                    ("高山", egui::Color32::from_rgb(205, 110, 75)),
                    ("极高山", egui::Color32::from_rgb(150, 65, 60)),
                    ("平坝", egui::Color32::from_rgb(85, 185, 235)),
                ]
            };
            let lg_h = items.len() as f32 * 18.0 + 12.0;
            let lg_rect = egui::Rect::from_min_size(
                egui::pos2(rect.left() + 10.0, rect.bottom() - lg_h - 10.0),
                egui::vec2(126.0, lg_h),
            );
            if !items.is_empty() {
            painter.rect_filled(
                lg_rect,
                6.0,
                egui::Color32::from_rgba_unmultiplied(255, 255, 255, 242),
            );
            painter.rect_stroke(lg_rect, 6.0, egui::Stroke::new(1.0, theme::STROKE), egui::StrokeKind::Inside);
            for (i, (name, col)) in items.iter().enumerate() {
                let p = lg_rect.left_top() + egui::vec2(10.0, 9.0 + i as f32 * 18.0);
                painter.rect_filled(egui::Rect::from_min_size(p, egui::vec2(12.0, 12.0)), 3.0, *col);
                painter.text(
                    p + egui::vec2(18.0, 6.0),
                    egui::Align2::LEFT_CENTER,
                    *name,
                    egui::FontId::proportional(12.0),
                    theme::TEXT,
                );
            }
            }
        } else if self.dem_loading {
            // DEM 预览加载中
            painter.text(
                rect.center(),
                egui::Align2::CENTER_CENTER,
                "正在读取 DEM 并生成山体阴影预览…",
                egui::FontId::proportional(15.0),
                theme::TEXT_SUB,
            );
        } else if self.running {
            // 运行中大进度
            let c = rect.center();
            let bar_rect =
                egui::Rect::from_center_size(c + egui::vec2(0.0, 40.0), egui::vec2(380.0, 10.0));
            painter.rect_filled(bar_rect, 5.0, egui::Color32::from_rgb(226, 232, 240));
            let fill = bar_rect.width() * (self.last_pct / 100.0).clamp(0.0, 1.0);
            if fill > 1.0 {
                painter.rect_filled(
                    egui::Rect::from_min_size(bar_rect.left_top(), egui::vec2(fill, 10.0)),
                    5.0,
                    theme::ACCENT,
                );
            }
            painter.text(
                c + egui::vec2(0.0, -46.0),
                egui::Align2::CENTER_CENTER,
                format!("{}%", self.last_pct as i32),
                egui::FontId::proportional(42.0),
                theme::ACCENT,
            );
            painter.text(
                c + egui::vec2(0.0, -4.0),
                egui::Align2::CENTER_CENTER,
                &self.last_stage,
                egui::FontId::proportional(16.0),
                theme::TEXT,
            );
            painter.text(
                c + egui::vec2(0.0, 76.0),
                egui::Align2::CENTER_CENTER,
                "可点击右上角取消",
                egui::FontId::proportional(12.0),
                theme::TEXT_DIM,
            );
        } else {
            // 空状态引导
            let c = rect.center();
            let badge_c = c + egui::vec2(0.0, -110.0);
            painter.rect_filled(
                egui::Rect::from_center_size(badge_c, egui::vec2(72.0, 72.0)),
                16.0,
                theme::ACCENT.gamma_multiply(0.10),
            );
            let tri = [
                badge_c + egui::vec2(-20.0, 18.0),
                badge_c + egui::vec2(0.0, -16.0),
                badge_c + egui::vec2(18.0, 18.0),
            ];
            painter.add(egui::Shape::convex_polygon(
                tri.to_vec(),
                theme::ACCENT.gamma_multiply(0.85),
                egui::Stroke::NONE,
            ));
            painter.circle_filled(badge_c + egui::vec2(10.0, -20.0), 5.0, theme::WARN);
            painter.text(
                c + egui::vec2(0.0, -30.0),
                egui::Align2::CENTER_CENTER,
                "从 DEM 划分 8 类地形部位",
                egui::FontId::proportional(19.0),
                theme::TEXT,
            );
            let steps = [
                ("1", "导入 DEM", "float32 GeoTIFF, 投影坐标系"),
                ("2", "调整参数", "盆地宽度/坡度阈值等, 均可默认"),
                ("3", "运行划分", "输出 8 类部位 + 地貌亚类栅格"),
            ];
            for (i, (n, t, s)) in steps.iter().enumerate() {
                let p = c + egui::vec2(-210.0 + i as f32 * 170.0, 40.0);
                painter.circle_filled(p + egui::vec2(0.0, -14.0), 11.0, theme::ACCENT);
                painter.text(
                    p + egui::vec2(0.0, -14.0),
                    egui::Align2::CENTER_CENTER,
                    *n,
                    egui::FontId::proportional(12.0),
                    egui::Color32::from_rgb(10, 12, 16),
                );
                painter.text(
                    p + egui::vec2(0.0, 8.0),
                    egui::Align2::CENTER_CENTER,
                    *t,
                    egui::FontId::proportional(13.5),
                    theme::TEXT,
                );
                painter.text(
                    p + egui::vec2(0.0, 28.0),
                    egui::Align2::CENTER_CENTER,
                    *s,
                    egui::FontId::proportional(10.0),
                    theme::TEXT_DIM,
                );
            }
            painter.text(
                c + egui::vec2(0.0, 110.0),
                egui::Align2::CENTER_CENTER,
                "点击左侧「运行划分」开始 · 示例 DEM 已随包附带(sample/)",
                egui::FontId::proportional(12.0),
                theme::TEXT_DIM,
            );
        }
    }

    // ---------- 统计卡带 ----------
    fn stats_strip(&mut self, ui: &mut egui::Ui) {
        if self.stats.is_empty() {
            return;
        }
        ui.add_space(6.0);
        egui::Frame::new()
            .fill(theme::BG_CARD)
            .corner_radius(8.0)
            .inner_margin(egui::Margin::symmetric(10, 8))
            .show(ui, |ui| {
                let max_area = self
                    .stats
                    .iter()
                    .map(|s| s.area_km2)
                    .fold(0f64, f64::max)
                    .max(1e-9);
                ui.horizontal(|ui| {
                    for s in &self.stats {
                        let w = (ui.available_width() - 30.0) / self.stats.len() as f32;
                        egui::Frame::new()
                            .inner_margin(egui::Margin::symmetric(6, 6))
                            .show(ui, |ui| {
                                ui.set_min_width(w - 12.0);
                                ui.horizontal(|ui| {
                                    ui.label(egui::RichText::new("■").color(s.color).small());
                                    ui.label(egui::RichText::new(s.name).small());
                                });
                                ui.label(
                                    egui::RichText::new(format!("{:.1} km²", s.area_km2))
                                        .monospace()
                                        .strong(),
                                );
                                let frac = (s.area_km2 / max_area).clamp(0.02, 1.0);
                                let (r, _) = ui
                                    .allocate_exact_size(egui::vec2(w - 26.0, 4.0), egui::Sense::hover());
                                ui.painter().rect_filled(r, 2.0, theme::BG_INPUT);
                                ui.painter().rect_filled(
                                    egui::Rect::from_min_size(
                                        r.left_top(),
                                        egui::vec2(r.width() * frac as f32, 4.0),
                                    ),
                                    2.0,
                                    s.color,
                                );
                            });
                    }
                });
            });
        ui.add_space(2.0);
    }
}

impl eframe::App for App {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        // 持续重绘: egui 默认按需重绘, 文件对话框/resize 事件竞争会留下
        // 未重绘的陈旧缓冲区(表现为局部黑色条带); 每帧请求重绘彻底消除。
        ctx.request_repaint();
        self.poll_worker(ctx);
        self.poll_dem(ctx);

        // 顶栏
        egui::TopBottomPanel::top("top")
            .frame(egui::Frame::new().fill(theme::BG_PANEL).inner_margin(egui::Margin::symmetric(12, 8)).stroke(egui::Stroke::new(1.0, theme::STROKE)))
            .show(ctx, |ui| {
                ui.horizontal(|ui| {
                    let (r, _) = ui.allocate_exact_size(egui::vec2(30.0, 30.0), egui::Sense::hover());
                    ui.painter().rect_filled(r, 8.0, theme::ACCENT_DIM);
                    ui.painter().text(
                        r.center(),
                        egui::Align2::CENTER_CENTER,
                        "TP",
                        egui::FontId::proportional(15.0),
                        egui::Color32::WHITE,
                    );
                    ui.vertical(|ui| {
                        ui.set_min_width(180.0);
                        ui.label(egui::RichText::new("TerraPos").strong().size(17.0));
                        ui.label(
                            egui::RichText::new(format!(
                                "地形部位划分工具 · 西南区 8 部位 · 构建 {}",
                                env!("BUILD_STAMP")
                            ))
                            .small()
                            .color(theme::TEXT_SUB),
                        );
                    });
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        status_pill(ui, self.running, &self.last_stage, self.last_pct);
                        if self.running
                            && ui
                                .add(egui::Button::new("取消").fill(theme::BG_INPUT).corner_radius(6.0))
                                .clicked()
                        {
                            self.cancel.store(true, std::sync::atomic::Ordering::Relaxed);
                        }
                    });
                });
            });

        // 左侧参数
        egui::SidePanel::left("params")
            .exact_width(348.0)
            .frame(egui::Frame::new().fill(theme::BG_PANEL).inner_margin(egui::Margin::symmetric(8, 8)))
            .resizable(false)
            .show(ctx, |ui| {
                egui::ScrollArea::vertical().show(ui, |ui| {
                    self.params_panel(ui);
                });
            });

        // 底部日志
        egui::TopBottomPanel::bottom("log")
            .height_range(egui::Rangef::new(96.0, 200.0))
            .frame(egui::Frame::new().fill(theme::BG_PANEL).inner_margin(egui::Margin::symmetric(8, 6)))
            .show(ctx, |ui| {
                if !self.running {
                    if let Some(s) = &self.result_summary {
                        egui::CollapsingHeader::new(egui::RichText::new("面积统计报告").small())
                            .default_open(true)
                            .show(ui, |ui| {
                                egui::ScrollArea::vertical().max_height(90.0).show(ui, |ui| {
                                    ui.monospace(s);
                                });
                            });
                    }
                }
                egui::ScrollArea::vertical().stick_to_bottom(true).show(ui, |ui| {
                    for (l, is_err) in &self.log {
                        ui.monospace(
                            egui::RichText::new(l)
                                .color(if *is_err { theme::ERR } else { theme::TEXT_SUB }),
                        );
                    }
                });
            });

        // 中央: 统计卡带 + 预览
        egui::CentralPanel::default()
            .frame(egui::Frame::new().fill(theme::BG_WIN).inner_margin(egui::Margin::same(6)))
            .show(ctx, |ui| {
                if !self.running {
                    self.stats_strip(ui);
                }
                self.preview_panel(ui);
            });

        if self.request_start {
            self.request_start = false;
            if !self.params.dem_path.is_empty() && !self.running {
                self.start_run(ctx);
            }
        }
    }
}

fn colorize(class: &[u8]) -> Vec<u8> {
    let mut v = Vec::with_capacity(class.len() * 4);
    for &c in class {
        let col = theme::CLASS_COLORS
            .iter()
            .find(|(k, _, _)| *k == c)
            .map(|x| x.2)
            .unwrap_or(egui::Color32::from_rgb(12, 13, 16));
        v.extend_from_slice(&[col.r(), col.g(), col.b(), 255]);
    }
    v
}

fn colorize_sub(class: &[u8]) -> Vec<u8> {
    let table = [
        (1u8, [180u8, 230, 150]),
        (2, [110, 195, 110]),
        (3, [250, 225, 130]),
        (4, [235, 170, 90]),
        (5, [205, 110, 75]),
        (6, [150, 65, 60]),
        (7, [85, 185, 235]),
    ];
    let mut v = Vec::with_capacity(class.len() * 4);
    for &c in class {
        let col = table
            .iter()
            .find(|(k, _)| *k == c)
            .map(|x| x.1)
            .unwrap_or([12, 13, 16]);
        v.extend_from_slice(&[col[0], col[1], col[2], 255]);
    }
    v
}

// ============================== 参数详解数据 ==============================
struct ParamDoc {
    name: &'static str,
    default: &'static str,
    /// 该参数控制什么(地貌含义)
    what: &'static str,
    /// 调整方向与影响
    effect: &'static str,
}

const PARAM_DOCS: &[(&str, egui::Color32, &[ParamDoc])] = &[
    (
        "六级坡位与海拔分区",
        theme::SEC[2],
        &[
            ParamDoc {
                name: "六级坡位",
                default: "HAND 频率分位自适应",
                what: "坡位栅格分六级：1 山谷、2 坡下、3 平坡、4 坡中、5 坡上、6 山脊。分级阈值取 HAND 的 10%/35%/85%/95% 频率分位，随区域自适应；坡度小于 3° 的缓坡面判为平坡，坡度 20° 以上的陡坡上调一档。",
                effect: "无需调整；分级随数据自动适应区域地形起伏幅度。",
            },
            ParamDoc {
                name: "海拔分区",
                default: "500 / 800 / 1200 m",
                what: "叠加用的海拔分区：1 <500m(丘陵)、2 500-800m、3 800-1200m、4 ≥1200m。坝子落在海拔 1/2 区为宽谷盆地、第 3 区为山间盆地。",
                effect: "固定口径，与参考脚本一致，不作为界面参数。",
            },
            ParamDoc {
                name: "坡位叠加分档",
                default: "S1,2=下 S3,4=中 S5,6=上",
                what: "最终组合时六级坡位的归档规则：山谷+坡下归「下部」，平坡+坡中归「中部」，坡上+山脊归「上部」，再与丘陵(海拔一区)/山地相乘得到部位 3-8。",
                effect: "固定口径，与参考脚本一致，不作为界面参数。",
            },
            ParamDoc {
                name: "阶地修正",
                default: "坡度<2° 且 HAND<80m",
                what: "同时满足时判为坡下——对应规则中坡下包含低阶地、高阶地、坡麓、河漫滩与谷底排水线的口径。",
                effect: "收紧（坡度阈值更小/HAND 更小）则阶地更多归入坡中；放宽则河谷两侧平缓带更多划入坡下。",
            },
        ],
    ),
    (
        "山间盆地",
        theme::SEC[1],
        &[
            ParamDoc {
                name: "识别体系（河谷低平带+对象检验）",
                default: "legacy 同参数",
                what: "坝子 = 河流沿岸的冲积平地。主干与您参考脚本 generate_500_area 同参数:                       DEM 河网锚定河谷(替代地类河流矢量), 500m 缓冲内 且 与河高差<5m                       且 坡度<5° 且 起伏<5m; 程序增益: 碎片桥接重组 + 对象级内部起伏                       检验(排除大起伏混入) + 填洞平滑。",
                effect: "hhgq 实测 87.2km²/1749 个(您参数版 65.8km²/935 个), 自然河谷覆盖更全。",
            },
            ParamDoc {
                name: "河网阈值",
                default: "1 km²",
                what: "DEM 提取河网的汇流面积门槛, 决定坝子沿哪一级河谷分布。                       等效您原方法的地类河流矢量(河流/湖泊/水库/沟渠)。",
                effect: "减小 → 支谷沟谷沿岸也出坝子(更全面更碎)；增大 → 只沿干流大坝子。",
            },
            ParamDoc {
                name: "河流缓冲",
                default: "500 m",
                what: "距最近河流的水平距离上限(对齐 legacy buffer_distance=500)。",
                effect: "增大 → 坝子向谷坡扩展；减小 → 紧贴水线。",
            },
            ParamDoc {
                name: "河流高程差上限",
                default: "5 m",
                what: "与最近河流的高差上限(HAND), 即与谷底同高(对齐 elevation_threshold=5)。",
                effect: "增大 → 更高的阶地/台地并入；减小 → 只留贴水线窄带。",
            },
            ParamDoc {
                name: "坡度上限",
                default: "5°",
                what: "逐像元坡度门槛(5m 层 Horn, 对齐 slope_threshold=5)。",
                effect: "增大 → 缓坡带并入；减小 → 只留平地。",
            },
            ParamDoc {
                name: "局部起伏上限",
                default: "5 m",
                what: "5×5 窗(25m)内高差上限(对齐 relief_threshold=5, 5m 层语义)。",
                effect: "增大 → 微起伏滩地并入；减小 → 只留极平区。",
            },
            ParamDoc {
                name: "最低面积",
                default: "100 亩",
                what: "坝子保留的面积下限(归并与平滑后按连通域计)。不足该面积的                       坝子整体转出：海拔一区(<500m)归丘陵下部，其余归山地坡下                       ——不再以碎斑形式留在坝子内部。",
                effect: "增大 → 只留大坝子(小坝转坡下)；减小 → 小坝子保留。",
            },
            ParamDoc {
                name: "碎斑归并",
                default: "100 m / 30 亩",
                what: "坝子域(桥接闭运算圈定)内的小面积非坝子碎斑(细碎坡下/坡中)                       整体并入坝子，使坝子完整连片(对齐 find_and_fill_hole 语义)。                       半径 0 = 关闭归并。",
                effect: "增大 → 更大通道内的碎斑也并入(坝子更整)；                         碎斑面积上限增大 → 更大碎斑也吞并。",
            },
            ParamDoc {
                name: "内部起伏上限",
                default: "15 m",
                what: "对象内高程 P95−P95 差的上限——程序增益判据, 防止整体高差大的                       地块(如包含陡坡的假平地)混入坝子。",
                effect: "增大 → 容纳更大倾角；减小 → 只要近水平地面。",
            },
            ParamDoc {
                name: "碎片桥接 / 平滑距离",
                default: "50 m / 50 m",
                what: "桥接: 被沟渠道路切碎的坝子重组(对齐 smoothing_buffer);                       平滑: 定型后闭运算消锯齿(对齐 smoothing_distance)。",
                effect: "增大 → 形态更连片圆滑；减小 → 更贴判据细节。",
            },
        ],
    ),
    (
        "丘陵与山地",
        theme::SEC[3],
        &[
            ParamDoc {
                name: "丘陵海拔上限",
                default: "500 m",
                what: "海拔低于该值为丘陵(再按起伏度分低丘/高丘)，\
                       不低于则为山地(按绝对高程分低山/中山/高山/极高山)。",
                effect: "本区海拔 626~1714m 全部为山地；换区域时若存在\
                         低海拔丘陵带会自动出现丘陵类。",
            },
            ParamDoc {
                name: "亚类起伏度窗口",
                default: "2000 m",
                what: "低丘/高丘分界所用的起伏度计算窗口，\
                       覆盖一个完整丘包所需的邻域尺度。",
                effect: "增大 → 起伏度趋大、更多划为高丘。",
            },
            ParamDoc {
                name: "低丘起伏度上限",
                default: "200 m",
                what: "相对高差低于该值为低丘，200~500m 为高丘（规则图口径）。",
                effect: "增大 → 更多划为低丘。",
            },
        ],
    ),
    (
        "后处理",
        theme::SEC[4],
        &[
            ParamDoc {
                name: "众数滤波轮数",
                default: "1",
                what: "3×3 众数滤波消除椒盐状碎斑，同票时保留原值以保护窄长地物。",
                effect: "增大 → 更平滑但窄带可能被抹除；0 = 关闭。",
            },
            ParamDoc {
                name: "最小图斑",
                default: "400 像元(约 1 公顷@5m)",
                what: "小于该面积的图斑按最近邻并入周边大图斑。",
                effect: "增大 → 结果更整洁；减小 → 保留更多细碎图斑。",
            },
        ],
    ),
];

fn app_icon() -> Option<std::sync::Arc<egui::IconData>> {
    let png = include_bytes!("../assets/icon-256.png");
    let img = image::load_from_memory(png).ok()?.into_rgba8();
    let (w, h) = img.dimensions();
    Some(std::sync::Arc::new(egui::IconData {
        width: w,
        height: h,
        rgba: img.into_raw(),
    }))
}

fn install_panic_logger() {
    let exe = std::env::current_exe()
        .ok()
        .and_then(|p| p.parent().map(|d| d.join("TerraPos-错误日志.txt")))
        .unwrap_or_else(|| std::path::PathBuf::from("TerraPos-错误日志.txt"));
    std::panic::set_hook(Box::new(move |info| {
        use std::io::Write;
        let t = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_secs())
            .unwrap_or(0);
        if let Ok(mut f) = std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(&exe)
        {
            let _ = writeln!(f, "[unix {t}] panic: {info}");
        }
    }));
}

fn main() -> eframe::Result {
    install_panic_logger();
    let native = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1560.0, 940.0])
            .with_min_inner_size([1100.0, 700.0])
            .with_title("TerraPos 地形部位划分工具")
            .with_icon(app_icon().expect("app icon")),
        ..Default::default()
    };
    eframe::run_native(
        "TerraPos",
        native,
        Box::new(|cc| Ok(Box::new(App::new(cc)))),
    )
}
