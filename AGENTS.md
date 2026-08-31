# AGENTS.md（项目级）

## 语言与称呼

- 始终使用简体中文回复。
- 回复开头使用 "Manba" 称呼。

## 构建发布流程（必须严格执行）

**触发条件**：对 `rust/` 下源代码、`python/` 脚本、文档或图标做出任何功能修改或修复后，必须完整执行以下流程，同步远程仓库代码与 GitHub Releases 产品。版本号保持 `v0.0.1` 不变（除非 Manba 明确要求升版）。

1. **质量门**：
   ```bash
   cd rust
   cargo test --release          # 全部测试通过(当前 26 项)
   cargo clippy --release        # 零警告
   cargo build --release
   cargo run -p topo_core --example run_sample --release   # 端到端冒烟, 确认分类分布合理
   ```
2. **重打便携包**（仓库根目录执行）：
   - 目录 `TerraPos-v0.0.1-win64/`，内含（**不含 sample/ 目录**, Manba 已要求移除）：
     - `TerraPos.exe` ← `rust/target/release/topo_app.exe`
     - `README.txt`（含默认分割策略等要点, 不提 sample）
     - `RELEASE_NOTES.txt` ← `docs/RELEASE_NOTES-v0.0.1.md`
   - 打包：`powershell -NoProfile -Command "Compress-Archive -Path 'TerraPos-v0.0.1-win64' -DestinationPath 'dist\TerraPos-v0.0.1-win64.zip' -Force"`（本机无 zip 命令，用 PowerShell）
   - 打包后删除临时目录 `TerraPos-v0.0.1-win64/`
   - 校验：zip 内 `TerraPos.exe` 大小应随构建变化，勿复用旧 zip
3. **提交推送**：
   ```bash
   git add -A && git commit && git push origin main
   ```
4. **更新 release 资产**（clobber 覆盖，不新建 release）：
   ```bash
   gh release upload TerraPos-v0.0.1 dist/TerraPos-v0.0.1-win64.zip --clobber
   gh release edit TerraPos-v0.0.1 --draft=false   # 关键: 删过 tag 后 release 会变 Draft(公众不可见), 必须转正式
   ```
5. **验证**：`gh release view TerraPos-v0.0.1` 确认 `draft=false`、资产时间戳/大小已更新；`git log origin/main` 确认代码同步。

**注意**：
- exe 图标与 Windows GUI 子系统已配置（`topo_app/build.rs` + `#![windows_subsystem = "windows"]`），图标源文件在 `rust/topo_app/assets/`，修改图标后运行 `python scripts/make_icon.py` 重新生成，再走上述流程。
- 若 Manba 要求"重置远程记录/保持初始提交"：用 orphan 分支将历史 squash 为单提交后 `git push --force`，删除并重建 tag `TerraPos-v0.0.1`，再走第 4 步。
- 本机无 `zip` 命令；Python 用 `D:/worker_code/.venvgis/Scripts/python.exe`（含 GDAL/matplotlib/PIL）。
