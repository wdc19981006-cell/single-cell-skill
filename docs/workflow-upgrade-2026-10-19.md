# 2026-10-19 修改与验收记录

本次在现有模块上修改；未重建项目，未安装或更新 R/Python 包。代码保留在
`D:/CodexProjects/single-cell-skill`。默认结果根目录由
`runtime/output_config.json` 唯一定义为
`D:/Desktop/codex-result/single-cell-skill`，数据集名称只追加一次。

## 已实现

- Python/R 共用存储契约；下载、QC、正式分析、科学审计和运行日志接入新根目录。
  `--output-root` 支持独立存储，显式 `--root OLD_WORKSPACE` 兼容旧 results/data；
  `--input` 读取旧 RDS。非 GEO 名称通过路径和 Windows 设备名称检查。
- `--run-id` 将替代方案隔离到 `03_analysis/runs/ID`，不会自动发现其他运行的
  Harmony、聚类或 checkpoint。已有最终 RDS、源矩阵和确认回执不覆盖。
- QC 阈值、全局基因 >3 个保留细胞、DoubletFinder 核心及周期评分不变。
  增加逐样本规则命中、单一规则独占、条件重叠和去重删除总数的 CSV/PDF。
  超过 15% 的预检查确认机制保留。4% 仍明确为双细胞预期参数。
- DecontX 保留 celda 路线和 RNA，校正矩阵不取整、不按污染比例删细胞。
  增加 P90/P95/P99、UMI/nFeature 分布、Marker 变化、异常警告及只读复核入口。
  真空液滴背景仅接受带文库映射与证据的稀疏整数矩阵；不自动生成背景。
- 正式表达默认 decontX，默认仅回归 pMT；缺失/无效校正矩阵明确停止。
  Harmony 默认 sample；混杂记录限制并保存未整合结果，不虚构逐次人工确认。
  缺失、无效或单水平批次显式跳过；group 禁用。清除重跑输入的旧整合状态。
- 基线聚类后按风险执行 feature/count 对照；周期风险先复核增殖 Marker。
  每个实际候选固定输入哈希、细胞、assay、HVG、种子、PC 和分辨率，并执行
  PCA/Harmony/邻居图/聚类，输出子空间、ARI、转移、Marker、周期及样本/分组证据。
  无风险和已回归变量不做重复计算。建议不代表生物学正确性证明。
- 有重大差异时，保留当前正式方案需用户明确记录 `regression_scheme=baseline`
  及原话 `regression_statement`；采纳候选需新 run ID 和匹配比较证据，重新评估
  PC/分辨率。分辨率仍由用户选择，生成该值专属 Marker/QC/周期/组成图。
- 审计区分 PASS/WARNING/ACTION_REQUIRED/RESOLVED/FAILED/MISSING_EVIDENCE。
  已验证 finalize 与文件校验和一致后，历史分辨率待办标记 RESOLVED；已选择的
  回归取舍同步处理。重复警告汇总展示，当前待办不包含已解决记录。
- 各入口提供简短 JSON 摘要。Windows 状态原子替换短暂锁定时有限重试，持久失败
  保留旧状态和 pending；验证/发布中断可复用成功计算，原生崩溃仍为失败。
  未命名 R 图形设备保存在对应调用日志旁，避免仓库产生 Rplots.pdf。
- 三个 Skill、README 和配置说明区分常规、定向复核、开发调试模式，固定脚本承担计算。

## 实际验证范围

| 验证 | 结果与边界 |
|---|---|
| 定向 Python 回归 | 104 项通过，0 失败、0 跳过；覆盖下载兼容、输入保护、审计、路径隔离、人工决策和 Windows 恢复 |
| 固定 R 环境 | R 4.5.0；实际 Seurat/SeuratObject 5.4.0、Matrix 1.7.3、DoubletFinder 2.0.6；环境未变更 |
| 分阶段 QC 回归 | 实际 PCA/过滤/周期评分及图表生成通过；此项的 DoubletFinder/ambient API 为显式模拟，不冒称真实 DF 计算 |
| 真实 celda 调用 | 160 细胞小型合成对象实际调用 celda::decontX；小数 counts 的 Seurat 5 存取、RNA/metadata/细胞一致性通过 |
| 默认正式分析 | 同一真实校正矩阵实际执行 sample Harmony（含混杂提示）、聚类、0.3 专属图、审计与无重算续跑；用户选择均标注 TEST FIXTURE |
| 额外回归 | 自然 count 风险触发完整对照；feature 分支用明确标注的零诊断阈值 fixture 验证；周期分支用模拟 Marker 复核决定验证；候选均未覆盖基线 |
| 方案采纳与 RNA 替代 | 独立 run 中验证用户授权采纳的 PCA/PC 重评估，另建 RNA PCA；未声称完成采纳方案的真实数据完整分析 |
| 实际原始构建 | 新根目录下生成 250 基因 ×4 细胞 10x fixture，实际构建、重读验证、PDF、报告复用及拒绝覆盖测试通过；没有新 GEO 网络下载 |
| 图表验证 | 第一轮 59 文件/116 PDF 页面；最终路径修复后独立分析 46 文件/90 PDF 页面；真实只读复核 1 PDF/36 页面；逐页渲染与 PNG DPI 检查通过，并抽看关键图 |
| 摘要与状态 | 最终 8 份阶段 JSON 的输入/输出路径实际存在，临时 request 不列入正式产物；最新审计有 5 条历史决定正确 RESOLVED |
| Skill | 三个 SKILL.md 的 quick_validate.py 通过（Windows 下启用 Python UTF-8 模式） |

## GSE149614

只读复核已完成的历史 QC-3。新报告保存为：

`D:/Desktop/codex-result/single-cell-skill/GSE149614/02_QC/03_ambient_RNA/reviews/upgrade-review-20261019/`

`review_provenance.json` 记录输入 SHA256，结束后再次检查一致。GSM4505956 被结构通用
规则识别为 EXTREME_CONTAMINATION，其 P95 约 0.991；细胞保持不变。这只是异常诊断，
不是删除或判定过度校正的结论。历史正式分析对象、用户决定和审计未被覆盖。

未无条件重算 GSE149614，也未宣称其 decontX/RNA、额外 nFeature 回归、G2M 增殖解释
或最终 0.3 小簇的完整生物学复核已完成。新的真实正式分析及其人工选择仍需按独立 run
执行。真实空液滴 background 路线因没有对应背景证据，未在真实数据上运行。

## 复现入口

在 Git Bash 中从代码仓库调用 `/c/Python312/python.exe`：

```text
tests/run_upgrade_units.py --output <new-results>/TEST_upgrade/validation/unit_results.json
tests/run_upgrade_acceptance.py --fixture <small-actual-corrected.rds> --dataset TEST_unique_name
tests/run_download_layout_regression.py --output-root <isolated-collection-root>
runtime/r45/run_r45.py tests/test_analysis.R D:/CodexProjects/single-cell-skill
qc/review_ambient.py GSE149614 --input <existing-QC3.rds> --review-id <new-review-id>
```

验收入口使用新目录；不覆盖已完成 fixture。正式用法及方案决定 schema 见
[analysis/README.md](../analysis/README.md)。日志中的 future 包构建版本提示未引起失败，
没有据此升级包。开发结果提交本地 Git；本次未发起真实 GEO 下载审计上传。

## 修改文件清单

```text
.agents/skills/geo-single-cell-loader/SKILL.md
.agents/skills/geo-single-cell-loader/references/audit-policy.md
.agents/skills/geo-single-cell-loader/scripts/audit_run.py
.agents/skills/geo-single-cell-loader/scripts/build_manifest.py
.agents/skills/geo-single-cell-loader/scripts/build_seurat.R
.agents/skills/geo-single-cell-loader/scripts/common.py
.agents/skills/geo-single-cell-loader/scripts/download_processed.py
.agents/skills/geo-single-cell-loader/scripts/inspect_geo.py
.agents/skills/geo-single-cell-loader/scripts/prebuild_probe.py
.agents/skills/geo-single-cell-loader/scripts/run_confirmed.py
.agents/skills/geo-single-cell-loader/scripts/sample_info.py
.agents/skills/geo-single-cell-loader/scripts/seurat_common.R
.agents/skills/single-cell-analysis/SKILL.md
.agents/skills/single-cell-qc/SKILL.md
AGENTS.md
README.md
analysis/README.md
analysis/analysis_utils.R
analysis/clustering.R
analysis/default_config.json
analysis/harmony.R
analysis/pca_regression.R
analysis/regression_review.R
analysis/run_analysis.R
analysis/run_workflow.py
analysis/validate_analysis.R
docs/workflow-upgrade-2026-10-19.md
qc/ambient_review.R
qc/ambient_rna.R
qc/default_config.yaml
qc/qc_filter.R
qc/qc_precheck.R
qc/qc_utils.R
qc/review_ambient.R
qc/review_ambient.py
qc/run_qc.R
qc/run_workflow.py
runtime/output_config.json
runtime/project_paths.R
runtime/project_paths.py
runtime/r45/README.md
runtime/r45/r45_profile.R
runtime/r45/run_r45.py
runtime/r45/runtime.py
tests/run_analysis_sensitivity.py
tests/run_download_layout_regression.py
tests/run_upgrade_acceptance.py
tests/run_upgrade_units.py
tests/test_analysis.R
tests/test_analysis_workflow.py
tests/test_download_layout.py
tests/test_output_paths.py
tests/test_qc_workflow.py
tests/test_runtime.py
tests/test_upgrade.R
tests/test_workflow_audit.py
tests/verify_staged_figures.py
workflow_audit/render_audit.R
workflow_audit/run_audit.py
```
