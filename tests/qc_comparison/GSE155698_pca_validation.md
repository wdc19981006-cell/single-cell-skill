# GSE155698 PCA 自动选择验证

基线 GitHub main：`cec9f452e8fe704c187b2058eefd195987d1a4e9`。本次仅新增测试观察脚本、真实诊断表、图形和文档；没有修改生产 QC、PCA 方法、Loader、下载流程、R45 Runtime 或 regression 策略。

当前规则的 20 个样本选择 **9–22 PCs**，合并选择 **23 PCs / max 50**。本次真实运行没有重现全样本 2-PC 错误。

## 输入与版本身份

原始输入为 `data/GSE155698/seurat_raw.rds`，55,854 cells × 26,630 genes。SHA-256：`8b97d5f5019bdca4c08319581b9866f0e5c649643b550b0a38494327a4f1fad7`，与此前工作流记录一致。原始 RDS、历史 RDS、报告、group_confirmation 和用户决策均在前后核对哈希，未覆盖。

- 历史原代码：根目录 `seurat.qc.rds` / `seurat.all_qc.rds` 和原始 PDF。前者与当前 Cell QC 的 40,915 个细胞 ID 完全相同，但只有 26,086 个基因；当前保留 26,630 个基因。原代码额外过滤 544 个基因，最终按 DF 过滤。
- 已保存的 `qc/`：已属于修复后的 Skill，20 样本 PCs 为 9–22、合并 23；不能称作修复前 2-PC Skill。它作为此前修复版本的复现对照。
- 真正修复前的 GSE155698 Skill RDS/UMAP 未找到。因此增加 **legacy_pc_control：旧 PC 规则受控复现**，只在测试进程中恢复 `ec37df1` 的第一个 drop 和 >=90 条件，其他现代流程相同。它不是历史完整旧 Skill 产物。
- current：独立测试目录执行当前固定 QC。两组受控运行均复制同一 raw、复用原有明确用户决策、使用同一样本顺序与 seed=155698、相同阈值和 pMT-only regression。没有构造新的用户接受记录。

## 算法核对

RunPCA 上限为 min(50, 配置上限, cells-1, HVGs-1, genes-1)，最终 pc_used 不超过 max_pc；现有小样本保护保留。选择使用标准差占比 `stdev/sum(stdev)*100`，cumulative >90 且 individual <5 的首次位置，与最后一个 drop >0.1 的位置 +1 取较小值。单个有效条件保留，无有效条件才 fallback，最终至少 2 PCs。

标准方差解释率使用 `stdev^2/sum(stdev^2)*100`，只作独立诊断。此处分母为实际计算的前 max_pc 个 PC 的平方标准差之和，不是全部基因总方差。ElbowPlot 的纵轴明确为 PCA standard deviation。

已从完整导出谱独立重算 21 个 PCA 选择，全部与当前 R 方法一致。两组受控运行的 PCA 标准差谱最大绝对差为 0；选择规则的差异可与谱本身分开检查。

合并对象：variance_pick=41（历史命名，实际是标准差累计规则）、drop_pick=23、pc_used=23。所选标准差累计占比 66.6010%，诊断方差累计解释率 85.1198%。旧规则受控合并选择 2 PCs。

## 每个样本及合并 PCA

cells 是 Cell QC 后、DF 前实际进入 PCA 的细胞数；原始每样本细胞数、group、HVG 数和完整指标另见 [诊断 CSV](GSE155698/current_pca_diagnostics.csv)。所有样本和 ALL 的单张 ElbowPlot 均已保存 PDF/300-dpi PNG，PNG 和完整谱随报告提交。

| sample | cells | max_pc | pc_used | variance_pick | drop_pick | 累计 SD % | 累计方差 % | ElbowPlot |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| GSM4710689 | 866 | 50 | 13 | 43 | 13 | 43.828 | 66.334 | [查看](GSE155698/elbows/elbow_GSM4710689.png) |
| GSM4710690 | 1512 | 50 | 12 | 43 | 12 | 42.957 | 68.608 | [查看](GSE155698/elbows/elbow_GSM4710690.png) |
| GSM4710691 | 1518 | 50 | 16 | 42 | 16 | 52.396 | 77.575 | [查看](GSE155698/elbows/elbow_GSM4710691.png) |
| GSM4710692 | 1636 | 50 | 16 | 42 | 16 | 52.136 | 74.470 | [查看](GSE155698/elbows/elbow_GSM4710692.png) |
| GSM4710693 | 913 | 50 | 10 | 43 | 10 | 36.733 | 59.262 | [查看](GSE155698/elbows/elbow_GSM4710693.png) |
| GSM4710694 | 1388 | 50 | 11 | 43 | 11 | 38.913 | 61.082 | [查看](GSE155698/elbows/elbow_GSM4710694.png) |
| GSM4710695 | 5906 | 50 | 14 | 41 | 14 | 55.052 | 82.406 | [查看](GSE155698/elbows/elbow_GSM4710695.png) |
| GSM4710696 | 3365 | 50 | 16 | 42 | 16 | 55.505 | 80.163 | [查看](GSE155698/elbows/elbow_GSM4710696.png) |
| GSM4710697 | 3913 | 50 | 20 | 42 | 20 | 61.455 | 82.797 | [查看](GSE155698/elbows/elbow_GSM4710697.png) |
| GSM4710698 | 1165 | 50 | 22 | 43 | 22 | 58.949 | 74.903 | [查看](GSE155698/elbows/elbow_GSM4710698.png) |
| GSM4710699 | 3105 | 50 | 12 | 43 | 12 | 40.271 | 62.128 | [查看](GSE155698/elbows/elbow_GSM4710699.png) |
| GSM4710700 | 1515 | 50 | 14 | 43 | 14 | 45.174 | 66.522 | [查看](GSE155698/elbows/elbow_GSM4710700.png) |
| GSM4710701 | 1121 | 50 | 17 | 43 | 17 | 52.195 | 74.298 | [查看](GSE155698/elbows/elbow_GSM4710701.png) |
| GSM4710702 | 3644 | 50 | 17 | 42 | 17 | 54.577 | 77.926 | [查看](GSE155698/elbows/elbow_GSM4710702.png) |
| GSM4710703 | 614 | 50 | 13 | 43 | 13 | 42.920 | 63.463 | [查看](GSE155698/elbows/elbow_GSM4710703.png) |
| GSM4710704 | 1609 | 50 | 12 | 43 | 12 | 42.417 | 67.228 | [查看](GSE155698/elbows/elbow_GSM4710704.png) |
| GSM4710705 | 650 | 50 | 10 | 43 | 10 | 34.186 | 55.770 | [查看](GSE155698/elbows/elbow_GSM4710705.png) |
| GSM4710706 | 3642 | 50 | 16 | 41 | 16 | 57.338 | 81.980 | [查看](GSE155698/elbows/elbow_GSM4710706.png) |
| GSM4710707 | 2397 | 50 | 15 | 42 | 15 | 52.336 | 77.125 | [查看](GSE155698/elbows/elbow_GSM4710707.png) |
| GSM4710708 | 436 | 50 | 9 | 43 | 9 | 35.341 | 61.275 | [查看](GSE155698/elbows/elbow_GSM4710708.png) |
| ALL | 40915 | 50 | 23 | 41 | 23 | 66.601 | 85.120 | [查看](GSE155698/elbows/elbow_ALL.png) |

## 每样本 DoubletFinder

DoubletFinder 按 sample 独立运行；观察脚本在每次 paramSweep 时断言只有一个 sample，且无逐样本 UMAP。merged 不运行 DoubletFinder，因此 selected_pK 为不适用，DF/DF_adj 只汇总逐样本结果。

| sample | selected_pK | DF Singlet / Doublet | DF_adj Singlet / Doublet |
|---|---:|---:|---:|
| GSM4710689 | 0.02 | 831 / 35 | 835 / 31 |
| GSM4710690 | 0.01 | 1452 / 60 | 1457 / 55 |
| GSM4710691 | 0.01 | 1457 / 61 | 1463 / 55 |
| GSM4710692 | 0.18 | 1571 / 65 | 1576 / 60 |
| GSM4710693 | 0.02 | 876 / 37 | 880 / 33 |
| GSM4710694 | 0.02 | 1332 / 56 | 1338 / 50 |
| GSM4710695 | 0.005 | 5670 / 236 | 5690 / 216 |
| GSM4710696 | 0.29 | 3230 / 135 | 3242 / 123 |
| GSM4710697 | 0.21 | 3756 / 157 | 3768 / 145 |
| GSM4710698 | 0.24 | 1118 / 47 | 1122 / 43 |
| GSM4710699 | 0.07 | 2981 / 124 | 2991 / 114 |
| GSM4710700 | 0.13 | 1454 / 61 | 1463 / 52 |
| GSM4710701 | 0.3 | 1076 / 45 | 1082 / 39 |
| GSM4710702 | 0.02 | 3498 / 146 | 3510 / 134 |
| GSM4710703 | 0.02 | 589 / 25 | 593 / 21 |
| GSM4710704 | 0.01 | 1545 / 64 | 1550 / 59 |
| GSM4710705 | 0.28 | 624 / 26 | 628 / 22 |
| GSM4710706 | 0.28 | 3496 / 146 | 3506 / 136 |
| GSM4710707 | 0.01 | 2301 / 96 | 2310 / 87 |
| GSM4710708 | 0.21 | 419 / 17 | 421 / 15 |
| ALL | 不适用 | 39276 / 1639 | 39425 / 1490 |

## 同细胞集合 UMAP 与原始 counts 验证

current 与 legacy_pc_control 的过滤前细胞集合完全相同，共 40,915 cells × 26,630 genes，样本和 group 对齐。DF 标签改变 2,920 个细胞，DF_adj 标签改变 2,700 个细胞；全量过滤前交叉表另存 CSV。

| 版本 | 最终 cells | 最终 genes | 与 raw 共同坐标 counts 差异 | clean object |
|---|---:|---:|---:|---|
| original | 39276 | 26086 | 0 | FALSE |
| prior_fixed | 39425 | 26630 | 0 | TRUE |
| current | 39425 | 26630 | 0 | TRUE |
| legacy_pc_control | 39383 | 26630 | 0 | TRUE |

| 对照 vs current | 共同 cells | 仅对照 cells | 仅 current cells | 共同 genes | 仅对照 genes | 仅 current genes |
|---|---:|---:|---:|---:|---:|---:|
| original | 38952 | 324 | 473 | 26086 | 0 | 544 |
| prior_fixed | 39425 | 0 | 0 | 26630 | 0 | 0 |
| legacy_pc_control | 38054 | 1329 | 1371 | 26630 | 0 | 0 |

当前最终 RDS 与此前修复版的序列化 SHA-256 完全相同：`7459f377b990986ce1cda7dd83dab765aa6796ab0d62b197280db17f6b9645a7`。本次相同哈希证明 counts 和 metadata 的最终序列化结果均可复现。

当前 DF_adj Doublet 的 nFeature_RNA 中位数为 2901，Singlet 为 1407。这只作为表达量辅助检查，不代表有双细胞生物学真值标签。

所有 counts 比较均使用最多 256 个细胞一块的 sparse matrix，未将整体表达矩阵转为 dense。两组新对象的全部基因和保留细胞 counts 均与 raw 对应子集精确一致。cell-cycle 前后保留细胞集合一致；只评分，不按 Phase 删除、不回归 S/G2M。最终对象不保留 normalized data、VariableFeatures、PCA、UMAP、graph、neighbor 或临时 cluster 列。

过滤前图明确标识为 Doublet UMAP；过滤后图明确标识为 Singlet UMAP，并使用各自过滤前坐标的保留子集。测试实际核对坐标与细胞 ID，没有为 singlet 重新拟合 UMAP。

过滤后不同规则保留的集合不同，因此主要几何比较另取两组共同保留的 38,054 个细胞。不能把两组不同细胞集合的最终图直接当作 PCA 优劣证据。原代码的 PCA/UMAP 坐标矩阵未保存，且基因过滤/regression 与 Skill 不同，其 PDF 只作历史背景，不作纯 PCA 因果对照。

经渲染并人工检查，current 的 UMAP 没有旧 PC 规则受控图中的全局密集线状/蛇形退化。局部细长分支本身不自动视为错误；UMAP 形态不是判断 PCA 科学有效性的唯一依据。实现验证同时依赖独立谱重算、动态边界、相同输入/细胞集合、原始 counts 和完整 DF_adj 比较。该启发式验证不证明存在唯一生物学最优 PC 数。

当前正式规则未出现需要改变 regression 的异常，因此没有启动 nFeature_RNA+pMT 敏感性分析，也未替换 pMT-only 策略。

- [原代码 / 旧 PC 受控复现 / 当前：过滤前](GSE155698/pre_filter_three_versions.png)
- [原代码 / 旧 PC 受控复现 / 当前：实际过滤后集合](GSE155698/post_filter_three_versions.png)
- [过滤前相同细胞：sample](GSE155698/pre_filter_matched_sample.png)
- [过滤前相同细胞：DF_adj](GSE155698/pre_filter_matched_DF_adj.png)
- [过滤后共同细胞：sample](GSE155698/post_filter_common_sample.png)
- [过滤后共同细胞：DF_adj](GSE155698/post_filter_common_DF_adj.png)
- [20 样本 ElbowPlot 总览](GSE155698/sample_elbows_overview.png)
- [合并 ElbowPlot](GSE155698/elbows/elbow_ALL.png)

## 实际测试与未解决限制

- Python tests：86/86 PASS；现有 PCA elbow、动态上限、小样本、fallback、标准差/方差分离和同坐标子集测试 PASS，R exit 0。
- 两组真实完整 QC、独立最终对象验证均 exit 0；使用安装的 DoubletFinder 2.0.6 和 Seurat 5.4.0。未安装、升级或切换任何 R 包。
- 固定 R 4.5.0 Runtime。21 个导出谱的独立计算、全量稀疏 counts、metadata、DF_adj、cell-cycle 和 clean object 检查 PASS。
- 原始输入和全部受保护历史文件哈希前后相同。诊断不重建 raw，不更改 group，不覆盖既有 RDS。
- 未找到真正修复前 GSE155698 Skill 历史产物；受控旧 PC 规则复现不能替代其历史 provenance。原代码坐标未保存，无法在原代码坐标上重新绘制精确共同细胞子集。
- 生产代码没有发现需要修复的新问题；本次提交只增加真实验证工具、诊断图表和文档。

[测试证据](GSE155698/test_evidence.json) · [独立选择重算](GSE155698/independent_selection_review.json) · [对象/集合验证](GSE155698/validation.json)

## 复现

从仓库根目录，通过 `/c/Python312/python.exe` 执行 `tests/qc_comparison/pca_validation/run_validation.py GSE155698 SOURCE NEW_ROOT`。旧规则受控复现使用另一个 NEW_ROOT 并加 `--mode legacy_pc_control`。目录必须全新；已有最终对象不被覆盖。R 计算由当前 qc/run_workflow.py 调用固定 Runtime，测试观察脚本只导出诊断并调用正式 qc_run()。

全部新 RDS、每张 ElbowPlot 的 PDF 和完整诊断坐标保留在 `.runtime/pca-validation-GSE155698/`。源码提交仅包含小型诊断表/图及测试文档。
