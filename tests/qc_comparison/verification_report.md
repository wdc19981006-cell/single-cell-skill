# QC PCA 修复与真实数据一致性验证

日期：2026-10-16（Asia/Shanghai）。参考仓库 `wdc19981006-cell/data-QC` 固定版本：`3c11adff25e461758eaca3e67743a4afa3e7f28e`。源码基线：`ec37df1`。

四组历史对象已完整 readRDS；八个真实文件来自本地镜像，读取前后均与 Git LFS SHA-256/大小核对一致。未把 LFS pointer 当作 RDS。历史比较先完成，之后修改 QC，再仅复跑 GSE167297。

## 四组历史对象：原代码 A vs 旧 Skill B

| 数据集 | 细胞 A / B | 共同细胞 | 仅 A / 仅 B | 基因 A / B | 共同基因 | 仅 A / 仅 B 基因 | counts 不同元素 |
|---|---:|---:|---:|---:|---:|---:|---:|
| GSE149614 | 54802 / 54931 | 53100 | 1702 / 1831 | 24194 / 24863 | 24194 | 0 / 669 | 0 |
| GSE167297 | 22312 / 22374 | 21566 | 746 / 808 | 21317 / 22026 | 21317 | 0 / 709 | 0 |
| GSE183904 | 135301 / 135589 | 130886 | 4415 / 4703 | 25092 / 25583 | 25092 | 0 / 491 | 0 |
| GSE189903 | 116996 / 117193 | 113152 | 3844 / 4041 | 25156 / 25847 | 25156 | 0 / 691 | 0 |

counts 使用最多 512 个细胞一块的 sparse matrix 精确比较；从未将整体表达矩阵转为 dense。全量差集 ID 保存在各 `*_cell_differences.csv` 和 `*_gene_differences.csv`，逐 sample/group 数量、结构和 metadata 比较见 [object_comparison.csv](object_comparison.csv) 与 [object_comparison.md](object_comparison.md)。

- 四组均为 Seurat / RNA Assay5，实际 A 和 B 都只有 counts layer，VariableFeatures、reductions、graphs、neighbors 均为空。A 的 metadata 仍有临时聚类列；B 删除这些列并增加 doublet_status 和 cell-cycle 字段。
- **EXPECTED_DIFFERENCE**：A 用 DF 过滤，B 用 DF_adj 过滤；A 额外做基因出现细胞数 >3 的过滤；A 回归 nFeature_RNA 和 pMT，Skill 只回归 pMT。已证明 nCount/nFeature 差异完全来自非共有基因行，各自均与 sparse counts 列统计吻合。共同细胞的 sample、group、pMT、pRP、pHB 一致。
- **UNEXPECTED_DIFFERENCE**：实际上传的四个 A 对象均没有 S.Score、G2M.Score、Phase、CC.Difference，也没有 reference 代码最后阶段应保存的 normalized data/PCA。因此无法对 A 做 cell-cycle 分数一致性比较。这是历史产物与代码描述不一致，未伪造缺失字段，也未修改历史对象。
- **UNEXPECTED_DIFFERENCE（本轮修复）**：旧 PC 实现取第一个 elbow drop，GSE167297 全部样本只用 2 PCs。不同 retained cell IDs 不能仅用“DF vs DF_adj”解释；PC 错误、regression 策略和随机算法均可能影响分类。

| 数据集 | A 最终 DF=Singlet | B 最终 DF=Doublet 但 DF_adj=Singlet | B 最终 DF_adj=Singlet |
|---|---:|---:|---:|
| GSE149614 | 54802 | 129 | 54931 |
| GSE167297 | 22312 | 62 | 22374 |
| GSE183904 | 135301 | 288 | 135589 |
| GSE189903 | 116996 | 197 | 117193 |

这些是过滤后对象的统计，不能反推已丢弃细胞的全部 DF 标签；完整共同细胞交叉表另存 CSV。

## PCA 与绘图实现

选择策略仍使用标准差占比 `stdev/sum(stdev)*100`：cumulative >90 且单项 <5 的首次 PC，与最后一个 drop >0.1 的位置 +1 取较小值；单个条件有效时使用该值，无有效条件时才 fallback。动态上限不超过 50、cells-1、HVG-1、genes-1，并保留较小用户配置上限及原有小对象保护。

真实方差解释率 `stdev^2/sum(stdev^2)*100` 只作为诊断，分母是本次计算的 PCs。报告明确区分标准差占比和方差解释率；至少有 10 PCs 但只选 2–3 时记录 PCA_LOW_DIM_WARNING。本次真实运行没有此警告。

保持逐 sample DoubletFinder、pN=0.25、expected_rate=0.04、pMT-only regression、DF_adj 过滤。没有增加逐 sample UMAP，没有改变 cell-cycle 评分/回归策略。最终 RDS 不保留 PCA/UMAP。

QC violin 改为 Seurat VlnPlot，固定五指标顺序、3 列、独立且 before/after 一致的样本色、20×8 inches 基础尺寸。doublet UMAP 三 panel 共用坐标和轴范围，固定比例；Doublet 红色 #ea5455、Singlet 蓝色 #3f72af，Doublet 在上层。新增 singlet 图只使用过滤前 embedding 的实际保留细胞子集，NotEvaluated 不改名。cell-cycle 图形逻辑未改。

57 样本的合成布局检查发现并修复底部图例裁切：超过 25 样本时纵向排列图例，并为两 panel 图增加宽度；重新渲染后标签和图例完整。该条件不影响本次 14 样本的真实计算。

## GSE167297：原代码 A、旧 Skill B、修复 Skill C

| 指标 | A | B | C |
|---|---:|---:|---:|
| 最终细胞 | 22312 | 22374 | 22417 |
| 最终基因 | 21317 | 22026 | 22026 |
| Cell QC 后细胞 | 历史最终 RDS 不含该阶段记录 | 23243 | 23243 |
| 每 sample PCs | 历史 RDS 未保存 | 全部 2 | 8–23 |
| 合并 PCs | 历史 RDS 未保存 | 旧报告未记录 | 17 / max 50 |
| 最终过滤 | DF | DF_adj | DF_adj |

C 合并 variance_pick=42（标准差占比阈值），drop_pick=17，selected cumulative true variance=78.7108%。不会为了让累计方差达到 90% 或让图更好看而强行增加 PC。

| sample | B PC | C PC | B DF_adj doublets | C DF_adj doublets |
|---|---:|---:|---:|---:|
| GSM5101013 | 2 | 9 | 15 | 15 |
| GSM5101014 | 2 | 8 | 63 | 61 |
| GSM5101015 | 2 | 11 | 72 | 69 |
| GSM5101016 | 2 | 13 | 48 | 46 |
| GSM5101017 | 2 | 14 | 82 | 78 |
| GSM5101018 | 2 | 10 | 45 | 41 |
| GSM5101019 | 2 | 13 | 146 | 139 |
| GSM5101020 | 2 | 23 | 116 | 110 |
| GSM5101021 | 2 | 13 | 17 | 16 |
| GSM5101022 | 2 | 11 | 99 | 91 |
| GSM5101023 | 2 | 13 | 79 | 75 |
| GSM5101024 | 2 | 11 | 32 | 31 |
| GSM5101025 | 2 | 9 | 16 | 16 |
| GSM5101026 | 2 | 14 | 39 | 38 |

B/C 未调整 DF 的总 Doublet 数均为 931（按相同样本量和 4% 的期望数）；homotypic 调整后从 869 变为 826，故 C 多保留 43 个细胞。分类对象并非相同集合：B/C 共同 21,642，B-only 732，C-only 775；共同细胞中 DF 标签有 149 个不同。最终 DF_adj 均为 Singlet。具体 PC、阈值、真实累计方差和每 sample doublet 数见 [GSE167297_pc_doublet_comparison.csv](GSE167297_pc_doublet_comparison.csv)。

A/C 共同 22,076 个细胞、21,317 个基因；A-only 236 个细胞，C-only 341 个细胞；共同 counts 全部相等，见 [A/C 对照](fixed_vs_original/object_comparison.md)。B/C 基因完全相同，共同 counts 全部相等；C 全部细胞/基因 counts 与 raw 对应子集精确一致。

raw 输入为 27,073 个细胞，其 SHA-256 与旧 QC 运行记录完全一致。现有用户决定按原字节复用，无新增阈值/分组决策。sample/group/database 逐细胞不变；S.Score、G2M.Score、Phase、CC.Difference 均通过有效性和关系检查；没有按 Phase 删除细胞。C 仅有 RNA counts + metadata，无 HVG、reductions、graphs、neighbors 或临时聚类列。

raw SHA-256（历史、源文件、测试副本、运行前后相同）：`ab7c98333f65e3f42ff79f0b138b68dcc8720110e7cb1b0995534a918d9601b2`。C RDS SHA-256：`302f1528b3b7ac10100c4f24d88c8dda1539d761760e463ebeaa6f935fa7bae6`。

新产物仅位于 `.runtime/qc-pca-fixed/data/GSE167297/qc/`；八个历史 RDS 的 SHA-256 在结束后重新验证未变。未重新运行另外三个真实数据集。

## 图像检查

四组各五类图已经并排渲染并视觉检查。旧 Skill 的指标排列是字母序、单色且标签拥挤，四组 doublet UMAP 都出现细线/蛇形。原代码的 singlet 图在多样本图例下受挤压；旧 Skill 无对应图。部分原始 FeatureStatPlot 图例把红/蓝分别对应 Singlet/Doublet，与本次明确要求的命名配色相反；新版使用明确的命名映射，不复制该视觉问题。

GSE167297 C 的 merged UMAP 恢复成团块结构，没有旧版密集蛇形。该视觉观察与 PC 选择修复、相同 Cell QC 细胞量及 counts 不变共同构成验证，不能单独证明细胞类型或 DoubletFinder 分类的生物学真值。保留 pMT-only regression；没有据图形擅自恢复 nFeature_RNA regression。

所有新图七对 PDF/PNG 已解码/渲染，物理尺寸与 300 dpi 校验通过；逐一检查无标签裁切。cell-cycle 两图也已渲染检查，算法和样式保持原流程。

| 数据集 | before QC | after QC | doublet UMAP | doublet violin | singlet |
|---|---|---|---|---|---|
| GSE149614 | [查看](previews/GSE149614_before_QC.png) | [查看](previews/GSE149614_after_QC.png) | [查看](previews/GSE149614_doublet_umap.png) | [查看](previews/GSE149614_doublet_vlnplot.png) | [查看](previews/GSE149614_singlet_dimplot.png) |
| GSE167297 | [查看](previews/GSE167297_before_QC.png) | [查看](previews/GSE167297_after_QC.png) | [查看](previews/GSE167297_doublet_umap.png) | [查看](previews/GSE167297_doublet_vlnplot.png) | [查看](previews/GSE167297_singlet_dimplot.png) |
| GSE183904 | [查看](previews/GSE183904_before_QC.png) | [查看](previews/GSE183904_after_QC.png) | [查看](previews/GSE183904_doublet_umap.png) | [查看](previews/GSE183904_doublet_vlnplot.png) | [查看](previews/GSE183904_singlet_dimplot.png) |
| GSE189903 | [查看](previews/GSE189903_before_QC.png) | [查看](previews/GSE189903_after_QC.png) | [查看](previews/GSE189903_doublet_umap.png) | [查看](previews/GSE189903_doublet_vlnplot.png) | [查看](previews/GSE189903_singlet_dimplot.png) |

GSE167297 预览按 A / B / C 三列展示；其他三组只有 A / B。全分辨率真实新图保留在独立测试输出目录。

## 测试与执行证据

- Python unittest：86/86 PASS。
- 完整 run_regressions.py：PASS，含已有 Loader 回归、QC synthetic、真实安装 DoubletFinder、已有产物恢复校验。
- 最终源码的 focused QC：PCA、动态 max_pc、小样本、逐 sample DF、merged 图、DF_adj、violin、singlet 同坐标子集、NotEvaluated、cell-cycle、clean object、raw hash 全部 PASS。
- R 4.5.0 固定 Runtime：真实 QC 和独立验证均 exit 0；Seurat 5.4.0、SeuratObject 5.4.0、Matrix 1.7.3、DoubletFinder 2.0.6。未安装、更新或切换 R 包。
- 7 对新图的 PDF/PNG 与 300 dpi 验证 PASS；57 样本 synthetic 排版渲染 PASS；Skill 格式校验 PASS；git diff --check PASS。
- 首次历史比较遇到 C locale/中文路径，改用 SHA 一致的 ASCII 副本解决；初次 QC 回归在测试文件仍被编辑时解析失败，冻结后完整重跑通过；首次额外真实比较因缓存缺旧报告停止，补齐后只读验证通过。所有失败均保留日志，没有把失败退出当 PASS。

详见 [测试结果](test_results.json)、[真实验证](GSE167297_real_validation.json)、[新 QC 报告](GSE167297_fixed_qc_report.csv)、[图像校验](GSE167297_figure_validation.csv)、[数据保护](data_protection.json)。

复现命令（从仓库根目录，经固定 Python/R45；真实运行使用独立 root）：

```text
/c/Python312/python.exe -m unittest discover -s tests -p "test_*.py"
/c/Python312/python.exe D:/CodexProjects/single-cell-skill/tests/run_regressions.py
/c/Python312/python.exe D:/CodexProjects/single-cell-skill/tests/qc_comparison/run_qc_tests.py
/c/Python312/python.exe D:/CodexProjects/single-cell-skill/qc/run_workflow.py GSE167297 --root D:/CodexProjects/single-cell-skill/.runtime/qc-pca-fixed
```

真实输出已存在时上述 QC 命令只进行验证，不重建覆盖。源码提交不包含 RDS、表达矩阵或其他大型原始数据。审计为独立 QC，Loader discovery/download/build 均 NOT_APPLICABLE；十个允许的文本文件保存于审计仓库 `GSE167297/20261016-qc-pca-validation/`。
