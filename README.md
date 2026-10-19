# Single-cell workflow Skills

项目使用现有 **R 4.5.0**，分为下载、三步 QC、正式分析三个独立阶段。固定 R 脚本负责运算，Skill 负责调用、解释和需要用户判断的决定。集中配置 `runtime/output_config.json` 决定默认结果根目录，所有 Python 入口支持 `--output-root`，代码仓库保持原位。非 GEO 数据可使用经过安全校验的名称（如 Study_A）。新任务采用以下布局；历史 `data/<GSE>/` 原地兼容，已有确认回执不改写，原始数据和有效 RDS 不覆盖。

```text
D:/Desktop/codex-result/single-cell-skill/<GSE>/
├── 01_download/               raw、seurat_raw.rds、metadata、下载统计/日志/PDF
├── 02_QC/
│   ├── 01_basic_QC/           seurat_QC-1.rds、细胞/基因过滤记录、前后图、配置/日志
│   ├── 02_doublet_cellcycle/  seurat_QC-2.rds、DF/DF_adj、周期评分、图表/配置/日志
│   └── 03_ambient_RNA/        seurat_QC-3.rds、decontX assay、污染评估/比较
└── 03_analysis/
    ├── seurat_pca.rds
    ├── seurat_harmony.rds     仅在设计检查允许且实际运行 Harmony 时生成
    ├── seurat_analysis.rds    用户选定已审阅分辨率后生成
    ├── 01_pca_regression/
    ├── 02_harmony/
    ├── 03_clustering/
    └── 99_workflow_audit/
```

已有数据可直接请求“对 GSE… 进行质控”；已有 QC 对象可请求“对 GSE… 做 PCA、Harmony 设计评估和聚类审阅”。技术人员独立启动示例（Git Bash）：

```bash
/c/Python312/python.exe D:/CodexProjects/single-cell-skill/qc/run_workflow.py GSE123 --precheck
/c/Python312/python.exe D:/CodexProjects/single-cell-skill/qc/run_workflow.py GSE123 --step all
/c/Python312/python.exe D:/CodexProjects/single-cell-skill/qc/run_workflow.py GSE123 --step qc3 --input D:/path/seurat_QC-2.rds
/c/Python312/python.exe D:/CodexProjects/single-cell-skill/analysis/run_workflow.py GSE123 --input D:/path/seurat_QC-3.rds
/c/Python312/python.exe D:/CodexProjects/single-cell-skill/analysis/run_workflow.py GSE123 --step finalize --resolution-decision D:/path/decision.json
/c/Python312/python.exe D:/CodexProjects/single-cell-skill/workflow_audit/run_audit.py GSE123
```

默认阈值保持不变；QC-1 新增保留细胞中表达数 >3 的基因过滤，并重算当前 counts 对应的 nCount/nFeature。QC-2 按实际 capture library 执行 DoubletFinder，只按 `DF_adj` 去双细胞，细胞周期只评分。QC-3 不按固定污染阈值删除细胞，也不覆盖 RNA counts。正式分析默认用 decontX 和 pMT 回归，原始 RNA/QC metadata 保留；缺失或无效校正矩阵明确停止。完成聚类后仅按实际风险运行 nFeature/nCount 对照；细胞周期对照先复核增殖 Marker。Harmony 默认按 sample，混杂时记录解释限制并保存未整合结果，不把 group 当校正变量，也不重复询问默认字段。分辨率扫描复用同一 UMAP，最终标签需用户明确选择。

替代分析使用 `--run-id decontx-v1`，结果写入该数据集的 `03_analysis/runs/decontx-v1/`；历史 `results/` 不移动、不删除，使用 `--input` 显式读取，或以 `--root OLD_WORKSPACE` 恢复旧布局。常规执行优先读取阶段 JSON 摘要，只有实际异常才读取日志；定向复核使用固定脚本，开发调试只运行受影响测试。详见 [正式分析配置与运行说明](analysis/README.md)。

科学审计读取既有统计和警告，不无条件重跑分析；缺失证据明确列出，已有报告在证据变化后版本化保存。它与下载阶段私有仓库的十文件审计分开。包缺失只停止相关步骤；不安装包、不升级 R。

实现细节和本次真实验证见 [三阶段升级记录](docs/staged-workflow-validation.md)。以下说明保留下载器的输入结构、确认和审计契约。

在 Codex 中打开 `single-cell-skill` 项目根目录，然后直接输入：

```text
下载单细胞数据 GSE231993
```

也支持“下载 GSE231993”“处理 GSE231993”“读取 GSE231993”“分析单细胞数据 GSE231993”“帮我下载并整理 GSE229413”，以及英文 “Download and prepare GSE231993”。普通用户不需要输入 Skill 名称、阶段名称，也不需要手动运行 Python/R 脚本。

1. Skill 搜索 GEO Series、全部 GSM、关联论文和公开样本 metadata。
2. 返回完整样本属性：sample、author_sample、tissue、disease、source_type、sample_description；有可靠证据时补充 patient、specimen、cohort、treatment。
3. 告诉用户“group 尚未创建，请确认分组方式。”确认前不批量下载表达矩阵。
4. 用户确认后，保存分组记录、生成并校验 manifest、自动下载源表达文件。
5. 创建并独立重新读取验证 Seurat 对象，更新结果说明。
6. 新下载结果保存在 `D:/Desktop/codex-result/single-cell-skill/<GSE>/01_download/`；历史 `data/<GSE>/` 原地兼容。Loader 以 `COMPLETE_RAW` 结束并询问是否继续 QC；未获用户确认不会自动质控。

## GEO MCP-first Stage A

每一次新的 Stage A discovery 都必须首先实际尝试用户级全局 MCP `geo`。关键词检索依次使用 `search_geo`、候选 GSE 的 `get_geo_info`、`list_geo_files`；已知 GSE 直接使用 `get_geo_info` 和 `list_geo_files`。只有完整 Series 结果确实缺少个别 GSM 细节时才调用 `get_geo_sample`，不会对全部 GSM 逐个重复请求。只有 accession 正确、`complete=true`、`sample_count` 与 GSM inventory 一致、source URL 存在，且文件具有明确的 owner/source level/candidate format 时，才接受 MCP 结果。MCP 只提供公开 metadata 与文件 URL，分组、下载、manifest 和 Seurat 构建仍全部由本 Skill 负责。

MCP 成功后，Stage A 直接用同一套 structured public facts 写入 `inspection.json`、`sample_report.csv` 和 `sample_info.txt`，不会再运行 `inspect_geo.py` 重复查询。只有真实工具不可用、启动/请求错误或 MCP 结果无法补全时才 fallback，并必须通过 `--fallback-reason` 记录实际失败原因。`inspection.json` 用 `discovery_method`、`mcp_tools_used`、`mcp_complete`、起止时间和秒数标识来源；MCP 为 `geo_mcp`，fallback 为 `inspect_geo_fallback`。不能仅因工具尚未调用或老会话未暴露工具就声称 MCP 不可用。

Stage A 还会执行 single-cell modality gate：明确的 microarray 或 bulk RNA-seq 会在表达矩阵下载和 Seurat 构建前停止；只有用户明确要求“即使不是单细胞也下载文件”时，才允许普通文件下载，且不会进入 Seurat 构建。

## Pooled/shared matrix handling

manifest 中多个 sample 可以共享一个 pooled expression matrix，但必须共享同一个经过确认的 `cell_map_path`。构建器用包含路径、格式、counts source、分隔符、方向、feature/drop columns 和 assay 的 read signature 检查读取语义；相同路径的配置不一致立即停止。每个唯一物理表达输入只调用一次 reader，每个 shared cell map 也只读取一次。

共享矩阵读取后，cell map 的 cell 集合必须与矩阵列精确一致，sample 集合必须与该输入对应的 manifest samples 精确一致。构建器按矩阵列顺序一次性生成 `sample_barcode`，保留一个完整 pooled sparse matrix；不会按 sample 重复读取、不会拆成多个完整矩阵再拼回。多个 pooled 或独立输入在 feature 严格对齐后才合并，最终仍只调用一次 `CreateSeuratObject(min.cells=3, min.features=200)`。

## Input Router 与标准输出

| Input format | Reader | Mapping strategy |
|---|---|---|
| 10x trio | prebuild 结构检查后 `Seurat::Read10X()` | 每个 physical trio 检查/读取一次；缺 gene name 时可用唯一 feature ID 修复 |
| 10x H5 | schema 检查后 `Seurat::Read10X_h5()` | 单样本前缀或显式 cell map |
| H5AD / H5AD.gz | 优先 `zellkonverter::readH5AD(reader="R")`，否则仅使用显式配置的 Python anndata | 必须明确 `X`、`raw` 或 counts layer；多样本使用 cell map |
| Seurat RDS | `readRDS()` 后只提取指定 RNA counts | 丢弃作者对象其余状态，重新建对象 |
| CSV / TSV / TXT / TXT.GZ | 小矩阵 `data.table::fread()`；dense 估算达到物理 RAM 25% 时 Python streaming | prebuild 记录 genes/cells、dense bytes、RAM 和 reader 选择原因 |
| pooled/shared | 复用底层 reader | 整个矩阵读取一次、cell map 读取一次，不按样本拆矩阵 |
| expression + separate metadata | expression 走自身 reader | metadata mapper 用明确 cell/sample key 映射 |
| ZIP / TAR / TGZ | 只提取 manifest 指定 exact member | 提取结果返回真实格式 reader，不建立 archive reader |
| FASTQ / SRA | V1 不支持定量 | 停止并报告需要 Cell Ranger/STARsolo/alevin |

Stage B 的所有 reader 都返回同一 contract，并在 `assert_counts()` 后形成 genes × cells 的 raw integer sparse counts。所有 input 完成 feature reconciliation 后，只调用一次 `CreateSeuratObject()`。

**Standard output：**无论来源是 10x、H5AD、RDS、text 或 pooled matrix，新任务最终写为 `D:/Desktop/codex-result/single-cell-skill/<GSE>/01_download/seurat_raw.rds`（历史目录不迁移）：一个 clean Seurat 对象、一个 RNA raw counts layer、标准 sample/cell metadata；不继承作者 PCA、UMAP、SCT/integrated assay、graph、neighbor、cluster 或 normalized expression。

## 独立单细胞 QC

已有原始 RDS 直接调用 `qc/run_workflow.py`，不重复下载。先 precheck；低细胞样本、物种/阈值冲突和预计删除比例 >15% 等情况需用户决定。持久化真实原话和实际受影响 sample ID，之后复用已有决定。默认包括 QC-1/2/3；`--step` 支持独立恢复，`--input` 可复用已有步骤输入。低细胞样本经用户同意保留但跳过 DF 时标记 NotEvaluated，不冒称 Singlet。

三步结果、阈值、DecontX counts 保留方式见文首及 [QC Skill](.agents/skills/single-cell-qc/SKILL.md)。下载完成仍需用户同意进入 QC；后续分析另由 [正式分析 Skill](.agents/skills/single-cell-analysis/SKILL.md) 调用。所有图使用 `qc/plot_utils.R` 的统一尺寸、主题、PDF 与 300 DPI PNG，复用相同 embedding 的坐标。

## R Runtime

正式环境为 Windows 11 + R 4.5.0，安装位置 `D:/R/R-4.5.0`；唯一环境定义位于 `runtime/r45/config.py`。所有 Loader、Raw Builder、QC、测试及临时 R 验证都经过 `runtime/r45/run_r45.py`，用户不需要手动调用 Rscript.exe。Python 负责 GEO discovery/download/routing；R 负责 Seurat build/validation/QC。

启动器固定 R 和 library，只在子进程设置 profile、`CLI_NO_THREAD=1` 和 `LC_ALL=C`。分析期间使用真实架构，仅退出 cleanup 阶段选择已验证的 cli workaround。正常 R 错误保留原始非零退出码；139 / Windows access violation 记录为 `R_NATIVE_CRASH`，即使输出 PASS 也不会被视为成功。没有 R 环境 fallback，不安装、升级或删除包，不修改系统设置。

`--healthcheck` 输出 `.runtime/r45/health.json`；`--route qc1|qc|ambient|analysis_pca|harmony|clustering|text|10x_h5|h5ad_native` 只要求对应的额外包。每个独立 workflow 至少检查一次。每次 R 调用的环境、阶段、脚本、stdout/stderr、退出码和时长写入数据集目录的 `.runtime/r45/`（开发临时检查仍可写仓库 `.runtime/`）。临时验证使用 `--expr "sessionInfo()"`。旧 `qc/r450_rscript.sh` 保留为兼容入口。

下载完成、原始对象完成、QC 完成分别记录为 `DOWNLOAD_COMPLETE`、`COMPLETE_RAW`、`COMPLETE_QC`。中断恢复先校验已有 bytes/provenance 或最终 RDS；不覆盖 final、不重复构建已完成大矩阵。`.workflow/state.json` 在审计清理前记录 checkpoint，原始 manifest/receipt 和状态随十文件审计保存；清理后可通过 `audit_run.py resume GSE` 恢复。已有 QC final 使用独立验证恢复。详见 [Runtime 说明](runtime/r45/README.md) 与 [历史退出问题记录](docs/r450-cli-exit-139.md)。

在样本资料确实支持这些类别时，用户可回复：

```text
健康人为 HC；UC患者非炎症组织为 UC_control；UC炎症组织为 UC。
```

这是确认方式示例，不是预设分组。group 只能由用户决定，不能根据 disease、文件名、Tumor/Normal 字样自动产生。PDAC 的 uninvolved 样本仍有 PDAC 疾病背景，不能因此自动归入 Tumor 或 Healthy。

## 一个 GSE，一个目录

新目录树见文首。第一阶段保留源文件 `raw/`、`seurat_raw.rds`、`sample_info.txt`、`group_confirmation.json`、必要的 `cell_map/`，以及 `metadata.csv`、`download_summary.csv`、`download_report.pdf`、`download.log`。这里 raw 表示未经本流程修改的下载源表达文件，不一定是 FASTQ。

`.workflow/` 暂存发现、回执、manifest、下载来源和运行证据。真实下载/构建结束后将十个指定文本审计文件推送到 [single-cell-skill-audit](https://github.com/wdc19981006-cell/single-cell-skill-audit) 的 `<GSE>/<run-id>/`；成功后只清理临时工作目录，保留阶段结果。失败或中断同样上传；上传失败则保留临时文件。矩阵和 RDS 不进入审计仓库。科学审计 PDF/CSV 保留在本地 `03_analysis/99_workflow_audit/`。Git 忽略 results 和 data 下运行数据；历史数据包保持原布局。

## 数据集运行示例

GSE149614 的 pooled count 矩阵与 21 个 GSM 映射流程见 [运行示例](examples/gse149614/README.md)。示例脚本只使用用户确认的 group；运行数据保存在本地 `data/GSE149614/`，不进入 Git 历史。

## 已发布数据

[GSE231993 v2：修正全局 gene filtering，UC 8 个样本、HC 4 个样本](https://github.com/wdc19981006-cell/single-cell-skill/releases/tag/gse231993-uc8-hc4-v2-20260914)（2026-09-14）：完整数据包包含 36 个 GEO 源表达文件、`seurat_raw.rds`、`sample_info.txt` 和本次运行的 `.workflow` 核心记录，保留 `data/GSE231993/` 布局，附有 SHA256 校验清单。

| 用户确认的 group | 样本 | 最终细胞数 |
|---|---|---:|
| UC | GSM7307094–GSM7307101，8 个样本 | 37,967 |
| HC | GSM7307102–GSM7307105，4 个样本 | 22,698 |

此处 UC 按用户要求同时包含作者的 UC-self control 和 UC 炎症样本，原始来源描述仍完整保留。v2 对象共 60,665 个细胞、25,953 个基因、305,606,810 UMI 和 58,767,679 个非零表达项；仅使用约定的构建阈值，未进行额外 QC、归一化或下游分析。构建、序列化重读及独立 metadata 验证均通过。

[旧版 Release](https://github.com/wdc19981006-cell/single-cell-skill/releases/tag/gse231993-uc8-hc4-20260914) 保留为历史记录，其中 60,665 cells × 23,183 genes 的 RDS 在每个 GSM 内分别执行了 `min.cells=3`，现已由 v2 取代。样本身份、cell barcode、group 和 metadata 映射没有改变。

## sample_info.txt

完整检索后即创建，不必等分组确认。第一行是：

```text
STATUS: WAITING_FOR_GROUP_CONFIRMATION
```

文件包括研究标题/说明、测序类型、组织、候选 processed 格式、总样本数，以及逐样本属性。尚未由证据确定的内容明确标注；整套资料没有的 optional 字段不会显示大量 NA。公开事实经 Skill 审阅后同步更新 TXT。

确认分组后状态为 `GROUP_CONFIRMED`，逐样本增加 group，并增加 GROUP SUMMARY。创建、序列化并验证原始对象成功后状态为 `COMPLETE_RAW`，附加 PROCESSING SUMMARY：总样本/细胞/基因数、sample/group 细胞数、输出位置、实际读取格式和 reader、counts source、Seurat 版本、创建时间和 warnings。QC 的状态独立记录在 `02_QC/` 及步骤报告中，不会把 Loader 状态改成 QC 完成。

构建失败时尽可能写入 `BUILD_FAILED` 和 `Failure:` 原因，保留已经下载的数据。确认前仍保持等待状态。部分开发检索保存在 `.workflow/development/`，不能覆盖完整报告，也不能用于后续分组确认。

## 数据与验证规则

manifest 是 metadata 和文件映射的唯一真源。最终 metadata 必须包含 `database/sample/tissue/disease/source_type/group`。可选 patient/specimen/cohort/treatment 每个非空值需要可靠 evidence。barcode 在合并前添加 sample 前缀；按稳定 cell/sample key 映射，验证 `rownames(meta.data) == Cells(seurat)`、样本集合和全部 manifest metadata。

确认回执绑定 manifest 内容 hash；修改后重新确认。下载记录保存 URL、SHA256、bytes、UTC timestamp；重复运行先验证后复用，checksum、来源或 archive member 不一致立即停止，绝不覆盖。归档成员也有独立 SHA256 校验。中断后可利用逐文件保存的 provenance 继续。`files_json=[]` 仅供明确审阅的本地输入，不能用来跳过下载校验。

支持 10x 三联（gzip/plain/genes.tsv）、10x H5、H5AD/H5AD.gz、明确方向的 TXT/CSV/TSV counts、Seurat RDS counts，以及 manifest 精确指定的 tar/zip 成员。text route 使用 `data.table::fread`，并继续强制 delimiter、orientation、feature_column、drop_columns 和整数 counts 校验；`data.table` 缺失不会影响非 text reader。同一样本优先 filtered，多份 ambiguous filtered 则停止。H5AD 不会交给 Read10X_h5；normalized matrix 不能当 raw counts。pooled matrix 必须有精确 cell-to-sample mapping。路径必须是 repository-relative、没有 `..`、不越出所属 GSE；expression 在 `raw/`，映射文件在 `.workflow/`。

多个样本先验证 feature 集合、统一顺序，并在原始稀疏 counts 矩阵层面合并，再统一执行一次 `CreateSeuratObject(counts=counts, min.cells=3, min.features=200)`。因此 `min.cells` 针对整个 GSE 数据集，而不是在每个 GSM 内分别过滤。feature 集合真正不同时停止，不能静默取交集或补零取并集。`min.features` 过滤后按最终 cell key 映射 sample，保持 `orig.ident = sample`。

这些构建阈值会过滤部分细胞/基因，summary 如实记录。没有额外 QC，不执行 NormalizeData、SCTransform、FindVariableFeatures、ScaleData、PCA/UMAP、邻居图、聚类、DoubletFinder、细胞注释、差异或富集分析。只有 FASTQ/SRA 的数据需要重新定量，超出当前范围。

## Advanced / Developer usage

以下命令供维护和诊断使用，日常工作由 Skill 内部执行。脚本从自身位置解析 repository root，Python 支持 `--root`；R 显式接收 repository root。

依赖：Python 3.10+（检索、manifest、下载只用标准库）；R + Seurat、SeuratObject、Matrix、jsonlite；H5 需要 hdf5r。H5AD 优先 zellkonverter、SingleCellExperiment、SummarizedExperiment；缺失时显式设置 `GEO_SINGLE_CELL_PYTHON` 指向已验证安装 anndata/numpy/scipy 的解释器。不会自动安装环境。合成 fixture 另需 h5py、pandas。先定位实际 Python/R，不使用 Windows Store 的 python/python3 占位程序。

```bash
ROOT=D:/CodexProjects/single-cell-skill
PY=/c/Python312/python.exe
SCRIPTS="$ROOT/.agents/skills/geo-single-cell-loader/scripts"
WORKFLOW="$ROOT/results/GSE231993/01_download/.workflow"
"$PY" "$SCRIPTS/audit_run.py" start GSE231993
"$PY" "$ROOT/runtime/r45/run_r45.py" --healthcheck
# 仅在真实 geo MCP 调用失败时使用，并记录实际原因：
"$PY" "$SCRIPTS/inspect_geo.py" GSE231993 --fallback-reason "geo/get_geo_info tool error: <actual error>"
# Skill 审阅证据、完善 sample_report.csv 与 inspection.json 后：
"$PY" "$SCRIPTS/sample_info.py" GSE231993
# 用户确认 group，并保存 group_confirmation.json 后：
"$PY" "$SCRIPTS/build_manifest.py" --report "$WORKFLOW/sample_report.csv" --confirmation "$WORKFLOW/group_confirmation.json" --output "$WORKFLOW/sample_manifest.csv"
"$PY" "$SCRIPTS/build_manifest.py" --validate "$WORKFLOW/sample_manifest.csv"
"$PY" "$SCRIPTS/run_confirmed.py" GSE231993
```

测试：

```bash
ROOT=D:/CodexProjects/single-cell-skill
PY=/c/Python312/python.exe
"$PY" -m unittest discover -s "$ROOT/tests" -p 'test_*.py' -v
"$PY" "$ROOT/tests/runtime/run_tests.py" --real-fixture
"$PY" "$ROOT/tests/run_regressions.py"
```

合成数据使用保留测试编号 `GSE999999999`，分组批准明确标注 SIMULATED；不对应真实临床组。`tests/run_regressions.py` 每次创建独立合成 workspace，不覆盖已有测试或真实 RDS。旧布局数据需按 GSE 迁移并比较大小/hash；保留原始 manifest/receipt，不自动修改 hash 来认可新路径。无法确定归属的文件保留并报告。

详见 [metadata contract](.agents/skills/geo-single-cell-loader/references/metadata-schema.md)、[format routing](.agents/skills/geo-single-cell-loader/references/format-routing.md) 和 [实际验证报告](tests/validation-report.md)。公开元数据仍需证据审阅，任意作者私有格式或超内存矩阵可能需要额外适配。Skill 已启用隐式调用，但是否自动选中取决于 Codex 会话的 Skill 发现；自动选择不由 Python 单元测试保证。
