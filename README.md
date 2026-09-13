# GEO single-cell Seurat loader skill

Windows 11 项目级 Codex Skill：搜索 GEO、检查样本来源、等待用户确认分组、下载 processed counts、创建并验证 `seurat_raw.rds`。不执行正式单细胞下游分析。

Skill 位置：`.agents/skills/geo-single-cell-loader/`。在 Codex 中打开整个 `single-cell-skill` 文件夹作为 workspace/repository root；不要把 Skill 放进全局目录或 `.codex/skills`。脚本默认从自身位置解析 repository root，支持含空格路径，不包含电脑用户名或固定盘符。

## 工作目录

部署时用 PowerShell 读取真实桌面：

```powershell
$Desktop = [Environment]::GetFolderPath('Desktop')
if (-not $Desktop) { throw '无法解析当前用户真实桌面；请在用户本机会话重试。' }
$Project = Join-Path $Desktop 'single-cell-skill'
Set-Location $Project
```

```text
single-cell-skill/
├── .agents/skills/geo-single-cell-loader/
│   ├── SKILL.md
│   ├── agents/openai.yaml
│   ├── scripts/
│   │   ├── inspect_geo.py
│   │   ├── build_manifest.py
│   │   ├── download_processed.py
│   │   ├── build_seurat.R
│   │   ├── validate_seurat.R
│   │   ├── common.py
│   │   ├── seurat_common.R
│   │   └── check_dependencies.R
│   └── references/
├── datasets/       # 下载数据与已验证 cell map，绝不放进 Skill
├── manifests/      # 样本报告、分组确认、正式 manifest
├── output/         # 每个 GSE 的 RDS 和 run_summary.txt
├── logs/           # 下载日期、SHA256、运行记录
├── tests/          # 合成 fixture 生成器和离线测试
├── README.md
└── .gitignore
```

## 两阶段工作流

**阶段 A：先搜索、后报告。** 优先调用可用的 bio-server GEO MCP 工具；不可用时运行本项目的 HTTPS SOFT 查询。读取 Series、GSM、相关论文以及公开 supplementary sample/clinical metadata。脚本保留原始事实和来源，Codex 负责依据公开证据规范化组织、疾病和材料类型；空白和歧义不能靠猜测填充。

返回所有样本的 `sample、author_sample、tissue、disease、source_type、sample_description`，说明测序类型、processed 格式、raw/filtered 与独立 metadata 的存在情况。**group 尚未创建，请确认分组方式。** 不下载整套表达矩阵。开发用 `--max-samples` 的部分报告不能伪装成完整样本报告。

**阶段 B：用户明确决定 group。** 保存用户原话及精确 sample→group 映射，构建正式 manifest 和内容校验回执后才允许批量下载。用户改变任何映射，需重新确认。receipt 只是可追溯的确认记录，不是身份认证或证据真实性证明。

输入示例：

> 分析 GSE231993，先告诉我所有样本的来源和疾病状态。

确认示例（只在实际样本报告支持这些类别时使用）：

> HC样本group设为HC，UC-self control设为UC_control，UC炎症组织设为UC。继续下载并创建Seurat。

另一种输入：

> 下载并读取 GSE229413

这仍然必须先返回 sample report，不能立即开始大型下载。示例不构成任何真实数据的预设分组。

## metadata 规则

| 强制字段 | 含义 |
|---|---|
| database | 同一对象内一致的 GSE 来源 |
| sample | GSM 优先；或有公开映射的稳定测序样本 ID |
| tissue | 解剖组织；不能写 Tumor/Normal |
| disease | 供体疾病背景；PDAC uninvolved tissue 仍是 PDAC |
| source_type | Tissue / PBMC / Whole_Blood / Organoid / Cultured_TIL 等材料类型 |
| group | 仅由用户确认的分析分组 |

`patient/specimen/cohort/treatment` 只在资料可靠时添加，并有逐样本 evidence；不能猜 P01/donor3/T01 是患者 ID。部分缺失用 NA，并在 summary 报告。辅助 `author_sample/sample_description` 保留在 report/manifest。

manifest 是唯一映射依据。混合样本矩阵还需 manifest 显式引用的 `cell_map_path`，其中 `cell,sample` 必须与整个输入矩阵和相关 manifest 样本精确对应。加载时绝不从文件名重新推断 biological metadata。barcode 在合并前加 sample 前缀；缺失、重复、错配、合并后样本消失均报错。metadata 行顺序必须和 Cells 完全一致。

## 支持格式

- 标准 10x 三联文件，gzip/plain/legacy genes.tsv；规范化临时文件后用 Read10X。
- 10x HDF5：验证 HDF5 schema 后 Read10X_h5；同一样本默认 filtered，多个 filtered 则停止消歧。
- H5AD / H5AD.gz：zellkonverter 原生 R reader → SingleCellExperiment → 显式选择的 counts assay → Seurat；支持显式 raw alternative experiment。绝不使用 Read10X_h5，不静默使用 normalized X。
- TXT/CSV/TSV count matrix（可 gzip）：先确认方向、delimiter、ID 列、非表达列，再读取。拒绝无法解释的非数值列。
- Seurat RDS：检查对象类型及 counts layer，重建 counts-only Seurat；其他 RDS 类型明确报错。
- tar/zip：由 manifest 指定精确文件成员，安全地按成员释放；不递归盲解压。

具体 schema、下载 JSON 字段、pooled 数据和 counts provenance 见 Skill 的 [metadata-schema](.agents/skills/geo-single-cell-loader/references/metadata-schema.md) 和 [format-routing](.agents/skills/geo-single-cell-loader/references/format-routing.md)。

## 依赖与命令行

Python 3.10+；GEO、manifest 和下载脚本仅用标准库。开发 fixture 额外需要 `numpy、h5py`；H5AD fallback 需要 `anndata、numpy、scipy`。先定位实际解释器；不要使用 Windows Store 的 python/python3 占位程序。若使用 fallback，显式设置 `$env:GEO_SINGLE_CELL_PYTHON = $PythonExe`。

R、Seurat、SeuratObject、Matrix、jsonlite 是核心依赖；hdf5r 用于 H5；zellkonverter、SingleCellExperiment、SummarizedExperiment 用于 H5AD。检查同时列出常用 dplyr、data.table 版本，本实现不依赖它们进行 metadata join。

```powershell
# 用实际安装路径设置，不要假定 PATH 中的 R 就是装有 Seurat 的版本。
$PythonExe = (Get-Command py -ErrorAction Stop).Source  # 或设置已核实的 Python 可执行文件
$RscriptExe = (Get-Command Rscript -ErrorAction Stop).Source
$Scripts = '.agents/skills/geo-single-cell-loader/scripts'
& $RscriptExe "$Scripts/check_dependencies.R"
& $PythonExe "$Scripts/inspect_geo.py" GSE231993
# 先在 Codex 中完成证据审阅、样本报告和用户明确分组。
& $PythonExe "$Scripts/build_manifest.py" --report manifests/GSE231993_sample_report.csv --confirmation manifests/GSE231993_group_confirmation.json --output manifests/GSE231993_sample_manifest.csv
& $PythonExe "$Scripts/build_manifest.py" --validate manifests/GSE231993_sample_manifest.csv
& $PythonExe "$Scripts/download_processed.py" manifests/GSE231993_sample_manifest.csv
& $RscriptExe "$Scripts/build_seurat.R" . manifests/GSE231993_sample_manifest.csv
& $RscriptExe "$Scripts/validate_seurat.R" . manifests/GSE231993_sample_manifest.csv output/GSE231993/seurat_raw.rds
```

在选定的 R 环境安装缺包：

```r
install.packages(c("Seurat", "dplyr", "data.table", "Matrix", "jsonlite", "hdf5r", "BiocManager"))
BiocManager::install(c("zellkonverter", "SingleCellExperiment"), ask=FALSE, update=FALSE)
```

Windows 源码安装需要匹配 R 的 Rtools，并把其 `usr/bin` 加到**当前安装进程** PATH；无需修改系统设置。缺少 zellkonverter 时可在明确的 Python 环境安装 `anndata numpy scipy` 并设置 `GEO_SINGLE_CELL_PYTHON`。缺少 H5AD 读取依赖只阻止 H5AD 路线，不影响项目文件创建或其他格式。原生 R reader 出现 warning 或 fallback 发现 normalized matrix 时停止。每次运行记录实际 reader 和包版本。

在某些精简的 Windows 非交互进程里，标准 `PROCESSOR_ARCHITECTURE` 变量可能缺失。`cli` 3.6.6 的卸载代码会因此触发访问冲突；本项目只在该变量缺失时，根据 R 架构为当前 R 进程补齐它，不修改系统环境。见 [r-lib/cli issue #375](https://github.com/r-lib/cli/issues/375) 和候选修复 [PR #838](https://github.com/r-lib/cli/pull/838)。

## 输出

`output/<GSE>/seurat_raw.rds` 和 `run_summary.txt`；`manifests/<GSE>_sample_manifest.csv`、confirmation 回执；`logs/<GSE>_download.json`。summary 包含日期、格式、读取方法/count source、细胞/基因/样本数量、sample/group 计数、组织/疾病/材料类型、可选字段、输入输出、warnings、sessionInfo。已有正式 RDS 不会被静默覆盖。

## 测试

```powershell
& $PythonExe -m unittest discover -s tests -p 'test_*.py' -v
& $PythonExe tests/make_fixtures.py
& $RscriptExe tests/test_seurat.R .
& $RscriptExe "$Scripts/build_seurat.R" . manifests/GSE999999999_sample_manifest.csv
& $RscriptExe "$Scripts/validate_seurat.R" . manifests/GSE999999999_sample_manifest.csv output/GSE999999999/seurat_raw.rds
```

合成数据用保留测试编号 `GSE999999999`，不查询 GEO，不对应真实生物学组。生成器只重写该测试 manifest；重复运行端到端之前需明确清理/归档其旧输出。测试实际读取 10x/H5/H5AD/text/RDS，覆盖错误分组、排序、normalized X、barcode 和错误映射。

真实测试类型：[GSE231993](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE231993)、[GSE202051](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202051)、[GSE211644](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE211644)、[GSE229413](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE229413)。实际执行范围和结果见 [test-datasets.md](.agents/skills/geo-single-cell-loader/references/test-datasets.md) 与 [validation-report.md](tests/validation-report.md)。大型真实数据不重复下载，真实 group 仍须用户决定。

## V1 边界

不跑 Cell Ranger/FASTQ/SRA 重新定量，不做额外 QC、NormalizeData、SCTransform、变异基因、ScaleData、PCA/UMAP、邻居图、聚类、DoubletFinder、注释、差异或富集分析。按要求使用 CreateSeuratObject(min.cells=3, min.features=200)，这些创建阈值会删除部分细胞/基因并记录，除此以外不做过滤。

不保证任意作者私有格式、异常 H5AD、混合 assay 或超内存数据能直接读取；歧义时停止。公开元数据缺失时，需要补充证据才能继续。Stage A 不是完全无人监督的临床信息推断器。Git 忽略 datasets/output/logs/manifests、RDS/H5/H5AD/FASTQ/SRA 和凭据；只上传代码、文档及小型测试生成器。
