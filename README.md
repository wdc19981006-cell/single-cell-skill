# GEO single-cell loader Skill

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
6. 全部结果保存在 `data/<GSE>/`。

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

**Standard output：**无论来源是 10x、H5AD、RDS、text 或 pooled matrix，最终均写为 `data/<GSE>/seurat_raw.rds`：一个 clean Seurat 对象、一个 RNA raw counts layer、标准 sample/cell metadata；不继承作者 PCA、UMAP、SCT/integrated assay、graph、neighbor、cluster 或 normalized expression。

在样本资料确实支持这些类别时，用户可回复：

```text
健康人为 HC；UC患者非炎症组织为 UC_control；UC炎症组织为 UC。
```

这是确认方式示例，不是预设分组。group 只能由用户决定，不能根据 disease、文件名、Tumor/Normal 字样自动产生。PDAC 的 uninvolved 样本仍有 PDAC 疾病背景，不能因此自动归入 Tumor 或 Healthy。

## 一个 GSE，一个目录

```text
single-cell-skill/
├── .agents/skills/geo-single-cell-loader/
├── data/
│   ├── .gitkeep
│   └── GSE231993/
│       ├── raw/
│       ├── seurat_raw.rds
│       ├── sample_info.txt
│       ├── group_confirmation.json
│       └── cell_map/             # 仅在重建对象需要时
├── tests/
├── README.md
└── .gitignore
```

运行完成后关注这些本地文件：

| 内容 | 用途 |
|---|---|
| `raw/` | 从 GEO 下载且未经本 Skill 修改的源 expression 文件。这里的 raw 不一定表示测序 FASTQ raw reads；GEO 的 processed counts 也放在这里。10x 三联/H5 按样本存放，pooled H5AD 可以共用一个源文件。 |
| `seurat_raw.rds` | 仅从明确的 counts 构建并验证的 Seurat 对象。已有文件不会静默覆盖。 |
| `sample_info.txt` | UTF-8 纯文本说明，记录研究、样本、用户分组、运行状态及最终统计。 |
| `group_confirmation.json` | 用户确认的精确分组映射。 |
| `cell_map/` | 仅在 pooled 矩阵或独立映射文件重建对象时保留。 |

`.workflow/` 只在运行期间暂存发现、确认回执、manifest、下载来源和审计文件。每次真实 GSE 运行结束会将 10 个审计文件推送到独立私有仓库 [single-cell-skill-audit](https://github.com/wdc19981006-cell/single-cell-skill-audit) 的 `<GSE>/<run-id>/`；上传成功后清理本地临时文件。失败或中断也上传；上传失败则保留临时文件以便重试。表达矩阵、RDS、H5/H5AD 不进入审计仓库。格式转换使用系统 temporary directory 并清理。Git 忽略 `data/*`，只保留 `data/.gitkeep`；运行数据和 TXT 默认只保存在本地。用户明确要求分享时，可将数据包上传到本仓库的 Releases，数据不进入 Git 历史。

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

确认分组后状态为 `GROUP_CONFIRMED`，逐样本增加 group，并增加 GROUP SUMMARY。创建、序列化并验证成功后状态为 `COMPLETE`，附加 PROCESSING SUMMARY：总样本/细胞/基因数、sample/group 细胞数、输出位置、实际读取格式和 reader、counts source、Seurat 版本、创建时间和 warnings。

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

```powershell
# 将变量设置为本机已验证的 Python 与 Rscript 可执行文件绝对路径。
$Scripts = '.agents/skills/geo-single-cell-loader/scripts'
$Workflow = 'data/GSE231993/.workflow'
& $PythonExe "$Scripts/audit_run.py" start GSE231993
& $RscriptExe "$Scripts/check_dependencies.R"
# 仅在真实 geo MCP 调用失败时使用，并记录实际原因：
& $PythonExe "$Scripts/inspect_geo.py" GSE231993 --fallback-reason "geo/get_geo_info tool error: <actual error>"
# Skill 审阅证据、完善 sample_report.csv 与 inspection.json 后：
& $PythonExe "$Scripts/sample_info.py" GSE231993
# 用户确认 group，并保存 group_confirmation.json 后：
& $PythonExe "$Scripts/build_manifest.py" --report "$Workflow/sample_report.csv" --confirmation "$Workflow/group_confirmation.json" --output "$Workflow/sample_manifest.csv"
& $PythonExe "$Scripts/build_manifest.py" --validate "$Workflow/sample_manifest.csv"
$env:GEO_RSCRIPT = $RscriptExe
& $PythonExe "$Scripts/run_confirmed.py" GSE231993
```

测试：

```powershell
& $PythonExe -m unittest discover -s tests -p 'test_*.py' -v
& $PythonExe tests/make_fixtures.py
$env:GEO_SINGLE_CELL_PYTHON = $PythonExe
& $RscriptExe tests/test_seurat.R .
& $RscriptExe tests/test_global_min_cells.R .
& $RscriptExe tests/test_pooled_inputs.R .
& $PythonExe tests/run_end_to_end.py --rscript $RscriptExe
```

合成数据使用保留测试编号 `GSE999999999`，分组批准明确标注 SIMULATED；不对应真实临床组。重复端到端测试前显式归档已有测试 RDS，构建脚本不会覆盖。旧布局数据需按 GSE 迁移并比较大小/hash；保留原始 manifest/receipt，不自动修改 hash 来认可新路径。无法确定归属的文件保留并报告。

详见 [metadata contract](.agents/skills/geo-single-cell-loader/references/metadata-schema.md)、[format routing](.agents/skills/geo-single-cell-loader/references/format-routing.md) 和 [实际验证报告](tests/validation-report.md)。公开元数据仍需证据审阅，任意作者私有格式或超内存矩阵可能需要额外适配。Skill 已启用隐式调用，但是否自动选中取决于 Codex 会话的 Skill 发现；自动选择不由 Python 单元测试保证。
