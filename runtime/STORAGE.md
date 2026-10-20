# RDS 存储与恢复协议

`runtime/storage_config.json` 独立于科学参数，默认 `compact`。它不改变过滤、DoubletFinder、DecontX、回归、PC 选择、Harmony 或聚类算法。

```text
DATASET/
  01_download/seurat_raw.rds
  02_QC/03_ambient_RNA/seurat_QC-3.rds
  03_analysis/seurat_analysis.rds
  03_analysis/runs/ID/seurat_analysis.rds  # 用户采用的独立方案
  .work/storage_policy.json
  .work/checkpoints/
    registry.json
    qc/{seurat_QC-1.rds,seurat_QC-2.rds,state/,formal_validation.json}
    analysis/baseline/{seurat_pca.rds,seurat_unintegrated.rds,
      seurat_harmony.rds,seurat_cluster_review.rds,regression_review.rds,
      count/seurat_sensitivity.rds,state/,formal_validation.json}
    analysis/runs/ID/...
  storage_report.csv
  storage_summary.json
```

图、表、配置、日志和用户决定保留在原阶段目录。request、pending publication、runtime completion、恢复 staging 和 manifest 均位于 `.work/checkpoints`。R45 进程审计仍在 `.runtime/r45`；下载十文件审计不变。

## 模式与引用

QC/分析均支持 `--storage-mode compact|full`、可重复的 `--retain-checkpoint STAGE_OR_BASENAME`；设置保存在 `.work/storage_policy.json`，恢复时不必重复指定。显式保留列表替代旧列表。阶段名包括 `qc1`、`qc2`、`pca`、`harmony`、`cluster-review`、`regression-review`；`analysis/runs/ID:pca` 限定一个运行；`seurat_sensitivity.rds` 保留候选对象。full 保留全部中间对象；发布成功后仍释放 staging 硬链接别名。

```bash
/c/Python312/python.exe qc/run_workflow.py Study_A --storage-mode compact --retain-checkpoint qc2
/c/Python312/python.exe analysis/run_workflow.py Study_A --storage-mode full
```

用户提供 raw Seurat RDS 可用 `qc/run_workflow.py Study_A --step qc1 --input ABSOLUTE_RDS`。记录绝对路径、SHA256、大小及禁止删除标记；后续 QC 自动使用引用，缺失/变化即停止，不复制、不移动、不删除源文件。分析 `--input` 同样记录引用，独立方案有独立作用域。

## 生命周期

- `INCOMPLETE`：执行/验证失败，包括创建 RDS 后进程崩溃；保留恢复文件，禁止清理。
- `COMPLETE_RETAINED`：成功完成且对象保留，可以作为输入。
- `COMPLETE_CLEANED`：计算已完成，已有独立读回验证的正式下游对象；旧检查点只作为历史证据，不能当作可读输入。

QC-3 经独立进程校验后才可清理 QC-1/2。分析等待分辨率期间保留复核必需对象；最终对象保存、独立校验及 Marker/绘图 API 检查全部成功后，才清理同一次分析的临时对象。任何保存/验证错误保留检查点。full、指定保留列表和 INCOMPLETE 阻止对应清理。

登记表记录路径、大小、SHA256、阶段、作用域、正式/临时属性和生命周期。清理匹配同一作用域的正式对象，并重新核实其 SHA256。先将 manifest 的对象从 `outputs` 移入 `cleaned_outputs`、写入下游凭据和清理状态，再记录删除意图，最后 unlink；中断后可幂等继续，不无条件忽略丢失文件。数据集写入使用 OS 锁，进程退出自动释放。

清理后 `all` 验证正式对象、参数、引用、manifest 和下游凭据，不依赖已删除对象。显式请求已清理的分析阶段时，从这次分析记录的 QC 引用创建 `baseline-rebuild-...`/`ID-rebuild-...` 独立运行，重建必要前驱；原分析 `last_rebuild_result.json` 记录新位置。原正式结果不覆盖，也不从其他 run 搜索中间对象。QC-1/2 重建同样隔离至 `.work/checkpoints/qc/rebuilds/ID`，报告写入 `02_QC/rebuilds/ID`，继承 canonical 用户决定；`02_QC/last_rebuild_result.json` 给出 `--run-id ID` 恢复命令。原 QC-3 和原图表/报告不覆盖。

## 无损精简与验证

`seurat_storage.R` 在正式保存边界统计 assay/layer/reduction/graph 的实际 `object.size()`，写入 `object_storage_inventory.csv`。诊断值存在包含关系，不能直接相加视为 RDS 大小。compact 分析对象删除 `scale.data`、与 canonical counts/data 完全 identical 的其他 layer、与正式 `umap` 数值/细胞顺序完全相同的 UMAP 别名、候选 `cluster_res.*`/`analysis_snn_res.*` 字段、已完成用途的图/邻居缓存。

保留 RNA counts、decontX fractional counts、各 assay 标准化 data、PCA、Harmony、唯一 unintegrated UMAP、正式 UMAP、QC/周期/污染/样本/分组、正式 cluster 和参数/决定记录。完整分辨率分配仍保存为 `cluster_assignments.csv`；对象保留候选列名、分辨率及正式标签指纹。不删除未知 assay/其他 reduction，不改变矩阵类型、精度或数值。

精简前后逐项比较表达层指纹和受保护 metadata；`validate_formal.R` 通过独立 R45 进程读回检查 SHA256。分析对象还实际构建 FeaturePlot、对子集运行 FindMarkers、验证注释字段赋值。API 检查不产生生物学结论；不足两个可测试群时仅在副本上构造临时 API 测试分组。正常运行通过小型 JSON 决定策略，无须 Codex 加载大对象选择清理文件。

## 旧结果与空间报告

已有 QC/分析 RDS 且无新策略登记的数据集不会被接管；传 `--storage-mode compact` 也不授予旧文件删除权限。历史 `data/DATASET/qc` 同样受保护。旧结果实际清理需先展示清单/依赖并获得用户对具体文件的确认；本工具没有绕过确认的旧结果清理入口。

```bash
/c/Python312/python.exe runtime/storage.py D:/path/old_dataset --output-dir D:/path/inventory
```

盘点只向指定目录写小报告，不读取 RDS 内容，不修改源数据集。CSV 包含路径、大小、类型、阶段、正式/临时属性、状态、安全清理标记、待确认候选、依赖、前后总字节数和节省量；未知旧依赖明确标 UNKNOWN。候选不是已验证可删除文件。摘要 JSON 提供总量、状态计数、引用和报告路径。统计为逻辑文件大小。before/after 是最后一次完成清理事务的实测快照（含当时 metadata 更新）；current_total_bytes 是生成报告时的现存总量；saved_bytes 累计安全删除的检查点字节。cleanup_events 保留逐次实测快照。硬链接/文件系统压缩可能影响物理分配。报告本身不计入数据总量。

## 验证入口

`tests/test_storage.py` 覆盖检查点可用、QC3→分析、暂停恢复、失败保留、full/指定保留、清理中断、清理后启动、隔离重建、旧结果保护及输入引用。`tests/test_seurat_storage.R` 用真实 Seurat 5 测试 counts/标准化层/metadata 精确保真与绘图/Marker/注释 API。`tests/run_storage_acceptance.py --output NEW_DIRECTORY` 运行实际 R45 合成数据 PCA/Harmony/多分辨率/回归对照/选择/精简/清理/恢复/审计验收。模拟选择只在 SYNTHETIC 测试数据中使用。
