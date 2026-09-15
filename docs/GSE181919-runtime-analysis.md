# GSE181919 运行过程与慢速原因分析

日期：2026-09-15（Asia/Shanghai）

## 结论

本次耗时的主要原因不是网络下载，而是 GEO 将稀疏的单细胞 UMI counts 发布为超宽的、以零为主的制表符文本矩阵。压缩文件只有 127,878,601 bytes，但解压后为 2,179,824,249 bytes，矩阵包含 20,000 genes × 54,239 cells，即 1,084,780,000 个数值字段。`base::read.table` 必须先将这些字段解析成密集的 R data frame/matrix，再转换为稀疏矩阵，实测私有内存峰值约 17–21 GB。

原构建器还有一个只在 pooled matrix 场景暴露的性能问题：manifest 每个样本一行，而 37 行都指向同一个 UMI 文件；旧循环会调用 `read_counts()` 37 次，相当于计划重复解析同一个 2.18 GB 文本 37 次。本次运行在观察到第二次读取开始后主动中止，并增加了按 `local_path` 缓存；修复后物理矩阵只读取一次。

## 输入数据

- GEO accession: [GSE181919](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181919)
- 研究：头颈鳞癌从邻近非肿瘤组织、白斑、原发癌到淋巴结转移的单细胞转录组
- 样本：37 个组织样本，23 位患者
- 用户确认的 group：NL 9、LP 4、CA 20、LN 4
- 表达文件：`GSE181919_UMI_counts.txt.gz`
- 条形码文件：`GSE181919_Barcode_metadata.txt.gz`
- 表达文件 SHA-256：`cc5e83a48afa4ec48f2b2f64d449354124cb54976c7a9755f96a759543e5c6d6`
- 表达文件压缩大小：127,878,601 bytes
- 解压大小：2,179,824,249 bytes
- 矩阵维度：20,000 genes × 54,239 cells

运行数据保存在本地 `data/GSE181919/`，由 `.gitignore` 排除；GitHub 仅保存代码、测试和本复盘，不保存原始矩阵或 266 MiB RDS。

## 运行过程

1. 通过全局 GEO MCP 获取完整 Series、37 个 GSM 和两个 Series-level supplementary files。
2. 下载 339,614-byte 条形码元数据并检查内容。公开文件表头只有 8 个标签，但每个数据行有 9 列；缺失的是最前面的 barcode 标签。
3. 对 UMI gzip 文件做 2 MiB header probe。矩阵表头有 54,239 个 cell ID，数据行有 54,240 列，因此第一列是没有标题的 gene ID。
4. 将条形码元数据末尾的 `-<library>` 转为矩阵使用的 `.<library>`；54,239 个转换后条形码与 UMI 表头完全同序匹配，无缺失、无额外、无重复。
5. 根据用户确认建立 NL/LP/CA/LN group map、37-row manifest 和显式 `cell_map.csv`，再下载完整表达文件并记录 URL、UTC 时间、bytes 和 SHA-256。
6. 完整 gzip 流校验通过；20,000 行 feature ID 均唯一，每个表达行宽度一致。
7. 第一次 R 构建读取了整张宽表，随后因 R 将无标题 feature 字段暴露为非数值 `row.names` 列而安全停止；没有生成最终 RDS。
8. 增加显式 `feature_column=__row_names__` 路由。第二次构建成功完成一次文本到稀疏矩阵的转换，但随后开始为下一 manifest 行重复读取同一 pooled matrix；在未生成 RDS 前主动中止。
9. 增加共享矩阵缓存和共享 reader 配置一致性检查。最终运行只读取一次物理矩阵，从缓存中按 cell map 提取 37 个样本，加 GSM 前缀后在 counts 层合并。
10. 全数据统一调用一次 `CreateSeuratObject(min.cells=3, min.features=200)`。内存中验证、pending RDS 保存、序列化回读验证和独立 `validate_seurat.R` 均通过。

## 实测阶段与瓶颈

| 阶段 | 观察结果 | 是否主要瓶颈 |
|---|---|---|
| 下载 121.95 MiB gzip | 约几十秒；SHA-256 与字节数记录成功 | 否 |
| gzip 完整解压扫描 | 数秒；2.18 GB、20,001 行（含表头） | 否 |
| `read.table` 解析 10.85 亿字段 | 数分钟；内存逐步升至约 17–21 GB | 是 |
| 密集矩阵转稀疏矩阵 | 内存峰值后明显回落 | 是 |
| 旧版按 37 manifest rows 重读同一文件 | 理论上将最慢阶段放大约 37 倍；第二次读取开始后被中止 | 严重缺陷，已修复 |
| 按 cell map 拆分和重新合并 | 修复缓存后在单次解析结果上完成 | 次要 |
| RDS 保存、回读和独立验证 | 额外几十秒到约一分钟，属于必要完整性检查 | 否 |

上述时间为交互运行期间的近似观察值，不是基准测试结果。内存值来自 Windows 进程的 Private Memory/Working Set 采样。

## 本次代码修复

### 1. 缺失 feature 表头

文本 reader 新增显式 `feature_column=__row_names__` 路由。只有在人工核验确认“header 比数据行少一列”时才能使用；reader 同时处理 R 自动创建的 `row.names` 列和实际 row names，普通带标签的 ID 列仍走原逻辑。

### 2. pooled matrix 只读一次

`build_seurat.R` 现在按 `local_path` 缓存 `read_counts()` 结果。同一路径被多个 manifest 样本共享时，先验证 `file_type`、`count_source`、方向、分隔符、feature column 等 reader 设置完全一致，再从缓存按显式 cell map 取出各样本。

这项修复消除了 37 次重复文本解析；最终构建实际只读取一个物理 source matrix。

### 3. 回归测试

R reader 测试增加两种缺失 feature 表头布局：

- R 将 feature ID 放进实际 row names；
- R 将其暴露为合成 `row.names` 列。

两者均要求输出保持正确的 gene/cell 名称与矩阵维度。

## 最终结果

- Seurat object：20,000 genes × 54,239 cells
- 样本：37，全部非空
- CA：23,088 cells
- LN：8,204 cells
- LP：6,527 cells
- NL：16,420 cells
- RDS 大小：279,066,168 bytes（266.14 MiB）
- RDS SHA-256：`d0ce5e6fcf4a750618cf750ae484a58e6e0f897d48ae2284e11f39267dc4f603`
- 构建策略：全部 counts 合并后统一应用 `min.cells=3, min.features=200`
- 未执行额外 QC、归一化、降维、聚类、注释或差异分析

## 仍可继续优化的部分

目前 pooled matrix 重复读取问题已经解决，但首个 `read.table` 仍需要密集解析。下一步最有价值的优化是增加流式 text-to-sparse reader：逐行读取文本，直接构造稀疏坐标或分块写入临时 Matrix Market/稀疏格式，避免同时保留大型 data frame、dense matrix 和 sparse matrix。该优化需要专门的性能与数值等价测试，本次没有为了追求速度而绕过 counts、feature、cell map 或序列化验证。
