# R 4.5.0：QC 测试完成后 Rscript 返回 139

> HISTORICAL：本文记录早期诊断。当前所有 R 入口已统一为
> [`runtime/r45/run_r45.py`](../runtime/r45/run_r45.py)，环境定义仅在
> `runtime/r45/config.py`。本文提及的旧安装试验不属于当前工作流。

记录日期：2026-09-27。此文档记录本机排查及仓库内的兼容处理，不代表 Windows/R/`cli` 的通用修复。

## 现象与环境

- 指定 Rscript：`D:/R/R-4.5.0/bin/Rscript.exe`，R 4.5.0 UCRT。
- `.libPaths()`：`D:/R/R-4.5.0/library`；Seurat 5.4.0，DoubletFinder 2.0.6。
- DoubletFinder 的真实计算已打印 `PASS`，但原本的 Rscript 在脚本结束时返回 139。不能把“已打印 PASS”当作进程成功退出，也不能据此认定 DoubletFinder 缺失。

用指定 Rscript 执行的最小复现如下。`library()` 调用成功且 `PASS` 已打印；139 发生在脚本完成后的退出清理阶段，而不是加载包时。

| `--vanilla -e` 内容 | 观察到的退出码 |
|---|---:|
| `1+1` | 0 |
| `library(ggplot2); print('PASS')` | 139 |
| `library(SeuratObject); print('PASS')` | 139 |
| `library(Seurat); print('PASS')` | 139 |
| `library(DoubletFinder); print('PASS')` | 139 |
| `library(cli); print('PASS')` | 139 |

上述结果说明问题不限于 QC 分析代码或 DoubletFinder。`cli` 会作为这些包的依赖加载；单独 `loadNamespace("cli")` 后调用 `unloadNamespace("cli")` 也会在卸载过程中崩溃。本机 `cli` 3.6.6 的包构建版本提示为 R 4.5.3，但这不是已证实的根因。

## 排查过程

1. 检查了 R 版本、`.libPaths()`、相关包版本及 `R_HOME`、`R_LIBS`、`R_LIBS_USER`、`PATH`。没有因为系统 `PATH` 中其他 R 环境缺包而改用其他 R。
2. 试过进程内清理 `PATH` 和调整 locale，退出码仍为 139。
3. 仅在项目忽略的 `data/` 诊断目录中，用 R 4.5.0 分别从源码隔离安装 `cli` 3.6.6、3.6.5、3.6.4；三个版本仍复现 139。没有覆盖或删除 `D:/R/R-4.5.0/library` 中的包。
4. 检查 `cli` 3.6.6 源码的 `src/thread.c`：`CLI_NO_THREAD` 阻止创建计时线程；Windows 退出清理函数在检测到 `PROCESSOR_ARCHITECTURE=ARM64` 时跳过线程取消，否则可能执行 `pthread_cancel`。仅设置 `CLI_NO_THREAD=1` 仍返回 139；在退出清理前切换该分支则返回 0。因此最可能的问题位于本机 `cli` 原生线程清理路径，但尚未证明底层究竟是 `cli`、Rtools 运行时还是其他 DLL 的缺陷。

## 仓库内处理

[统一启动器](../runtime/r45/run_r45.py)设置项目内 `R_PROFILE_USER` 和进程局部的 `CLI_NO_THREAD=1`。[退出 profile](../runtime/r45/r45_profile.R)只在退出清理时将当前进程的 `PROCESSOR_ARCHITECTURE` 临时设为 `ARM64`；分析计算期间保留真实架构。未处理异常仍以状态码 1 退出。启动时缺少架构环境变量会按 `R.version$arch` 补齐真实架构，避免退出处理递归报错。

调用方式（Git Bash，在仓库根目录）：

```bash
bash qc/r450_rscript.sh qc/qc_precheck.R D:/CodexProjects/single-cell-skill GSE123456
bash qc/r450_rscript.sh qc/run_qc.R D:/CodexProjects/single-cell-skill GSE123456
```

启动脚本拒绝 `--vanilla` 和 `--no-init-file`，因为它们会禁用必需的 profile。若指定 Rscript 不存在则停止；不会自动切换 R、安装包或修改 QC 分析流程。

当前处理覆盖本项目 Loader、QC、测试和临时 R 命令。主动在分析中卸载 `cli` 不在退出保护范围内；若进程仍崩溃，记录 `R_NATIVE_CRASH` 并停止，保留已完成 artifact。此 profile 不应设置为机器全局配置。

## 验证

- `bash tests/test_r450_launcher.sh`：正常退出 0、故意触发 R 错误退出 1、`--vanilla` 被拒绝，均通过。
- `bash qc/r450_rscript.sh tests/test_qc.R .`：synthetic QC、用户决策和 clean output 测试通过，进程退出 0。
- `bash qc/r450_rscript.sh tests/test_qc_real_df.R .`：安装的 DoubletFinder 真正执行 `paramSweep`、pK 选择和两次调用，通过且退出 0；这不是对真实 GSE 的 QC 运行。
- Python 测试：69 项通过，1 项跳过。

关联实现提交：`0be8af8082304343680a6b0c0789c38e8efd6e98`。
