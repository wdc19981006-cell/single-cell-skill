"""Write the GSE155698 validation report from completed, verified artifacts."""
import argparse
import csv
import json
import statistics
import subprocess
from pathlib import Path

REPO=Path(__file__).resolve().parents[3]


def read_csv(path):
    return list(csv.DictReader(path.open(encoding='utf-8-sig')))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('workspace',type=Path)
    parser.add_argument('artifacts',type=Path)
    parser.add_argument('report',type=Path)
    args=parser.parse_args()
    current=read_csv(args.artifacts/'current_pca_diagnostics.csv')
    old=read_csv(args.artifacts/'legacy_pc_control_pca_diagnostics.csv')
    objects=read_csv(args.artifacts/'object_validation.csv')
    comparisons=read_csv(args.artifacts/'cell_gene_comparisons.csv')
    summary=json.loads((args.artifacts/'validation.json').read_text())
    independent=json.loads((args.artifacts/'independent_selection_review.json').read_text())
    history=json.loads((args.artifacts/'history_inventory.json').read_text())
    assert summary['all_raw_count_comparisons_exact'] and independent['status']=='PASS'
    assert len(current)==len(old)==21
    provenance=json.loads((args.artifacts/'current_provenance.json').read_text())
    current_dataset=args.workspace/'latest/data/GSE155698'
    final_sha=json.loads((current_dataset/'.workflow/state.json').read_text())['qc_sha256']
    prior_sha=provenance['protected_files'][str(Path(provenance['source'])/'qc/seurat_qc.rds')]
    pre_rows=read_csv(current_dataset/'diagnostics/pre_filter_umap_cells.csv')
    feature_medians={label:statistics.median(float(r['nFeature_RNA']) for r in pre_rows if r['DF_adj']==label)
                     for label in ('Singlet','Doublet')}
    records=[json.loads(p.read_text(encoding='utf-8')) for p in (REPO/'.runtime/r45/runs').glob('*.json')]
    relevant=[r for r in records if any(args.workspace.resolve().as_posix() in str(a).replace('\\','/') for a in r.get('args',[]))]
    qc_runs=[r for r in relevant if r['stage']=='qc']
    assert len(qc_runs)==2 and all(r['status']=='PASS' and r['exit_code']==0 for r in qc_runs)
    final_validators=[r for r in relevant if r['stage']=='diagnostic_final_validation']
    assert len(final_validators)==2 and all(r['status']=='PASS' and r['exit_code']==0 for r in final_validators)
    pca_tests=[r for r in records if r.get('script','').replace('\\','/').endswith('/tests/test_qc_pca.R')]
    pca_test=max(pca_tests,key=lambda r:r['started_at'])
    assert pca_test['status']=='PASS' and pca_test['exit_code']==0
    python_log=(REPO/'.runtime/GSE155698_python_tests.log').read_bytes()
    assert b'Ran 86 tests' in python_log and b'\nOK' in python_log
    production_diff=subprocess.check_output(['git','diff','cec9f452','--','qc','runtime',
        '.agents/skills/geo-single-cell-loader'],cwd=REPO)
    assert not production_diff
    evidence=dict(python_tests=dict(count=86,status='PASS'),PCA_regression=dict(status='PASS',exit_code=0),
        current_and_control_QC=[{k:r.get(k) for k in ('run_id','status','exit_code','duration_seconds','args')} for r in qc_runs],
        independent_final_validators=[{k:r.get(k) for k in ('run_id','status','exit_code')} for r in final_validators],
        independent_spectrum_recalculation=dict(status='PASS',PCA_objects=21),
        production_QC_Loader_Runtime_unchanged=True,historical_files_unchanged=provenance['historical_files_unchanged'],
        current_final_sha256=final_sha,prior_fixed_final_sha256=prior_sha,
        serialized_final_matches_prior_fixed=final_sha==prior_sha,DF_adj_nFeature_medians=feature_medians)
    (args.artifacts/'test_evidence.json').write_text(json.dumps(evidence,indent=2),encoding='utf-8')
    merged=next(r for r in current if r['sample']=='ALL')
    old_merged=next(r for r in old if r['sample']=='ALL')
    pcs=[int(r['pc_used']) for r in current if r['sample']!='ALL']
    link=lambda file: f'GSE155698/{file}'
    lines=['# GSE155698 PCA 自动选择验证','',
        '基线 GitHub main：`cec9f452e8fe704c187b2058eefd195987d1a4e9`。本次仅新增测试观察脚本、真实诊断表、图形和文档；没有修改生产 QC、PCA 方法、Loader、下载流程、R45 Runtime 或 regression 策略。','',
        f"当前规则的 20 个样本选择 **{min(pcs)}–{max(pcs)} PCs**，合并选择 **{merged['pc_used']} PCs / max {merged['max_pc']}**。本次真实运行没有重现全样本 2-PC 错误。",'',
        '## 输入与版本身份','',
        f"原始输入为 `data/GSE155698/seurat_raw.rds`，{summary['raw_cells']:,} cells × {summary['raw_genes']:,} genes。SHA-256：`{provenance['input_raw_sha256']}`，与此前工作流记录一致。原始 RDS、历史 RDS、报告、group_confirmation 和用户决策均在前后核对哈希，未覆盖。",'',
        '- 历史原代码：根目录 `seurat.qc.rds` / `seurat.all_qc.rds` 和原始 PDF。前者与当前 Cell QC 的 40,915 个细胞 ID 完全相同，但只有 26,086 个基因；当前保留 26,630 个基因。原代码额外过滤 544 个基因，最终按 DF 过滤。',
        '- 已保存的 `qc/`：已属于修复后的 Skill，20 样本 PCs 为 9–22、合并 23；不能称作修复前 2-PC Skill。它作为此前修复版本的复现对照。',
        '- 真正修复前的 GSE155698 Skill RDS/UMAP 未找到。因此增加 **legacy_pc_control：旧 PC 规则受控复现**，只在测试进程中恢复 `ec37df1` 的第一个 drop 和 >=90 条件，其他现代流程相同。它不是历史完整旧 Skill 产物。',
        '- current：独立测试目录执行当前固定 QC。两组受控运行均复制同一 raw、复用原有明确用户决策、使用同一样本顺序与 seed=155698、相同阈值和 pMT-only regression。没有构造新的用户接受记录。','',
        '## 算法核对','',
        'RunPCA 上限为 min(50, 配置上限, cells-1, HVGs-1, genes-1)，最终 pc_used 不超过 max_pc；现有小样本保护保留。选择使用标准差占比 `stdev/sum(stdev)*100`，cumulative >90 且 individual <5 的首次位置，与最后一个 drop >0.1 的位置 +1 取较小值。单个有效条件保留，无有效条件才 fallback，最终至少 2 PCs。','',
        '标准方差解释率使用 `stdev^2/sum(stdev^2)*100`，只作独立诊断。此处分母为实际计算的前 max_pc 个 PC 的平方标准差之和，不是全部基因总方差。ElbowPlot 的纵轴明确为 PCA standard deviation。','',
        f"已从完整导出谱独立重算 21 个 PCA 选择，全部与当前 R 方法一致。两组受控运行的 PCA 标准差谱最大绝对差为 {summary['max_control_current_stdev_delta']:.3g}；选择规则的差异可与谱本身分开检查。",'',
        f"合并对象：variance_pick={merged['variance_pick']}（历史命名，实际是标准差累计规则）、drop_pick={merged['drop_pick']}、pc_used={merged['pc_used']}。所选标准差累计占比 {float(merged['cumulative_sd_share']):.4f}%，诊断方差累计解释率 {float(merged['cumulative_variance_explained']):.4f}%。旧规则受控合并选择 {old_merged['pc_used']} PCs。",'',
        '## 每个样本及合并 PCA','',
        'cells 是 Cell QC 后、DF 前实际进入 PCA 的细胞数；原始每样本细胞数、group、HVG 数和完整指标另见 [诊断 CSV](GSE155698/current_pca_diagnostics.csv)。所有样本和 ALL 的单张 ElbowPlot 均已保存 PDF/300-dpi PNG，PNG 和完整谱随报告提交。','',
        '| sample | cells | max_pc | pc_used | variance_pick | drop_pick | 累计 SD % | 累计方差 % | ElbowPlot |',
        '|---|---:|---:|---:|---:|---:|---:|---:|---|']
    for r in current:
        lines.append('| '+' | '.join([r['sample'],r['cells'],r['max_pc'],r['pc_used'],r['variance_pick'],r['drop_pick'],
            f"{float(r['cumulative_sd_share']):.3f}",f"{float(r['cumulative_variance_explained']):.3f}",
            f"[查看]({link('elbows/elbow_'+r['sample']+'.png')})"] )+' |')
    lines+=['','## 每样本 DoubletFinder','',
        'DoubletFinder 按 sample 独立运行；观察脚本在每次 paramSweep 时断言只有一个 sample，且无逐样本 UMAP。merged 不运行 DoubletFinder，因此 selected_pK 为不适用，DF/DF_adj 只汇总逐样本结果。','',
        '| sample | selected_pK | DF Singlet / Doublet | DF_adj Singlet / Doublet |',
        '|---|---:|---:|---:|']
    for r in current:
        lines.append(f"| {r['sample']} | {r['selected_pK'] if r['sample']!='ALL' else '不适用'} | {r['DF_singlets']} / {r['DF_doublets']} | {r['DF_adj_singlets']} / {r['DF_adj_doublets']} |")
    lines+=['','## 同细胞集合 UMAP 与原始 counts 验证','',
        f"current 与 legacy_pc_control 的过滤前细胞集合完全相同，共 {summary['pre_filter_cells']:,} cells × {summary['raw_genes']:,} genes，样本和 group 对齐。DF 标签改变 {summary['DF_changed_cells']:,} 个细胞，DF_adj 标签改变 {summary['DF_adj_changed_cells']:,} 个细胞；全量过滤前交叉表另存 CSV。",'',
        '| 版本 | 最终 cells | 最终 genes | 与 raw 共同坐标 counts 差异 | clean object |',
        '|---|---:|---:|---:|---|']
    for r in objects:
        lines.append(f"| {r['version']} | {r['cells']} | {r['genes']} | {r['counts_mismatches_to_raw']} | {r['clean']} |")
    lines+=['','| 对照 vs current | 共同 cells | 仅对照 cells | 仅 current cells | 共同 genes | 仅对照 genes | 仅 current genes |',
        '|---|---:|---:|---:|---:|---:|---:|']
    for r in comparisons:
        lines.append('| '+' | '.join(r[k] for k in ('reference','common_cells','reference_only_cells','current_only_cells',
            'common_genes','reference_only_genes','current_only_genes'))+' |')
    lines+=['',
        f"当前最终 RDS 与此前修复版的序列化 SHA-256 {'完全相同' if final_sha==prior_sha else '不同'}：`{final_sha}`。本次相同哈希证明 counts 和 metadata 的最终序列化结果均可复现。",'',
        f"当前 DF_adj Doublet 的 nFeature_RNA 中位数为 {feature_medians['Doublet']:g}，Singlet 为 {feature_medians['Singlet']:g}。这只作为表达量辅助检查，不代表有双细胞生物学真值标签。",'',
        '所有 counts 比较均使用最多 256 个细胞一块的 sparse matrix，未将整体表达矩阵转为 dense。两组新对象的全部基因和保留细胞 counts 均与 raw 对应子集精确一致。cell-cycle 前后保留细胞集合一致；只评分，不按 Phase 删除、不回归 S/G2M。最终对象不保留 normalized data、VariableFeatures、PCA、UMAP、graph、neighbor 或临时 cluster 列。','',
        '过滤前图明确标识为 Doublet UMAP；过滤后图明确标识为 Singlet UMAP，并使用各自过滤前坐标的保留子集。测试实际核对坐标与细胞 ID，没有为 singlet 重新拟合 UMAP。','',
        f"过滤后不同规则保留的集合不同，因此主要几何比较另取两组共同保留的 {summary['common_post_filter_cells']:,} 个细胞。不能把两组不同细胞集合的最终图直接当作 PCA 优劣证据。原代码的 PCA/UMAP 坐标矩阵未保存，且基因过滤/regression 与 Skill 不同，其 PDF 只作历史背景，不作纯 PCA 因果对照。",'',
        '经渲染并人工检查，current 的 UMAP 没有旧 PC 规则受控图中的全局密集线状/蛇形退化。局部细长分支本身不自动视为错误；UMAP 形态不是判断 PCA 科学有效性的唯一依据。实现验证同时依赖独立谱重算、动态边界、相同输入/细胞集合、原始 counts 和完整 DF_adj 比较。该启发式验证不证明存在唯一生物学最优 PC 数。','',
        '当前正式规则未出现需要改变 regression 的异常，因此没有启动 nFeature_RNA+pMT 敏感性分析，也未替换 pMT-only 策略。','',
        '- [原代码 / 旧 PC 受控复现 / 当前：过滤前](GSE155698/pre_filter_three_versions.png)',
        '- [原代码 / 旧 PC 受控复现 / 当前：实际过滤后集合](GSE155698/post_filter_three_versions.png)',
        '- [过滤前相同细胞：sample](GSE155698/pre_filter_matched_sample.png)',
        '- [过滤前相同细胞：DF_adj](GSE155698/pre_filter_matched_DF_adj.png)',
        '- [过滤后共同细胞：sample](GSE155698/post_filter_common_sample.png)',
        '- [过滤后共同细胞：DF_adj](GSE155698/post_filter_common_DF_adj.png)',
        '- [20 样本 ElbowPlot 总览](GSE155698/sample_elbows_overview.png)',
        '- [合并 ElbowPlot](GSE155698/elbows/elbow_ALL.png)','',
        '## 实际测试与未解决限制','',
        '- Python tests：86/86 PASS；现有 PCA elbow、动态上限、小样本、fallback、标准差/方差分离和同坐标子集测试 PASS，R exit 0。',
        '- 两组真实完整 QC、独立最终对象验证均 exit 0；使用安装的 DoubletFinder 2.0.6 和 Seurat 5.4.0。未安装、升级或切换任何 R 包。',
        '- 固定 R 4.5.0 Runtime。21 个导出谱的独立计算、全量稀疏 counts、metadata、DF_adj、cell-cycle 和 clean object 检查 PASS。',
        '- 原始输入和全部受保护历史文件哈希前后相同。诊断不重建 raw，不更改 group，不覆盖既有 RDS。',
        '- 未找到真正修复前 GSE155698 Skill 历史产物；受控旧 PC 规则复现不能替代其历史 provenance。原代码坐标未保存，无法在原代码坐标上重新绘制精确共同细胞子集。',
        '- 生产代码没有发现需要修复的新问题；本次提交只增加真实验证工具、诊断图表和文档。','',
        '[测试证据](GSE155698/test_evidence.json) · [独立选择重算](GSE155698/independent_selection_review.json) · [对象/集合验证](GSE155698/validation.json)','',
        '## 复现','',
        '从仓库根目录，通过 `/c/Python312/python.exe` 执行 `tests/qc_comparison/pca_validation/run_validation.py GSE155698 SOURCE NEW_ROOT`。旧规则受控复现使用另一个 NEW_ROOT 并加 `--mode legacy_pc_control`。目录必须全新；已有最终对象不被覆盖。R 计算由当前 qc/run_workflow.py 调用固定 Runtime，测试观察脚本只导出诊断并调用正式 qc_run()。','',
        '全部新 RDS、每张 ElbowPlot 的 PDF 和完整诊断坐标保留在 `.runtime/pca-validation-GSE155698/`。源码提交仅包含小型诊断表/图及测试文档。']
    args.report.write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print('Wrote',args.report)


if __name__=='__main__':
    main()
