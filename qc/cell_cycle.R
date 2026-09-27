qc_cycle_features <- function(object) {
  reference <- Seurat::cc.genes
  available <- rownames(object)
  matched <- function(names) {
    indexes <- match(toupper(names),toupper(available))
    unique(available[stats::na.omit(indexes)])
  }
  list(s=matched(reference$s.genes),g2m=matched(reference$g2m.genes))
}

qc_score_cycle <- function(object,paths,report) {
  features <- qc_cycle_features(object)
  if(length(features$s)<5L || length(features$g2m)<5L)
    stop("Cell-cycle gene coverage is insufficient: S=",length(features$s),"; G2M=",length(features$g2m),". Review species/gene identifiers")
  qc_report_add(report,"CELL_CYCLE",metric="S_features_matched",value=length(features$s))
  qc_report_add(report,"CELL_CYCLE",metric="G2M_features_matched",value=length(features$g2m))
  normalized <- Seurat::NormalizeData(object,verbose=FALSE)
  scored <- Seurat::CellCycleScoring(normalized,s.features=features$s,g2m.features=features$g2m,set.ident=FALSE)
  scored$CC.Difference <- scored$S.Score-scored$G2M.Score
  qc_plot_cycle(scored,paths,report)
  for (phase in c("G1","S","G2M")) qc_report_add(report,"CELL_CYCLE",metric=paste0(phase,"_cells"),value=sum(scored$Phase==phase))
  # Rebuild from raw counts to remove the temporary normalized data layer and any analysis state.
  counts <- qc_counts(object)
  metadata <- scored@meta.data[colnames(counts),,drop=FALSE]
  final <- Seurat::CreateSeuratObject(counts=counts,min.cells=0,min.features=0,meta.data=metadata)
  if(length(final@reductions) || length(final@graphs) || length(final@neighbors)) stop("QC object retained downstream analysis state")
  if(!identical(qc_counts(final),counts)) stop("Raw counts changed during cell-cycle scoring")
  final
}
