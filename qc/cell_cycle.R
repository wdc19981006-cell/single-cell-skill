qc_score_cycle <- function(object,paths,report) {
  features <- qc_cycle_features(object)
  if(length(features$s)<5L || length(features$g2m)<5L)
    stop("Cell-cycle gene coverage is insufficient: S=",length(features$s),"; G2M=",length(features$g2m),". Review species/gene identifiers")
  qc_report_add(report,"CELL_CYCLE",metric="S_features_matched",value=length(features$s))
  qc_report_add(report,"CELL_CYCLE",metric="G2M_features_matched",value=length(features$g2m))
  normalized <- Seurat::NormalizeData(object,verbose=FALSE)
  scored <- Seurat::CellCycleScoring(normalized,s.features=features$s,g2m.features=features$g2m,set.ident=FALSE)
  scored$CC.Difference <- scored$S.Score-scored$G2M.Score
  # This object is for spatial diagnostics only; the final object below is
  # reconstructed from the untouched raw counts and scored metadata.
  plot_object <- qc_prepare_cycle_plot_object(scored)
  qc_plot_cycle(plot_object,paths,report)
  rm(plot_object); gc()
  for (phase in c("G1","S","G2M")) qc_report_add(report,"CELL_CYCLE",metric=paste0(phase,"_cells"),value=sum(scored$Phase==phase))
  # Rebuild from raw counts to remove the temporary normalized data layer and any analysis state.
  counts <- qc_counts(object)
  metadata <- scored@meta.data[colnames(counts),,drop=FALSE]
  final <- Seurat::CreateSeuratObject(counts=counts,min.cells=0,min.features=0,meta.data=metadata)
  if(length(final@reductions) || length(final@graphs) || length(final@neighbors)) stop("QC object retained downstream analysis state")
  if(!identical(qc_counts(final),counts)) stop("Raw counts changed during cell-cycle scoring")
  final
}

qc_prepare_cycle_plot_object <- function(scored) {
  plot_object <- Seurat::NormalizeData(scored,verbose=FALSE)
  plot_object <- Seurat::FindVariableFeatures(plot_object,nfeatures=min(2000L,nrow(plot_object)-1L),verbose=FALSE)
  max_pc <- min(20L,ncol(plot_object)-1L,length(Seurat::VariableFeatures(plot_object))-1L,nrow(plot_object)-1L)
  if (max_pc<2L) stop("Insufficient cells/features for cell-cycle UMAP")
  plot_object <- Seurat::ScaleData(plot_object,verbose=FALSE)
  plot_object <- Seurat::RunPCA(plot_object,npcs=max_pc,verbose=FALSE)
  plot_object <- Seurat::RunUMAP(plot_object,dims=seq_len(max_pc),verbose=FALSE)
  plot_object
}
