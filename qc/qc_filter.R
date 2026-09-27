qc_filter_cells <- function(object,cfg,report) {
  meta <- object@meta.data
  thresholds <- cfg$cell_qc
  keep <- meta$nFeature_RNA >= thresholds$nFeature_min & meta$nFeature_RNA <= thresholds$nFeature_max &
    meta$nCount_RNA >= thresholds$nCount_min & meta$nCount_RNA <= thresholds$nCount_max &
    meta$pMT <= thresholds$pMT_max & meta$pHB <= thresholds$pHB_max & meta$pRP <= thresholds$pRP_max
  for (sample in unique(as.character(meta$sample))) {
    selected <- meta$sample==sample
    qc_report_add(report,"CELL_QC",sample,"input_cells",sum(selected))
    qc_report_add(report,"CELL_QC",sample,"after_qc_cells",sum(selected & keep))
    qc_report_add(report,"CELL_QC",sample,"removed_cells",sum(selected & !keep))
  }
  if(!any(keep)) stop("Cell QC removed every cell; review thresholds")
  remaining <- table(as.character(meta$sample[keep]))
  if(!setequal(names(remaining),unique(as.character(meta$sample)))) stop("Cell QC emptied one or more samples; requires user decision")
  # Subset columns only; do not apply a second min.cells gene filter.
  object[,rownames(meta)[keep]]
}
