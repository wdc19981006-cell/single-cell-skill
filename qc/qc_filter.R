qc_write_filter_reasons <- function(meta,cfg,paths,removed_samples=character()) {
  reasons <- qc_filter_reasons(meta,cfg$cell_qc,removed_samples)
  utils::write.csv(reasons,file.path(paths$qc,"qc_filter_reasons.csv"),row.names=FALSE)
  size <- sc_plot_dimensions("diagnostics")
  grDevices::pdf(file.path(paths$qc,"qc_filter_reasons.pdf"),width=size[1],height=size[2])
  on.exit(grDevices::dev.off(),add=TRUE)
  for(i in seq_len(nrow(reasons))) {
    fields <- c(paste0(names(qc_cell_filter_rules(meta,cfg$cell_qc)),"_hits"),
      "multiple_conditions","total_removed_cells","user_removed_sample_cells")
    frame <- data.frame(reason=factor(fields,levels=rev(fields)),cells=as.numeric(unlist(reasons[i,fields],use.names=FALSE)))
    print(ggplot2::ggplot(frame,ggplot2::aes(x=reason,y=cells))+ggplot2::geom_col(fill="#427A9A")+
      ggplot2::coord_flip()+theme_sc()+ggplot2::labs(x=NULL,y="Cells",title=paste("Filter reasons:",reasons$sample[i]),
        subtitle="Rule hits overlap; total_removed_cells counts each removed cell once"))
  }
  invisible(reasons)
}

qc_filter_cells <- function(object,cfg,report) {
  meta <- object@meta.data
  thresholds <- cfg$cell_qc
  keep <- qc_cell_keep(meta,thresholds)
  for (sample in unique(as.character(meta$sample))) {
    selected <- meta$sample==sample
    qc_report_add(report,"CELL_QC",sample,"input_cells",sum(selected))
    qc_report_add(report,"CELL_QC",sample,"after_qc_cells",sum(selected & keep))
    qc_report_add(report,"CELL_QC",sample,"removed_cells",sum(selected & !keep))
  }
  if(!any(keep)) stop("Cell QC removed every cell; review thresholds")
  remaining <- table(as.character(meta$sample[keep]))
  if(!setequal(names(remaining),unique(as.character(meta$sample)))) stop("Cell QC emptied one or more samples; requires user decision")
  filtered <- object[,rownames(meta)[keep]]
  counts <- qc_counts(filtered)
  expressed <- Matrix::rowSums(counts>0)
  genes <- expressed>3L
  if(!any(genes)) stop("Gene QC removed every gene; review counts and retained-cell coverage")
  qc_report_add(report,"GENE_QC",metric="input_genes",value=nrow(object))
  qc_report_add(report,"GENE_QC",metric="retained_genes",value=sum(genes))
  qc_report_add(report,"GENE_QC",metric="removed_genes",value=sum(!genes))
  qc_report_add(report,"GENE_QC",metric="rule",value="Matrix::rowSums(counts > 0) > 3")
  qc_refresh_metrics(filtered[rownames(counts)[genes],],preserve=TRUE)
}
