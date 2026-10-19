# DecontX estimates contamination from raw RNA counts. No cell removal is applied.
qc_run_ambient <- function(object,paths,report,api=NULL,seed=42L,cfg=NULL) {
  if("decontX" %in% names(object@assays)) stop("Existing decontX assay preserved; use read-only ambient review")
  cfg <- qc_ambient_config(cfg)
  original <- qc_counts(object)
  original_meta <- object@meta.data
  object <- qc_capture_libraries(object)
  if(any(Matrix::colSums(original)<=0)) stop("DECONTX_INPUT: zero-count cell found; preserve QC-2 and review input")
  if(is.null(api)) {
    needed <- c("celda","SingleCellExperiment","SummarizedExperiment")
    missing <- needed[!vapply(needed,requireNamespace,logical(1),quietly=TRUE)]
    if(length(missing)) stop("DECONTX_DEPENDENCY: missing ",paste(missing,collapse=", "),"; QC-2 remains available")
    api <- list(run=function(sce,batch,seed,...) celda::decontX(sce,batch=batch,seed=seed,verbose=FALSE,...),
                counts=function(sce) celda::decontXcounts(sce))
  }
  sce <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=original),
    colData=S4Vectors::DataFrame(object@meta.data))
  background <- qc_ambient_background(cfg$background_manifest,object)
  arguments <- list(sce=sce,batch=as.character(object$library_id),seed=seed)
  if(!is.null(background)) arguments <- c(arguments,list(background=background$counts,background_batch=background$batch))
  corrected <- do.call(api$run,arguments)
  contamination <- as.numeric(SummarizedExperiment::colData(corrected)$decontX_contamination)
  corrected_counts <- api$counts(corrected)
  if(!inherits(corrected_counts,"sparseMatrix")) corrected_counts <- Matrix::Matrix(corrected_counts,sparse=TRUE)
  if(!identical(dimnames(corrected_counts),dimnames(original)) ||
     any(!is.finite(corrected_counts@x)) || any(corrected_counts@x<0) ||
     length(contamination)!=ncol(object) || any(!is.finite(contamination) | contamination<0 | contamination>1))
    stop("DECONTX_RESULT: matrix alignment or contamination validation failed")
  object[["decontX"]] <- SeuratObject::CreateAssay5Object(counts=corrected_counts)
  object$decontX_contamination <- contamination
  object$decontX_original_UMI <- as.numeric(Matrix::colSums(original))
  object$decontX_corrected_UMI <- as.numeric(Matrix::colSums(corrected_counts))
  object$decontX_removed_UMI <- object$decontX_original_UMI-object$decontX_corrected_UMI
  SeuratObject::DefaultAssay(object) <- "RNA"
  if(!identical(qc_counts(object),original)) stop("DECONTX_RESULT: original RNA counts changed")
  if(!identical(object@meta.data[,names(original_meta),drop=FALSE],original_meta))
    stop("DECONTX_RESULT: original QC metadata changed: ",paste(names(original_meta)[!vapply(names(original_meta),
      function(field) identical(object@meta.data[[field]],original_meta[[field]]),logical(1))],collapse=","))
  object@misc$ambient_RNA <- list(method="celda::decontX",input_assay="RNA",input_layer="counts",
    output_assay="decontX",output_layer="counts",batch="library_id",seed=seed,
    cell_removal="none",default_analysis_assay="decontX",diagnostic_thresholds=cfg,
    background=if(is.null(background)) "not supplied" else background[c("manifest","manifest_md5","sources")],
    metric_source="QC metadata from RNA; decontX_*_UMI from corrected counts",package_version=as.character(utils::packageVersion("celda")))
  cells <- data.frame(cell_id=colnames(object),sample=as.character(object$sample),library_id=as.character(object$library_id),
    contamination=contamination,original_UMI=object$decontX_original_UMI,
    corrected_UMI=object$decontX_corrected_UMI,removed_UMI=object$decontX_removed_UMI)
  utils::write.csv(cells,file.path(paths$qc,"decontX_cells.csv"),row.names=FALSE)
  summary <- do.call(rbind,lapply(unique(cells$sample),function(sample) {
    x <- cells[cells$sample==sample,]
    data.frame(sample=sample,cells=nrow(x),mean_contamination=mean(x$contamination),median_contamination=stats::median(x$contamination),
      p90_contamination=as.numeric(stats::quantile(x$contamination,.9)),original_UMI=sum(x$original_UMI),
      corrected_UMI=sum(x$corrected_UMI),removed_UMI=sum(x$removed_UMI))
  }))
  utils::write.csv(summary,file.path(paths$qc,"contamination_by_sample.csv"),row.names=FALSE)
  utils::write.csv(data.frame(input_cells=ncol(object),output_cells=ncol(object),genes=nrow(object),
    original_UMI=sum(original),corrected_UMI=sum(corrected_counts),mean_contamination=mean(contamination),
    median_contamination=stats::median(contamination),source_assay="RNA",corrected_assay="decontX",cells_removed=0L),
    file.path(paths$qc,"decontX_summary.csv"),row.names=FALSE)
  distribution <- ggplot2::ggplot(cells,ggplot2::aes(x=sample,y=contamination,fill=sample))+
    ggplot2::geom_violin(scale="width")+ggplot2::geom_boxplot(width=.15,outlier.shape=NA,fill="white")+
    ggplot2::labs(x="Sample",y="Estimated contamination fraction",title="DecontX contamination; cells retained")+
    theme_sc()+ggplot2::theme(legend.position="none",axis.text.x=ggplot2::element_text(angle=50,hjust=1))
  save_sc_plot(distribution,"decontX_contamination","violin",paths$qc,ncol(object),length(unique(object$sample)),report)
  comparison <- ggplot2::ggplot(cells,ggplot2::aes(x=original_UMI,y=corrected_UMI,color=sample))+
    ggplot2::geom_point(size=get_umap_point_size(ncol(object)),alpha=.6)+
    ggplot2::geom_abline(slope=1,intercept=0,linetype=2,color="grey50")+
    ggplot2::labs(x="Original RNA counts per cell",y="DecontX corrected counts per cell",title="Expression counts before and after correction")+theme_sc()
  save_sc_plot(comparison,"decontX_comparison","diagnostics",paths$qc,ncol(object),length(unique(object$sample)),report)
  qc_report_add(report,"AMBIENT_RNA",metric="cells_removed",value=0L,note="Evaluation and correction only; no fixed contamination threshold")
  qc_report_add(report,"AMBIENT_RNA",metric="mean_contamination",value=mean(contamination))
  qc_review_ambient(object,paths,report,cfg)
  object
}
