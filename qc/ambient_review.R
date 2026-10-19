# Read-only diagnostics, also callable on an existing QC-3 without rerunning DecontX.
qc_ambient_config <- function(cfg=NULL) {
  defaults <- yaml::read_yaml(file.path(qc_script_dir,"default_config.yaml"))$ambient
  if(is.null(cfg)) cfg <- defaults else cfg <- utils::modifyList(defaults,cfg)
  for(key in c("mean_warning","median_warning","p95_warning","high_contamination","marker_retention_warning"))
    if(!is.numeric(cfg[[key]]) || length(cfg[[key]])!=1L || !is.finite(cfg[[key]]) ||
       cfg[[key]]<0 || cfg[[key]]>1) stop("Invalid ambient diagnostic threshold: ",key)
  cfg
}

qc_ambient_background <- function(manifest,object) {
  if(is.null(manifest)) return(NULL)
  entries <- utils::read.csv(manifest,stringsAsFactors=FALSE)
  if(!all(c("library_id","path","evidence") %in% names(entries)) ||
     anyDuplicated(entries$library_id) || any(qc_blank(entries$evidence)) ||
     !setequal(entries$library_id,unique(object$library_id)))
    stop("Background manifest must document genuine empty droplets for every capture library")
  blocks <- lapply(seq_len(nrow(entries)),function(i) {
    x <- readRDS(entries$path[i])
    if(!inherits(x,"sparseMatrix") || is.null(rownames(x)) || anyDuplicated(rownames(x)) ||
       !all(rownames(qc_counts(object)) %in% rownames(x)) ||
       !ncol(x) || any(!is.finite(x@x) | x@x<0 | x@x!=round(x@x)) ||
       is.null(colnames(x)) || anyDuplicated(colnames(x)) || any(qc_blank(colnames(x))) ||
       any(colnames(x) %in% colnames(object))) stop("Invalid or non-empty-droplet background matrix")
    x <- x[rownames(qc_counts(object)),,drop=FALSE]
    colnames(x) <- paste(entries$library_id[i],colnames(x),sep="::background::")
    x
  })
  list(counts=do.call(cbind,blocks),batch=rep(entries$library_id,vapply(blocks,ncol,integer(1))),
       manifest=normalizePath(manifest,winslash="/"),manifest_md5=unname(tools::md5sum(manifest)),
       sources=data.frame(entries,md5=unname(tools::md5sum(entries$path))))
}

qc_review_ambient <- function(object,paths,report,cfg=NULL) {
  cfg <- qc_ambient_config(cfg)
  original <- qc_counts(object)
  corrected <- SeuratObject::LayerData(object,assay="decontX",layer="counts")
  if(!inherits(corrected,"sparseMatrix") || !identical(dimnames(original),dimnames(corrected)) ||
     any(!is.finite(corrected@x) | corrected@x<0)) stop("Invalid corrected matrix in ambient review")
  contamination <- object$decontX_contamination
  if(length(contamination)!=ncol(object) || any(!is.finite(contamination) | contamination<0 | contamination>1))
    stop("Invalid contamination estimates in ambient review")
  cells <- data.frame(cell_id=colnames(object),sample=as.character(object$sample),
    contamination=contamination,original_UMI=as.numeric(Matrix::colSums(original)),
    corrected_UMI=as.numeric(Matrix::colSums(corrected)),nFeature_RNA=as.numeric(Matrix::colSums(original>0)))
  rows <- lapply(unique(cells$sample),function(sample) {
    x <- cells[cells$sample==sample,]; high <- x$contamination>=cfg$high_contamination
    q <- as.numeric(stats::quantile(x$contamination,c(.9,.95,.99)))
    data.frame(sample=sample,cells=nrow(x),mean_contamination=mean(x$contamination),
      median_contamination=stats::median(x$contamination),p90_contamination=q[1],p95_contamination=q[2],p99_contamination=q[3],
      high_cells=sum(high),high_median_UMI=if(any(high)) stats::median(x$original_UMI[high]) else NA_real_,
      high_median_nFeature=if(any(high)) stats::median(x$nFeature_RNA[high]) else NA_real_,
      original_UMI=sum(x$original_UMI),corrected_UMI=sum(x$corrected_UMI),
      removed_UMI=sum(x$original_UMI-x$corrected_UMI),
      retained_UMI_fraction=sum(x$corrected_UMI)/sum(x$original_UMI))
  })
  summary <- do.call(rbind,rows)
  utils::write.csv(summary,file.path(paths$qc,"contamination_by_sample.csv"),row.names=FALSE)
  warnings <- data.frame(sample=character(),category=character(),severity=character(),message=character(),recommendation=character())
  add <- function(sample,category,message) {
    warnings <<- rbind(warnings,data.frame(sample=sample,category=category,severity="WARNING",message=message,
      recommendation="Review contamination, UMI/nFeature and marker evidence; no cells were removed"))
  }
  for(i in seq_len(nrow(summary))) {
    x <- summary[i,]
    if(x$mean_contamination>=cfg$mean_warning || x$median_contamination>=cfg$median_warning || x$p95_contamination>=cfg$p95_warning)
      add(x$sample,"EXTREME_CONTAMINATION",paste("mean",signif(x$mean_contamination,3),"median",signif(x$median_contamination,3),
        "P95",signif(x$p95_contamination,3),"retained UMI",signif(x$retained_UMI_fraction,3)))
    if(any(cells$corrected_UMI[cells$sample==x$sample]<=0)) add(x$sample,"ZERO_CORRECTED_UMI","Corrected zero-count cells; formal analysis must stop")
  }
  marker_rows <- list()
  size <- sc_plot_dimensions("diagnostics")
  grDevices::pdf(file.path(paths$qc,"marker_correction_review.pdf"),width=size[1],height=size[2])
  tryCatch({
    for(sample in unique(cells$sample)) {
      index <- cells$sample==sample
      for(marker in cfg$markers) {
        present <- marker %in% rownames(original)
        before <- if(present) as.numeric(original[marker,index]) else NA_real_
        after <- if(present) as.numeric(corrected[marker,index]) else NA_real_
        retention <- if(present && sum(before)>0) sum(after)/sum(before) else NA_real_
        marker_rows[[length(marker_rows)+1L]] <- data.frame(sample=sample,marker=marker,
          status=if(present) "MEASURED" else "MARKER_ABSENT",mean_original=mean(before),mean_corrected=mean(after),
          detected_original=mean(before>0),detected_corrected=mean(after>0),retained_fraction=retention)
        if(is.finite(retention) && retention<cfg$marker_retention_warning)
          add(sample,"MARKER_LOSS",paste(marker,"retained count fraction",signif(retention,3),"; possible overcorrection requires biological review"))
      }
      x <- cells[index,]; x$diagnostic_stratum <- ifelse(x$contamination>=cfg$high_contamination,"High estimate","Other")
      panels <- lapply(c("original_UMI","nFeature_RNA"),function(field)
        ggplot2::ggplot(x,ggplot2::aes(x=diagnostic_stratum,y=.data[[field]],fill=diagnostic_stratum))+
          ggplot2::geom_boxplot()+theme_sc()+ggplot2::labs(title=paste(sample,field),x=NULL,y=field)+ggplot2::theme(legend.position="none"))
      print(patchwork::wrap_plots(panels,ncol=2))
    }
    markers <- if(length(marker_rows)) do.call(rbind,marker_rows) else data.frame(sample=character(),marker=character(),status=character())
    for(sample in unique(cells$sample)) {
      x <- markers[markers$sample==sample & markers$status=="MEASURED",,drop=FALSE]
      if(nrow(x)) print(ggplot2::ggplot(x,ggplot2::aes(x=mean_original,y=mean_corrected,label=marker))+
        ggplot2::geom_point()+ggplot2::geom_text(check_overlap=TRUE,vjust=-.5)+
        ggplot2::geom_abline(slope=1,intercept=0,linetype=2)+theme_sc()+
        ggplot2::labs(title=paste("Marker correction:",sample),x="Mean original counts",y="Mean corrected counts"))
      else { grid::grid.newpage(); grid::grid.text(paste(sample,": configured markers absent; supply tissue-relevant markers")) }
    }
  },finally=grDevices::dev.off())
  utils::write.csv(markers,file.path(paths$qc,"marker_correction_review.csv"),row.names=FALSE,na="")
  utils::write.csv(warnings,file.path(paths$qc,"contamination_warnings.csv"),row.names=FALSE)
  yaml::write_yaml(list(diagnostic_thresholds=cfg,cell_removal="none",metric_source="RNA counts; corrected metrics explicitly separate"),
    file.path(paths$qc,"ambient_review_config.yaml"))
  if(nrow(warnings)) qc_report_add(report,"AMBIENT_RNA",metric="diagnostic_warnings",value=nrow(warnings),status="WARNING")
  invisible(list(summary=summary,warnings=warnings))
}
