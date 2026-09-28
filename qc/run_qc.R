# Invoke through runtime/r45/run_r45.py.
qc_runner_source <- if(sys.nframe() >= 1L && !is.null(sys.frame(1)$ofile)) sys.frame(1)$ofile else sub("^--file=","",commandArgs()[grepl("^--file=",commandArgs())][1])
qc_runner_dir <- dirname(normalizePath(qc_runner_source,winslash="/",mustWork=TRUE))
for (name in c("qc_utils.R","plot_utils.R","qc_precheck.R","qc_filter.R","doublet_finder.R","cell_cycle.R"))
  source(file.path(qc_runner_dir,name))

qc_run <- function(root,gse,package_available=function(p) requireNamespace(p,quietly=TRUE),df_api=NULL) {
  paths <- qc_paths(root,gse)
  if(file.exists(paths$final)) stop("seurat_qc.rds exists; archive it explicitly before rebuilding")
  status <- qc_precheck(root,gse,package_available)
  if(status!="PASS") stop("QC precheck stopped: ",status,"; inspect ",paths$report)
  prior <- utils::read.csv(paths$report,stringsAsFactors=FALSE,check.names=FALSE)
  prior$metric[prior$section=="RUN" & prior$metric=="status"] <- "precheck_status"
  report <- qc_report_new()
  report$rows <- lapply(seq_len(nrow(prior)),function(i) prior[i,,drop=FALSE])
  started <- proc.time()[["elapsed"]]
  result <- tryCatch({
    object <- readRDS(paths$raw)
    counts <- qc_validate_input(object,gse)
    decision <- qc_decision(paths,unique(as.character(object$sample)))
    cfg <- qc_config(paths,decision)
    raw_cells <- ncol(object)
    removed <- decision$removed_samples %||% character()
    if(length(removed)) {
      qc_report_add(report,"SAMPLE",metric="removed_samples",value=paste(removed,collapse=";"),
        status="USER_CONFIRMED_REMOVAL",note=decision$user_statement)
      object <- object[,!object$sample %in% removed]
    }
    active_cells <- ncol(object)
    object <- qc_add_percentages(object,cfg$species)
    qc_plot_metrics(object,paths,"before_QC","Before cell QC",report)
    filtered <- qc_filter_cells(object,cfg,report)
    after_cell_qc <- ncol(filtered)
    qc_plot_metrics(filtered,paths,"after_QC","After cell QC",report)
    rm(object,counts); gc()
    actual_api <- if(is.null(df_api)) qc_default_doublet_api() else df_api
    doubled <- qc_run_doublets(filtered,cfg,decision,paths,report,actual_api)
    rm(filtered); gc()
    final <- qc_score_cycle(doubled$object,paths,report)
    allowed <- if(doubled$unevaluated>0L) c("Singlet","NotEvaluated") else "Singlet"
    if(any(!final$DF_adj %in% allowed)) stop("Final object contains disallowed DF_adj calls")
    if(any(final$DF_adj=="NotEvaluated" & final$doublet_status!="SKIPPED_LOW_CELL")) stop("Unassessed cells lack the approved skip marker")
    if(!all(c("pMT","pRP","pHB","S.Score","G2M.Score","Phase","CC.Difference","pANN","DF","DF_adj","doublet_status") %in% names(final@meta.data))) stop("Final QC metadata incomplete")
    if(!identical(rownames(final@meta.data),SeuratObject::Cells(final))) stop("Final cell/metadata alignment changed")
    for (field in c("sample","database","group")) if(any(qc_blank(final@meta.data[[field]]))) stop("Final metadata missing ",field)
    finish_status <- if(doubled$unevaluated>0L) "QC_COMPLETE_WITH_UNEVALUATED_DOUBLETS" else "COMPLETE_QC"
    fields <- list(raw_cells=raw_cells,after_cell_qc_cells=after_cell_qc,
      after_doublet_cells=ncol(final),removed_by_cell_qc=active_cells-after_cell_qc,
      removed_doublets=doubled$pre_filter_cells-ncol(final),final_cells=ncol(final),final_genes=nrow(final),
      unevaluated_doublet_cells=doubled$unevaluated,runtime_seconds=round(proc.time()[["elapsed"]]-started,3))
    for(field in names(fields)) qc_report_add(report,"RUN",metric=field,value=fields[[field]])
    if(doubled$unevaluated>0L) qc_report_add(report,"RUN",metric="unevaluated_samples",
      value=paste(unique(final$sample[final$DF_adj=="NotEvaluated"]),collapse=";"),
      status="USER_CONFIRMED_SKIP",note=decision$user_statement)
    qc_report_add(report,"RUN",metric="status",value=finish_status,
      status=if(finish_status=="COMPLETE_QC") "PASS" else finish_status)
    required <- as.vector(outer(c("before_QC","after_QC","doublet_umap","doublet_vlnplot",
                                   "cell_cycle_phase","cell_cycle_score"),c(".pdf",".png"),paste0))
    files <- file.path(paths$qc,required)
    if(any(!file.exists(files)) || any(file.info(files)$size<=0)) stop("Required QC PDF/PNG missing or empty")
    pending <- paste0(paths$final,".pending")
    saveRDS(final,pending)
    serialized <- readRDS(pending)
    if(!identical(qc_counts(serialized),qc_counts(final))) stop("Serialized QC counts changed")
    if(file.exists(paths$final)) stop("seurat_qc.rds appeared during QC; refusing to overwrite")
    if(!file.rename(pending,paths$final)) stop("Could not finalize seurat_qc.rds")
    qc_report_write(report,paths$report)
    finish_status
  },error=function(e) {
    report$rows <- Filter(function(row) !(row$section=="RUN" && row$metric=="status"),report$rows)
    qc_report_add(report,"RUN",metric="status",value="FAILED",status="FAILED_INPUT",note=conditionMessage(e))
    qc_report_write(report,paths$report)
    stop(e)
  })
  result
}

if (sys.nframe()==0L) {
  args <- commandArgs(trailingOnly=TRUE)
  if(length(args)!=2L) stop("Usage: run_r45.py qc/run_qc.R REPOSITORY_ROOT GSE")
  cat("QC run:",qc_run(args[1],args[2]),"\n")
}
