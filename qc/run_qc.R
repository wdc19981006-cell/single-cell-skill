# Invoke through runtime/r45/run_r45.py.
qc_runner_files <- unlist(lapply(sys.frames(),function(frame) frame$ofile),use.names=FALSE)
qc_runner_source <- if(length(qc_runner_files)) tail(qc_runner_files,1L) else sub("^--file=","",commandArgs()[grepl("^--file=",commandArgs())][1])
qc_runner_dir <- dirname(normalizePath(qc_runner_source,winslash="/",mustWork=TRUE))
for(name in c("qc_utils.R","plot_utils.R","qc_precheck.R","qc_filter.R","doublet_finder.R","cell_cycle.R","ambient_rna.R"))
  source(file.path(qc_runner_dir,name))

qc_md5 <- function(path) unname(tools::md5sum(path))
qc_stage_log <- function(paths,step,text) {
  cat(format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z"),text,"\n",file=file.path(paths$qc,paste0(step,".log")),append=TRUE)
}
qc_stage_validate <- function(object,step,gse,input=NULL,paths=NULL,legacy=FALSE) {
  counts <- qc_validate_input(object,gse)
  if(any(abs(object$nCount_RNA-as.numeric(Matrix::colSums(counts)))>1e-8) ||
     any(object$nFeature_RNA!=as.numeric(Matrix::colSums(counts>0)))) stop("QC checkpoint metrics disagree with current RNA counts")
  if(length(object@reductions) || length(object@graphs) || length(object@neighbors)) stop("QC checkpoint contains temporary analysis state")
  if(step %in% c("qc2","qc3")) {
    required <- c("pMT","pRP","pHB","S.Score","G2M.Score","Phase","CC.Difference","pANN","DF","DF_adj","doublet_status")
    if(!all(required %in% names(object@meta.data))) stop("QC-2 metadata incomplete")
    if(any(!object$DF_adj %in% c("Singlet","NotEvaluated"))) stop("QC checkpoint contains rejected DF_adj cells")
    if(any(object$DF_adj=="NotEvaluated" & object$doublet_status!="SKIPPED_LOW_CELL")) stop("Unassessed cells lack approved skip marker")
  }
  if(step=="qc3") {
    if(!"decontX" %in% names(object@assays) || !"decontX_contamination" %in% names(object@meta.data)) stop("DecontX checkpoint incomplete")
    corrected <- SeuratObject::LayerData(object,assay="decontX",layer="counts")
    if(!identical(dimnames(corrected),dimnames(counts)) || any(!is.finite(corrected@x) | corrected@x<0)) stop("Invalid DecontX corrected counts")
    if(any(!is.finite(object$decontX_contamination) | object$decontX_contamination<0 | object$decontX_contamination>1)) stop("Invalid DecontX contamination estimates")
    if(SeuratObject::DefaultAssay(object)!="RNA") stop("QC-3 must preserve RNA as default assay")
  } else if(!identical(names(object@assays),"RNA") || !identical(SeuratObject::Layers(object[["RNA"]]),"counts")) stop("QC-1/2 checkpoint must contain raw RNA counts only")
  provenance <- object@misc$qc_pipeline[[step]]
  if(!legacy && (is.null(provenance) || !identical(provenance$stage,step))) stop("QC checkpoint provenance missing")
  if(!is.null(input) && !legacy) {
    if(!file.exists(input) || !identical(qc_md5(input),provenance$input_md5)) stop("QC checkpoint input differs; refusing to replace existing valid output")
    parent <- qc_counts(readRDS(input))
    if(!all(rownames(counts) %in% rownames(parent)) || !all(colnames(counts) %in% colnames(parent)) ||
       !identical(counts,parent[rownames(counts),colnames(counts),drop=FALSE])) stop("QC checkpoint changed original RNA counts")
    if(step=="qc3" && !identical(counts,parent)) stop("QC-3 changed cells or original RNA counts")
  }
  if(!is.null(paths) && !legacy) {
    samples <- if(file.exists(paths$raw)) unique(as.character(readRDS(paths$raw)$sample)) else unique(as.character(object$sample))
    current_decision <- qc_decision(paths,samples)
    if(!identical(provenance$configuration,qc_config(paths,current_decision)) ||
       (!is.null(provenance$user_decision) && !identical(provenance$user_decision,current_decision)))
      stop("QC checkpoint configuration/decision differs; preserve existing result and use a separate output root")
    stems <- switch(step,qc1=c("before_QC","after_QC"),qc2=c("doublet_umap","doublet_vlnplot","singlet_dimplot","cell_cycle_phase","cell_cycle_score"),qc3=c("decontX_contamination","decontX_comparison"))
    files <- file.path(paths$qc,c(as.vector(outer(stems,c(".pdf",".png"),paste0)),"qc_report.csv",paste0(step,"_config.yaml"),paste0(step,".log")))
    if(step=="qc3") files <- c(files,file.path(paths$qc,"decontX_config.yaml"))
    if(any(!file.exists(files)) || any(file.info(files)$size<=0)) stop("QC checkpoint outputs missing or empty")
    expected_tables <- switch(step,qc1=c("cell_QC_all.csv","gene_filter_summary.csv","sample_QC_summary.csv"),
      qc2=c("doublet_cells_all.csv","doublet_summary.csv","cell_cycle_summary.csv","cell_cycle_cells.csv"),
      qc3=c("decontX_cells.csv","decontX_summary.csv","contamination_by_sample.csv"))
    for(table in expected_tables) {
      file <- file.path(paths$qc,table)
      if(!file.exists(file) || !nrow(utils::read.csv(file,check.names=FALSE))) stop("QC checkpoint table missing or empty: ",table)
    }
  }
  invisible(object)
}

qc_stage_save <- function(object,paths,step,input,cfg,report,status,decision=list()) {
  if(file.exists(paths$final)) stop("QC checkpoint already exists; refusing overwrite: ",paths$final)
  object@misc$qc_pipeline[[step]] <- list(stage=step,input=normalizePath(input,winslash="/"),input_md5=qc_md5(input),
    output=paths$final,configuration=cfg,user_decision=decision,R=R.version.string,completed_at=format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z"))
  qc_report_add(report,"RUN",metric="input_file",value=input)
  qc_report_add(report,"RUN",metric="input_md5",value=qc_md5(input))
  qc_report_add(report,"RUN",metric="final_cells",value=ncol(object))
  qc_report_add(report,"RUN",metric="final_genes",value=nrow(object))
  qc_report_add(report,"RUN",metric="status",value=status,status=status)
  qc_report_write(report,paths$report)
  yaml::write_yaml(list(step=step,input=input,input_md5=qc_md5(input),output=paths$final,parameters=cfg),file.path(paths$qc,paste0(step,"_config.yaml")))
  if(step=="qc3") yaml::write_yaml(object@misc$ambient_RNA,file.path(paths$qc,"decontX_config.yaml"))
  qc_stage_log(paths,step,paste("COMPLETE",status,"cells",ncol(object),"genes",nrow(object)))
  pending <- tempfile(pattern=paste0(basename(paths$final),"."),tmpdir=paths$qc,fileext=".pending")
  saveRDS(object,pending)
  verified <- readRDS(pending)
  qc_stage_validate(verified,step,paths$gse,input,paths)
  if(file.exists(paths$final) || !file.rename(pending,paths$final)) stop("Could not finalize QC checkpoint without overwrite")
  invisible(verified)
}

qc_basic_tables <- function(before,after,cfg,paths,report) {
  meta <- before@meta.data
  rules <- qc_cell_filter_rules(meta,cfg$cell_qc)
  table <- data.frame(cell_id=rownames(meta),meta,row.names=NULL,check.names=FALSE)
  table$retained <- rownames(meta) %in% colnames(after)
  table$failed_rules <- vapply(seq_len(nrow(meta)),function(i) paste(names(rules)[vapply(rules,function(x) x[i],logical(1))],collapse=";"),character(1))
  table$removal_reason <- ifelse(table$retained,"retained",ifelse(nzchar(table$failed_rules),table$failed_rules,"user_confirmed_sample_removal"))
  utils::write.csv(table,file.path(paths$qc,"cell_QC_all.csv"),row.names=FALSE,na="")
  selected <- colnames(after)
  expressed <- Matrix::rowSums(qc_counts(before)[,selected,drop=FALSE]>0)
  utils::write.csv(data.frame(gene=names(expressed),expressed_retained_cells=as.numeric(expressed),retained=expressed>3L),
    file.path(paths$qc,"gene_filter_summary.csv"),row.names=FALSE)
  summary <- do.call(rbind,lapply(unique(as.character(meta$sample)),function(sample) {
    input <- sum(meta$sample==sample); retained <- sum(after$sample==sample)
    data.frame(sample=sample,input_cells=input,retained_cells=retained,removed_cells=input-retained,
      cell_retention=retained/input,input_genes=nrow(before),retained_genes=nrow(after),gene_retention=nrow(after)/nrow(before))
  }))
  utils::write.csv(summary,file.path(paths$qc,"sample_QC_summary.csv"),row.names=FALSE)
}

qc_run_stage <- function(root,gse,step,package_available=function(p) requireNamespace(p,quietly=TRUE),df_api=NULL,ambient_api=NULL,input=NULL) {
  base <- qc_paths(root,gse); paths <- qc_stage_paths(base,step)
  if(is.null(input)) input <- switch(step,qc1=base$raw,qc2=qc_stage_paths(base,"qc1")$final,
    qc3=if(file.exists(base$final)) base$final else base$legacy_final)
  if(file.exists(paths$final)) {
    object <- readRDS(paths$final)
    qc_stage_validate(object,step,gse,input,paths)
    qc_stage_log(paths,step,"VALIDATED existing checkpoint; no processing repeated")
    return(object)
  }
  if(!file.exists(input)) stop("QC prerequisite missing: ",input,"; run the preceding step or supply --input")
  report <- qc_report_new(); qc_stage_log(paths,step,paste("START input",input))
  tryCatch({
    if(step=="qc1") {
      status <- qc_precheck(root,gse,package_available,step="qc1")
      if(status!="PASS") stop("QC precheck stopped: ",status,"; inspect ",base$report)
    }
    object <- readRDS(input); qc_validate_input(object,gse)
    samples <- if(file.exists(base$raw)) unique(as.character(readRDS(base$raw)$sample)) else unique(as.character(object$sample))
    decision <- qc_decision(base,samples); cfg <- qc_config(base,decision)
    qc_report_config(report,cfg)
    yaml::write_yaml(list(step=step,status="STARTED",input=input,input_md5=qc_md5(input),output=paths$final,parameters=cfg),
      file.path(paths$qc,paste0(step,"_config.yaml")))
    if(step=="qc1") {
      raw_cells <- ncol(object); raw_genes <- nrow(object)
      removed <- decision$removed_samples %||% character()
      object <- qc_add_percentages(object,cfg$species)
      qc_plot_metrics(object,paths,"before_QC","Before QC",report)
      active <- if(length(removed)) object[,!object$sample %in% removed] else object
      final <- qc_filter_cells(active,cfg,report)
      qc_basic_tables(object,final,cfg,paths,report)
      final <- qc_add_percentages(final,cfg$species)
      qc_plot_metrics(final,paths,"after_QC","After QC",report)
      qc_report_add(report,"RUN",metric="raw_cells",value=raw_cells)
      qc_report_add(report,"RUN",metric="raw_genes",value=raw_genes)
      qc_report_add(report,"RUN",metric="user_removed_samples",value=paste(removed,collapse=";"),status=if(length(removed)) "USER_CONFIRMED_REMOVAL" else "PASS")
      status <- "COMPLETE_QC1"
    } else if(step=="qc2") {
      qc_stage_validate(object,"qc1",gse,legacy=is.null(object@misc$qc_pipeline$qc1))
      object <- qc_capture_libraries(object)
      skipped <- decision$kept_low_cell_samples %||% character()
      for(library in unique(object$library_id)) {
        selected <- object$library_id==library
        samples <- unique(as.character(object$sample[selected]))
        if(sum(selected)<cfg$doublet$min_cells && !all(samples %in% skipped))
          stop("NEEDS_USER_DECISION: capture library ",library," has ",sum(selected)," cells; samples=",paste(samples,collapse=";"))
      }
      doubled <- qc_run_doublets(object,cfg,decision,paths,report,if(is.null(df_api)) qc_default_doublet_api() else df_api)
      final <- qc_score_cycle(doubled$object,paths,report)
      status <- if(doubled$unevaluated>0L) "QC_COMPLETE_WITH_UNEVALUATED_DOUBLETS" else "COMPLETE_QC"
    } else {
      qc_stage_validate(object,"qc2",gse,legacy=is.null(object@misc$qc_pipeline$qc2))
      final <- qc_run_ambient(object,paths,report,ambient_api)
      status <- "COMPLETE_QC3"
    }
    if(is.null(final@misc$qc_pipeline)) final@misc$qc_pipeline <- object@misc$qc_pipeline
    qc_stage_save(final,paths,step,input,cfg,report,status,decision)
  },error=function(e) {
    qc_stage_log(paths,step,paste("FAILED",conditionMessage(e)))
    qc_report_add(report,"RUN",metric="status",value="FAILED",status="FAILED",note=conditionMessage(e))
    qc_report_write(report,paths$report)
    stop(e)
  })
}

qc_run <- function(root,gse,package_available=function(p) requireNamespace(p,quietly=TRUE),df_api=NULL,step="all",input=NULL,ambient_api=NULL) {
  if(!step %in% c("all","qc1","qc2","qc3")) stop("Unknown QC step")
  base <- qc_paths(root,gse)
  if(step=="all" && !is.null(input)) stop("Use --input with one explicit QC step")
  steps <- if(step=="all") c("qc1","qc2","qc3") else step
  for(current in steps) {
    if(current=="qc2" && is.null(input) && !file.exists(qc_stage_paths(base,"qc1")$final))
      qc_run_stage(root,gse,"qc1",package_available,df_api)
    object <- qc_run_stage(root,gse,current,package_available,df_api,ambient_api,input)
  }
  records <- lapply(base$stage_dirs,function(directory) {
    file <- file.path(directory,"qc_report.csv")
    if(file.exists(file)) utils::read.csv(file,stringsAsFactors=FALSE,check.names=FALSE) else NULL
  })
  combined <- qc_report_new(); combined$rows <- Filter(Negate(is.null),records)
  qc_report_write(combined,base$report)
  if(tail(steps,1)=="qc3") "COMPLETE_QC3" else if(tail(steps,1)=="qc1") "COMPLETE_QC1" else
    if(any(object$DF_adj=="NotEvaluated")) "QC_COMPLETE_WITH_UNEVALUATED_DOUBLETS" else "COMPLETE_QC"
}

if(sys.nframe()==0L) {
  args <- commandArgs(trailingOnly=TRUE)
  if(!length(args) %in% 2:4) stop("Usage: run_r45.py qc/run_qc.R ROOT GSE [qc1|qc2|qc3|all] [INPUT_RDS]")
  cat("QC run:",qc_run(args[1],args[2],step=if(length(args)>2L) args[3] else "all",input=if(length(args)>3L) args[4] else NULL),"\n")
}
