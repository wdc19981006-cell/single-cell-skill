#!/usr/bin/env Rscript
qc_precheck_source <- if(!is.null(sys.frame(1)$ofile)) sys.frame(1)$ofile else sub("^--file=","",commandArgs()[grepl("^--file=",commandArgs())][1])
source(file.path(dirname(normalizePath(qc_precheck_source,winslash="/",mustWork=TRUE)),"qc_utils.R"))

qc_precheck <- function(root,gse,package_available=function(p) requireNamespace(p,quietly=TRUE)) {
  paths <- qc_paths(root,gse)
  report <- qc_report_new()
  started <- proc.time()[["elapsed"]]
  finish <- function(status,note="") {
    qc_report_add(report,"RUN",metric="status",value=status,status=if(status=="PASS") "PASS" else status,note=note)
    qc_report_add(report,"RUN",metric="precheck_seconds",value=round(proc.time()[["elapsed"]]-started,3))
    qc_report_write(report,paths$report)
    status
  }
  if (!file.exists(paths$raw)) return(finish("FAILED_INPUT","seurat_raw.rds missing"))
  missing <- qc_check_dependencies(report,package_available)
  if (length(missing)) return(finish("FAILED_INPUT","Required package missing; no packages were installed"))
  object <- tryCatch(readRDS(paths$raw),error=function(e)e)
  if (inherits(object,"error")) return(finish("FAILED_INPUT",conditionMessage(object)))
  result <- tryCatch(qc_validate_input(object,gse),error=function(e)e)
  if (inherits(result,"error")) return(finish("FAILED_INPUT",conditionMessage(result)))
  samples <- unique(as.character(object$sample))
  decision <- tryCatch(qc_decision(paths,samples),error=function(e)e)
  if (inherits(decision,"error")) return(finish("FAILED_INPUT",conditionMessage(decision)))
  cfg <- tryCatch(qc_config(paths,decision),error=function(e)e)
  if (inherits(cfg,"error")) return(finish("FAILED_INPUT",conditionMessage(cfg)))
  qc_report_config(report,cfg)

  needs <- character()
  if (!identical(cfg$species,"human") && !identical(cfg$species,"mouse")) needs <- c(needs,"Unsupported species")
  source_text <- paste(capture.output(print(unique(object@meta.data[,intersect(names(object@meta.data),c("source_type","modality","technology")),drop=FALSE]))),collapse=" ")
  if (grepl("snRNA|single[- ]nucle|nuclei",source_text,ignore.case=TRUE) &&
      is.null(decision$accept_default_conflicts) && is.null(decision$parameter_overrides)) {
    qc_report_add(report,"CELL_QC",metric="modality",value="snRNA-seq",status="NEEDS_USER_DECISION",note="Confirm species and QC thresholds for nuclei")
    needs <- c(needs,"snRNA-seq thresholds require confirmation")
  }
  patterns <- qc_patterns(cfg$species)
  features <- rownames(result)
  for (metric in names(patterns)) {
    hits <- sum(grepl(patterns[[metric]],features))
    qc_report_add(report,"CELL_QC",metric=paste0(metric,"_matched_features"),value=hits,
                  status=if(hits==0L) "NEEDS_USER_DECISION" else "PASS",
                  note=if(hits==0L) "Review species or gene naming; percentage would be zero" else "")
    if (hits==0L) needs <- c(needs,paste(metric,"features absent"))
  }
  object <- qc_add_percentages(object,cfg$species)
  metadata <- object@meta.data
  removed <- decision$removed_samples %||% character()
  kept <- decision$kept_low_cell_samples %||% character()
  for (sample in samples) {
    cells <- sum(metadata$sample==sample)
    group <- unique(as.character(metadata$group[metadata$sample==sample]))
    if (length(group)!=1L) return(finish("FAILED_INPUT",paste("Sample has mixed group:",sample)))
    qc_report_add(report,"SAMPLE",sample,"group",group,status="PASS")
    qc_report_add(report,"SAMPLE",sample,"input_cells",cells,status="PASS")
    if (sample %in% removed) {
      qc_report_add(report,"SAMPLE",sample,"decision","removed",status="USER_CONFIRMED_REMOVAL",note=decision$user_statement)
      next
    }
    if (cells < cfg$doublet$min_cells) {
      if (sample %in% kept) {
        qc_report_add(report,"SAMPLE",sample,"doublet_decision","NotEvaluated",status="USER_CONFIRMED_SKIP",
                      note=paste("group=",group,"; cells=",cells,"; ",decision$user_statement,sep=""))
      } else {
        qc_report_add(report,"SAMPLE",sample,"cell_count",cells,status="LOW_CELL_SAMPLE",
                      note=paste("group=",group,"; min_cells_for_doublet=",cfg$doublet$min_cells,sep=""))
        needs <- c(needs,paste0(sample," (",cells," cells; group=",group,")"))
      }
    }
  }
  active <- !metadata$sample %in% removed
  if (!any(active)) return(finish("FAILED_INPUT","All samples would be removed"))
  rules <- list(nFeature_min=metadata$nFeature_RNA < cfg$cell_qc$nFeature_min,
                nFeature_max=metadata$nFeature_RNA > cfg$cell_qc$nFeature_max,
                nCount_min=metadata$nCount_RNA < cfg$cell_qc$nCount_min,
                nCount_max=metadata$nCount_RNA > cfg$cell_qc$nCount_max,
                pMT_max=metadata$pMT > cfg$cell_qc$pMT_max,
                pHB_max=metadata$pHB > cfg$cell_qc$pHB_max,
                pRP_max=metadata$pRP > cfg$cell_qc$pRP_max)
  for (metric in names(rules)) {
    fraction <- mean(rules[[metric]][active])
    conflict <- fraction > 0.5 && !isTRUE(decision$accept_default_conflicts) &&
      is.null((decision$parameter_overrides %||% list())$cell_qc[[metric]])
    qc_report_add(report,"CELL_QC",metric=paste0(metric,"_fraction_outside"),value=round(fraction,4),
      status=if(conflict) "NEEDS_USER_DECISION" else "PASS",note=if(conflict) "More than half of active cells fail this default threshold" else "")
    if(conflict) needs <- c(needs,paste("Default",metric,"conflicts with data distribution"))
  }
  predicted_keep <- active & !Reduce(`|`,rules)
  for (sample in setdiff(samples,removed)) {
    count <- sum(predicted_keep & metadata$sample==sample)
    qc_report_add(report,"CELL_QC",sample,"predicted_after_qc_cells",count)
    if(count==0L) {
      qc_report_add(report,"CELL_QC",sample,"predicted_after_qc_cells",count,status="NEEDS_USER_DECISION",note="Default QC would empty this sample")
      needs <- c(needs,paste("Cell QC would empty",sample))
    } else if(count<cfg$doublet$min_cells && !sample %in% kept) {
      qc_report_add(report,"SAMPLE",sample,"predicted_after_qc_cells",count,status="LOW_CELL_SAMPLE",
                    note="DoubletFinder minimum would not be met after cell QC")
      needs <- c(needs,paste0(sample," (",count," cells after QC; group=",unique(metadata$group[metadata$sample==sample]),")"))
    }
  }
  if(length(needs)) {
    low <- needs[grepl("cells; group=",needs,fixed=TRUE)]
    if(length(low)) cat("以下样本细胞数不足",cfg$doublet$min_cells,"，不建议单独运行DoubletFinder：\n",paste(low,collapse="\n"),"\n请用户决定：删除样本；保留但跳过DoubletFinder；或修改最低细胞数阈值。\n",sep="")
    return(finish("NEEDS_USER_DECISION",paste(unique(needs),collapse="; ")))
  }
  finish("PASS")
}

if (sys.nframe()==0L) {
  args <- commandArgs(trailingOnly=TRUE)
  if(length(args)!=2L) stop("Usage: Rscript qc/qc_precheck.R REPOSITORY_ROOT GSE")
  status <- qc_precheck(args[1],args[2])
  cat("QC precheck:",status,"\n")
  if(status!="PASS") quit(status=1L)
}
