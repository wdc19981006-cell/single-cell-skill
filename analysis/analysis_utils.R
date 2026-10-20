if (!isTRUE(getOption("r45.runtime.active"))) stop("Use runtime/r45/run_r45.py")
analysis_dir <- dirname(normalizePath(sys.frame(1)$ofile,winslash="/",mustWork=TRUE))
analysis_repo <- dirname(analysis_dir)
source(file.path(analysis_repo,"runtime","seurat_storage.R"))
for (module in c("qc_utils.R","doublet_finder.R","plot_utils.R"))
  source(file.path(analysis_repo,"qc",module))

analysis_require <- function(packages) {
  for (package in packages) if(!requireNamespace(package,quietly=TRUE))
    stop("R_PACKAGE_MISSING: ",package,"; only this analysis route is stopped; no packages were changed")
}

analysis_config <- function(cfg) {
  scalar <- function(x,lower=0) is.numeric(x) && length(x)==1L && is.finite(x) && x>=lower
  for (field in c("seed","nfeatures","max_pcs"))
    if(!scalar(cfg[[field]],if(field=="seed") 0 else 3) || cfg[[field]]!=as.integer(cfg[[field]]))
      stop("Invalid analysis configuration: ",field)
  if(!scalar(cfg$risk_correlation) || cfg$risk_correlation>1) stop("Invalid risk_correlation")
  for (field in c("pc_cumulative_variance","pc_individual_variance","pc_elbow_drop"))
    if(!scalar(cfg[[field]]) || cfg[[field]]>100) stop("Invalid PCA selection threshold: ",field)
  if(!is.character(cfg$assay) || length(cfg$assay)!=1L || !nzchar(cfg$assay)) stop("Select one expression assay")
  cfg$sensitivity <- unlist(cfg$sensitivity,use.names=FALSE)
  cfg$sensitivity[cfg$sensitivity=="rna_total"] <- "count"
  if(any(!cfg$sensitivity %in% c("feature","count","cell_cycle"))) stop("Unknown targeted sensitivity scheme")
  for(field in c("cluster_eta_warning","similar_ari","similar_subspace","association_tolerance"))
    if(!scalar(cfg$regression_review[[field]]) || cfg$regression_review[[field]]>1) stop("Invalid regression review threshold: ",field)
  if(!is.logical(cfg$regression_review$auto_numeric) || length(cfg$regression_review$auto_numeric)!=1L ||
     !scalar(cfg$regression_review$comparison_resolution,.Machine$double.eps)) stop("Invalid regression review settings")
  cfg$resolutions <- as.numeric(unlist(cfg$resolutions,use.names=FALSE))
  if(length(cfg$resolutions)<2L || any(!is.finite(cfg$resolutions) | cfg$resolutions<=0) || anyDuplicated(cfg$resolutions))
    stop("At least two unique positive clustering resolutions required")
  cfg$resolutions <- sort(cfg$resolutions)
  if(!cfg$cluster_reduction %in% c("auto","pca","harmony")) stop("Unknown cluster_reduction")
  cfg$markers <- as.character(unlist(cfg$markers,use.names=FALSE))
  cfg$proliferation_markers <- as.character(unlist(cfg$proliferation_markers,use.names=FALSE))
  if(!is.null(cfg$regression_decision)) {
    d <- cfg$regression_decision
    if(!identical(d$confirmed_by,"user") || !d$scheme %in% c("feature","count","cell_cycle") ||
       !is.character(d$user_statement) || length(d$user_statement)!=1L || !nzchar(trimws(d$user_statement)))
      stop("A new formal regression requires the user's literal scheme decision")
  }
  if(!is.logical(cfg$harmony$enabled) || length(cfg$harmony$enabled)!=1L ||
     !is.logical(cfg$harmony$allow_sample) || length(cfg$harmony$allow_sample)!=1L)
    stop("Harmony enabled and allow_sample must be boolean")
  cfg$harmony$technical_fields <- as.character(unlist(cfg$harmony$technical_fields,use.names=FALSE))
  if("group" %in% c(cfg$harmony$batch_field,cfg$harmony$technical_fields))
    stop("Biological group cannot be used as a Harmony correction variable")
  cfg
}

analysis_counts <- function(object,assay="RNA") {
  if(!inherits(object,"Seurat") || !assay %in% names(object@assays)) stop("Missing Seurat assay: ",assay)
  layers <- SeuratObject::Layers(object[[assay]])
  if(!"counts" %in% layers || any(startsWith(layers,"counts.")))
    stop("A single unambiguous counts layer is required for ",assay,"; normalized values cannot substitute")
  counts <- SeuratObject::LayerData(object,assay=assay,layer="counts")
  if(!inherits(counts,"sparseMatrix") || any(dim(counts)<2L) || is.null(rownames(counts)) ||
      is.null(colnames(counts)) || anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts)) ||
      any(qc_blank(rownames(counts))) || any(qc_blank(colnames(counts))) ||
      any(!is.finite(counts@x)) || any(counts@x<0) ||
      any(Matrix::colSums(counts)<=0)) stop("Invalid sparse nonnegative counts or empty cells in ",assay)
  if(assay=="RNA" && any(counts@x!=round(counts@x))) stop("RNA counts must contain original nonnegative integer counts")
  if(!identical(colnames(counts),rownames(object@meta.data))) stop("Counts/metadata cell order differs")
  counts
}

analysis_metadata <- function(object) {
  if(!"sample" %in% names(object@meta.data) || any(qc_blank(object$sample)))
    stop("A nonblank sample field is required")
  if(!"pMT" %in% names(object@meta.data) || !is.numeric(object$pMT) || any(!is.finite(object$pMT)))
    stop("Finite pMT required for the baseline regression")
  invisible(TRUE)
}

analysis_write_csv <- function(frame,path) {
  for(field in names(frame)) if(is.character(frame[[field]])) Encoding(frame[[field]]) <- "unknown"
  utils::write.csv(frame,path,row.names=FALSE,na="")
}

analysis_warning_table <- function() data.frame(category=character(),severity=character(),
  message=character(),recommendation=character(),stringsAsFactors=FALSE)

analysis_warning <- function(table,category,message,recommendation,severity="WARNING")
  rbind(table,data.frame(category=category,severity=severity,message=message,
                        recommendation=recommendation,stringsAsFactors=FALSE))

analysis_save <- function(object,path) {
  if(file.exists(path)) stop("Existing RDS preserved: ",path)
  if(basename(path)=="seurat_analysis.rds") object <- sc_compact_formal(object,path,"analysis")
  saveRDS(object,path)
  rebuilt <- readRDS(path)
  if(!identical(colnames(rebuilt),colnames(object)) ||
     !identical(analysis_counts(rebuilt,"RNA"),analysis_counts(object,"RNA")))
    stop("Serialized analysis cells or original RNA counts changed")
  assay <- object@misc$workflow_analysis$assay
  if(!is.null(assay) && (!identical(analysis_counts(rebuilt,assay),analysis_counts(object,assay)) ||
      !identical(rebuilt@meta.data,object@meta.data))) stop("Serialized selected assay or metadata changed")
  invisible(path)
}

analysis_plot <- function(plot,name,type,out) {
  # All dimensions, fonts, point sizing and output dpi come from the QC utility.
  save_sc_plot(plot,name,type,out)
}

analysis_embedding_panel <- function(object,field,reduction,title=field) {
  embedding <- SeuratObject::Embeddings(object[[reduction]])
  frame <- data.frame(x=embedding[,1],y=embedding[,2],
    value=object@meta.data[rownames(embedding),field],row.names=rownames(embedding))
  limits <- sc_umap_limits(embedding)
  ggplot2::ggplot(frame,ggplot2::aes(x=x,y=y,color=value)) +
    ggplot2::geom_point(size=get_umap_point_size(nrow(frame))) +
    ggplot2::labs(title=title,color=field,x=paste0(toupper(reduction),"_1"),y=paste0(toupper(reduction),"_2")) +
    ggplot2::coord_fixed(xlim=limits$x,ylim=limits$y,expand=FALSE) + theme_sc_umap()
}

analysis_umap <- function(object,reduction,dims,name,seed) {
  if(ncol(object)<5L) stop("At least five cells are required for an informative UMAP review")
  Seurat::RunUMAP(object,reduction=reduction,dims=dims,reduction.name=name,
    reduction.key=paste0(gsub("_","",name),"_"),n.neighbors=min(30L,ncol(object)-1L),
    seed.use=seed,verbose=FALSE)
}

analysis_pca <- function(object,cfg,regress="pMT",features=NULL) {
  original_rna <- analysis_counts(object,"RNA")
  original_meta <- object@meta.data
  selected_counts <- analysis_counts(object,cfg$assay)
  if(!identical(dimnames(selected_counts),dimnames(original_rna))) stop("Selected assay must preserve original RNA gene/cell IDs and order")
  analysis_metadata(object)
  SeuratObject::DefaultAssay(object) <- cfg$assay
  # Clear diagnostic reductions/graphs from DoubletFinder and any prior input;
  # the formal PCA is always recomputed from the explicitly selected counts.
  object@reductions <- list()
  object@graphs <- list()
  object@neighbors <- list()
  old <- names(object@meta.data)[grepl("_snn_res\\.|^seurat_clusters$|^cluster_res\\.",names(object@meta.data))]
  for(field in old) object[[field]] <- NULL
  SeuratObject::Idents(object) <- object$sample
  object <- Seurat::NormalizeData(object,assay=cfg$assay,verbose=FALSE)
  if(is.null(features)) {
    object <- Seurat::FindVariableFeatures(object,assay=cfg$assay,nfeatures=min(cfg$nfeatures,nrow(object[[cfg$assay]])-1L),verbose=FALSE)
    features <- SeuratObject::VariableFeatures(object,assay=cfg$assay)
  } else {
    if(anyDuplicated(features) || !all(features %in% rownames(object[[cfg$assay]]))) stop("Fixed HVGs invalid")
    SeuratObject::VariableFeatures(object[[cfg$assay]]) <- features
  }
  max_pc <- min(cfg$max_pcs,ncol(object)-1L,length(features)-1L,nrow(object)-1L)
  if(max_pc<2L) stop("Fewer than two usable PCA dimensions; retain input and review dataset size")
  if(any(!regress %in% names(object@meta.data)) ||
     any(!vapply(object@meta.data[,regress,drop=FALSE],function(x) is.numeric(x) && all(is.finite(x)),logical(1))))
    stop("Regression variables must have finite numeric metadata")
  object <- Seurat::ScaleData(object,assay=cfg$assay,features=features,vars.to.regress=regress,verbose=FALSE)
  object <- Seurat::RunPCA(object,assay=cfg$assay,features=features,npcs=max_pc,seed.use=cfg$seed,verbose=FALSE)
  if(!identical(analysis_counts(object,"RNA"),original_rna)) stop("Formal PCA changed original RNA counts")
  fields <- intersect(c("pMT","nCount_RNA","nFeature_RNA","S.Score","G2M.Score","Phase"),names(original_meta))
  if(!identical(object@meta.data[,fields,drop=FALSE],original_meta[,fields,drop=FALSE])) stop("Formal analysis changed QC metadata")
  selection <- qc_pc_selection(object[["pca"]]@stdev,cfg,max_pc)
  object@misc$workflow_analysis <- list(assay=cfg$assay,regress=regress,
    pc_selection=selection,formal_pca=TRUE,seed=cfg$seed,features=features,
    regression_decision=cfg$regression_decision,
    metric_source="QC metadata from original RNA; expression from selected assay; corrected counts retain fractional values")
  object
}

analysis_regression_diagnostics <- function(object,cfg) {
  embedding <- SeuratObject::Embeddings(object[["pca"]])
  selected_pcs <- object@misc$workflow_analysis$pc_selection$pc_used
  variables <- c("pMT","nCount_RNA","nFeature_RNA","S.Score","G2M.Score")
  rows <- lapply(variables,function(variable) do.call(rbind,lapply(seq_len(ncol(embedding)),function(pc) {
    available <- variable %in% names(object@meta.data)
    x <- if(available) object@meta.data[[variable]] else rep(NA_real_,ncol(object))
    valid <- is.finite(x) & is.finite(embedding[,pc])
    estimable <- sum(valid)>=3L && stats::sd(x[valid])>0 && stats::sd(embedding[valid,pc])>0
    rho <- if(estimable) stats::cor(x[valid],embedding[valid,pc],method="spearman") else NA_real_
    data.frame(variable=variable,pc=pc,rho=rho,n=sum(valid),
      selected_pc=pc<=selected_pcs,
      risk=if(is.finite(rho)) pc<=selected_pcs && abs(rho)>=cfg$risk_correlation else FALSE,
      status=if(!available) "MISSING_METADATA" else if(!estimable) "NOT_ESTIMABLE" else "MEASURED")
  })))
  do.call(rbind,rows)
}

analysis_pc_selection_table <- function(selection) {
  fields <- c("max_pc","pc_used","variance_pick","drop_pick","selection_method",
              "selected_pc_cumulative_sd_share","selected_pc_cumulative_variance")
  data.frame(metric=fields,value=vapply(selection[fields],function(x) as.character(x),character(1)),
    note="Selection uses standard-deviation share; variance explained is diagnostic only")
}

analysis_finish <- function(step,status,out,summary,warnings,cfg=NULL,provenance=NULL,details=NULL) {
  analysis_write_csv(summary,file.path(out,paste0(gsub("-","_",step),"_summary.csv")))
  analysis_write_csv(warnings,file.path(out,paste0(gsub("-","_",step),"_warnings.csv")))
  jsonlite::write_json(list(step=step,status=status,input=provenance$input,
    outputs=list.files(out,full.names=TRUE),parameters=cfg,statistics=summary,unresolved_issues=warnings,
    needs_user_decision=status %in% c("NEEDS_RESOLUTION_SELECTION","USER_REVIEW_REQUIRED"),details=details,
    next_command=switch(step,pca="analysis/run_workflow.py DATASET --step harmony",
      harmony="analysis/run_workflow.py DATASET --step cluster-review",
      `cluster-review`="analysis/run_workflow.py DATASET --step regression-review",
      `regression-review`="Review comparisons then finalize with the user's resolution-decision JSON",
      finalize="workflow_audit/run_audit.py DATASET")),file.path(out,paste0(gsub("-","_",step),"_result.json")),
                        auto_unbox=TRUE,pretty=TRUE)
}
