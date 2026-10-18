analysis_validate_object <- function(object,kind="pca") {
  analysis_counts(object,"RNA")
  analysis_metadata(object)
  info <- object@misc$workflow_analysis
  if(!isTRUE(info$formal_pca) || !identical(info$regress,"pMT"))
    stop("Formal baseline PCA with pMT-only regression is required")
  analysis_counts(object,info$assay)
  if(!"pca" %in% names(object@reductions)) stop("Missing formal PCA reduction")
  pca <- SeuratObject::Embeddings(object[["pca"]])
  used <- info$pc_selection$pc_used
  if(!identical(rownames(pca),colnames(object)) || ncol(pca)<2L || any(!is.finite(pca)) ||
     length(used)!=1L || !is.finite(used) || used<2L || used>ncol(pca)) stop("Invalid PCA embedding or selected PCs")
  if(kind=="harmony") {
    if(!isTRUE(info$harmony_completed) || !"harmony" %in% names(object@reductions))
      stop("Harmony stage has no completed Harmony reduction")
    corrected <- SeuratObject::Embeddings(object[["harmony"]])
    if(!identical(rownames(corrected),colnames(object)) || ncol(corrected)<used || any(!is.finite(corrected)))
      stop("Invalid Harmony embedding")
  }
  if(kind %in% c("cluster-review","finalize")) {
    review <- info$cluster_review
    if(!is.list(review) || !all(review$columns %in% names(object@meta.data)) || length(review$resolutions)<2L)
      stop("Incomplete cluster resolution review")
    if(!"umap" %in% names(object@reductions) || !review$umap %in% names(object@reductions))
      stop("Review UMAP is unavailable")
    actual <- SeuratObject::Embeddings(object[["umap"]])
    reference <- SeuratObject::Embeddings(object[[review$umap]])
    # Seurat assigns a new reduction key/column names when aliasing a DimReduc;
    # cell IDs and numeric coordinates must remain exactly identical.
    if(!identical(rownames(actual),rownames(reference)) || !identical(unname(actual),unname(reference)))
      stop("Review UMAP coordinates differ from selected representation")
    for(column in review$columns) if(any(qc_blank(object@meta.data[[column]]))) stop("Missing review cluster assignments")
    if(kind=="cluster-review" && (isTRUE(review$selected) || "seurat_clusters" %in% names(object@meta.data)))
      stop("Review object contains formal labels before user selection")
    if(kind=="finalize") {
      decision <- info$resolution_decision
      if(!isTRUE(review$selected) || !identical(decision$confirmed_by,"user") ||
         is.null(decision$user_statement) || !nzchar(trimws(decision$user_statement)) ||
         !"seurat_clusters" %in% names(object@meta.data)) stop("Formal result lacks a valid user decision")
      index <- which(abs(review$resolutions-decision$resolution)<1e-10)
      if(length(index)!=1L || !identical(as.character(object$seurat_clusters),object@meta.data[[review$columns[index]]]))
        stop("Formal labels differ from the user's reviewed resolution")
    }
  }
  invisible(TRUE)
}

if(sys.nframe()==0L) {
  file_arg <- commandArgs()[grepl("^--file=",commandArgs())][1]
  directory <- dirname(normalizePath(sub("^--file=","",file_arg),winslash="/",mustWork=TRUE))
  source(file.path(directory,"analysis_utils.R"))
  args <- commandArgs(trailingOnly=TRUE)
  if(length(args)!=2L) stop("Usage: validate_analysis.R ARTIFACT KIND")
  analysis_validate_object(readRDS(args[1]),args[2])
  cat("PASS: analysis RDS roundtrip, original counts, metadata, dimensions and decision contract\n")
}
