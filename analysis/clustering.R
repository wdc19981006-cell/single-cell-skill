analysis_ari <- function(a,b) {
  if(length(a)!=length(b) || length(a)<2L) stop("ARI requires aligned clusterings with at least two cells")
  pairs <- function(x) x*(x-1)/2
  cross <- table(a,b)
  observed <- sum(pairs(cross))
  row_pairs <- sum(pairs(rowSums(cross)))
  col_pairs <- sum(pairs(colSums(cross)))
  expected <- row_pairs*col_pairs/pairs(sum(cross))
  denominator <- (row_pairs+col_pairs)/2-expected
  if(abs(denominator)<.Machine$double.eps^0.5)
    return(if(identical(outer(a,a,"=="),outer(b,b,"=="))) 1 else 0)
  (observed-expected)/denominator
}

analysis_clustree_plot <- function(frame) {
  # Explicit guide objects avoid ggraph's legacy string registration path,
  # which ggplot2 4 rejects as "Unknown guide: edge_colourbar".
  clustree::clustree(frame,prefix="cluster_res.") +
    ggplot2::guides(edge_colour=ggraph::guide_edge_colourbar(),edge_fill=ggraph::guide_edge_colourbar()) +
    theme_sc()
}

analysis_cluster_associations <- function(meta,columns,resolutions) {
  variables <- intersect(c("pMT","nCount_RNA","nFeature_RNA","S.Score","G2M.Score"),names(meta))
  do.call(rbind,lapply(seq_along(columns),function(i) do.call(rbind,lapply(variables,function(variable) {
    x <- meta[[variable]]
    valid <- is.finite(x)
    total <- sum((x[valid]-mean(x[valid]))^2)
    groups <- split(x[valid],meta[[columns[i]]][valid])
    between <- sum(vapply(groups,function(group) length(group)*(mean(group)-mean(x[valid]))^2,numeric(1)))
    data.frame(resolution=resolutions[i],variable=variable,n=sum(valid),
      eta_squared=if(total>0) between/total else NA_real_,
      note="Descriptive between-cluster variance fraction; biological review required before changing regression")
  }))))
}

analysis_run_cluster_review <- function(input,out,cfg,provenance) {
  analysis_require(c("Seurat","SeuratObject","jsonlite","yaml","ggplot2","patchwork","clustree"))
  object <- readRDS(input)
  analysis_validate_object(object,"pca")
  reduction <- cfg$cluster_reduction
  if(reduction=="auto") reduction <- if(isTRUE(object@misc$workflow_analysis$harmony_completed) &&
                                         "harmony" %in% names(object@reductions)) "harmony" else "pca"
  if(!reduction %in% names(object@reductions)) stop("Selected representation is missing: ",reduction)
  dims <- seq_len(object@misc$workflow_analysis$pc_selection$pc_used)
  object <- Seurat::FindNeighbors(object,reduction=reduction,dims=dims,
    graph.name=c("analysis_nn","analysis_snn"),verbose=FALSE)
  umap_name <- if(reduction=="harmony") "umap_harmony" else "umap_unintegrated"
  if(!umap_name %in% names(object@reductions)) object <- analysis_umap(object,reduction,dims,umap_name,cfg$seed)
  # Every resolution uses this exact saved coordinate matrix, including finalization.
  object[["umap"]] <- object[[umap_name]]
  columns <- paste0("cluster_res.",format(cfg$resolutions,trim=TRUE,scientific=FALSE))
  if(anyDuplicated(columns)) stop("Resolution formatting produced duplicate metadata columns")
  summaries <- list()
  composition <- list()
  metric_rows <- list()
  for(i in seq_along(cfg$resolutions)) {
    resolution <- cfg$resolutions[i]
    object <- Seurat::FindClusters(object,graph.name="analysis_snn",resolution=resolution,
      random.seed=cfg$seed,verbose=FALSE)
    object[[columns[i]]] <- as.character(SeuratObject::Idents(object))
    sizes <- table(object@meta.data[[columns[i]]])
    summaries[[i]] <- data.frame(resolution=resolution,cluster=names(sizes),cells=as.integer(sizes),
      fraction=as.numeric(sizes)/ncol(object),representation=reduction)
    fields <- data.frame(cluster=object@meta.data[[columns[i]]],sample=as.character(object$sample))
    if("group" %in% names(object@meta.data)) fields$group <- as.character(object$group)
    comp <- stats::aggregate(rep(1L,ncol(object)),fields,sum)
    names(comp)[ncol(comp)] <- "cells"
    comp$resolution <- resolution
    composition[[i]] <- comp
    metric_rows[[i]] <- data.frame(resolution=resolution,n_clusters=length(sizes),
      min_cluster_cells=min(sizes),max_cluster_fraction=max(sizes)/ncol(object),representation=reduction)
  }
  # FindClusters assigns identities as a side effect. These are review candidates;
  # no formal labels or seurat_analysis.rds exist before a real user selection.
  generated <- names(object@meta.data)[grepl("^analysis_snn_res\\.|^seurat_clusters$",names(object@meta.data))]
  for(field in generated) object[[field]] <- NULL
  SeuratObject::Idents(object) <- object$sample
  object@misc$workflow_analysis$cluster_review <- list(resolutions=cfg$resolutions,columns=columns,
    reduction=reduction,umap=umap_name,selected=FALSE,provenance=provenance)
  panels <- lapply(seq_along(columns),function(i)
    analysis_embedding_panel(object,columns[i],"umap",paste("Resolution",cfg$resolutions[i])))
  # One PDF contains every page of the same-coordinate review, with at most 9 panels/page.
  size <- sc_plot_dimensions("clustering_grid")
  grDevices::pdf(file.path(out,"cluster_resolution_review.pdf"),width=size[1],height=size[2])
  tryCatch(for(start in seq(1L,length(panels),by=9L))
    print(patchwork::wrap_plots(panels[start:min(length(panels),start+8L)],ncol=3)),
    finally=grDevices::dev.off())
  for(page in seq_len(ceiling(length(panels)/9L))) {
    start <- (page-1L)*9L+1L
    grid <- patchwork::wrap_plots(panels[start:min(length(panels),start+8L)],ncol=3)
    ggplot2::ggsave(file.path(out,paste0("cluster_resolution_review_",page,".png")),grid,
      width=size[1],height=size[2],units="in",dpi=300,limitsize=FALSE,bg="white")
  }
  review_frame <- object@meta.data[,columns,drop=FALSE]
  for(column in columns) review_frame[[column]] <- as.numeric(review_frame[[column]])
  tree <- analysis_clustree_plot(review_frame)
  analysis_plot(tree,"clustree","violin",out)
  analysis_write_csv(do.call(rbind,summaries),file.path(out,"cluster_summary.csv"))
  analysis_write_csv(do.call(rbind,composition),file.path(out,"cluster_composition.csv"))
  analysis_write_csv(do.call(rbind,metric_rows),file.path(out,"cluster_metrics.csv"))
  analysis_write_csv(data.frame(cell=colnames(object),object@meta.data[,c("sample",columns),drop=FALSE],check.names=FALSE),
    file.path(out,"cluster_assignments.csv"))
  analysis_write_csv(analysis_cluster_associations(object@meta.data,columns,cfg$resolutions),
    file.path(out,"cluster_regression_diagnostics.csv"))
  stability <- do.call(rbind,lapply(seq_len(length(columns)-1L),function(i) data.frame(
    resolution_from=cfg$resolutions[i],resolution_to=cfg$resolutions[i+1L],
    adjusted_rand_index=analysis_ari(object@meta.data[[columns[i]]],object@meta.data[[columns[i+1L]]]),
    note="Stability between adjacent resolution candidates on the same graph; does not select a final resolution")))
  analysis_write_csv(stability,file.path(out,"resolution_stability.csv"))
  yaml::write_yaml(list(resolutions=cfg$resolutions,reduction=reduction,pc_used=length(dims),
    seed=cfg$seed,formal_resolution=NULL,decision="USER_SELECTION_REQUIRED"),file.path(out,"clustering_config.yaml"))
  analysis_save(object,file.path(out,"seurat_cluster_review.rds"))
  warnings <- analysis_warning(analysis_warning_table(),"BIOLOGICAL_TRADEOFF",
    "Cluster resolution candidates require user review; no formal cluster labels have been assigned",
    "Review marker structure, cluster size, clustree and stability, then provide a resolution decision",severity="ACTION_REQUIRED")
  summary <- data.frame(metric=c("cells","representation","pc_used","candidate_resolutions","status"),
    value=c(ncol(object),reduction,length(dims),paste(cfg$resolutions,collapse=";"),"NEEDS_RESOLUTION_SELECTION"))
  analysis_finish("cluster-review","NEEDS_RESOLUTION_SELECTION",out,summary,warnings)
}

analysis_finalize <- function(input,out,cfg,provenance,decision) {
  analysis_require(c("Seurat","SeuratObject","jsonlite","yaml"))
  object <- readRDS(input)
  analysis_validate_object(object,"cluster-review")
  if(!is.list(decision) || !identical(decision$confirmed_by,"user") ||
     !is.character(decision$user_statement) || length(decision$user_statement)!=1L ||
     !nzchar(trimws(decision$user_statement)) || !is.numeric(decision$resolution) ||
     length(decision$resolution)!=1L || !is.finite(decision$resolution))
    stop("Final labels require a literal user selection and confirmed_by=user")
  review <- object@misc$workflow_analysis$cluster_review
  index <- which(abs(review$resolutions-decision$resolution)<1e-10)
  if(length(index)!=1L) stop("Chosen resolution was not present in the reviewed scan")
  labels <- object@meta.data[[review$columns[index]]]
  object$seurat_clusters <- factor(labels)
  SeuratObject::Idents(object) <- object$seurat_clusters
  object@misc$workflow_analysis$cluster_review$selected <- TRUE
  object@misc$workflow_analysis$resolution_decision <- decision
  object@misc$workflow_analysis$final_provenance <- provenance
  yaml::write_yaml(decision,file.path(out,"selected_resolution.yaml"))
  analysis_save(object,file.path(out,"seurat_analysis.rds"))
  summary <- data.frame(metric=c("cells","clusters","selected_resolution","confirmed_by"),
    value=c(ncol(object),length(unique(labels)),decision$resolution,"user"))
  analysis_finish("finalize","COMPLETE",out,summary,analysis_warning_table())
}
