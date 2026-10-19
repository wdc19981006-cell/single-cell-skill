# Fixed, targeted post-clustering comparisons. Nothing here replaces the baseline.
analysis_marker_structure <- function(object,labels,cfg,scheme) {
  data <- SeuratObject::LayerData(object,assay=cfg$assay,layer="data")
  markers <- unique(c(cfg$markers,cfg$proliferation_markers))
  do.call(rbind,lapply(markers,function(marker) {
    if(!marker %in% rownames(data)) return(data.frame(scheme=scheme,cluster=NA_character_,marker=marker,
      mean_expression=NA_real_,fraction_detected=NA_real_,status="MARKER_ABSENT"))
    expression <- as.numeric(data[marker,])
    do.call(rbind,lapply(sort(unique(labels)),function(label) {
      x <- expression[labels==label]
      data.frame(scheme=scheme,cluster=label,marker=marker,mean_expression=mean(x),
        fraction_detected=mean(x>0),status="MEASURED")
    }))
  }))
}

analysis_composition <- function(object,labels,scheme) {
  fields <- intersect(c("sample","group","Phase"),names(object@meta.data))
  do.call(rbind,lapply(fields,function(field) {
    x <- as.data.frame(table(cluster=labels,value=as.character(object@meta.data[[field]])))
    names(x)[3] <- "cells"; x$field <- field; x$scheme <- scheme; x
  }))
}

analysis_comparison_recommendation <- function(ari,overlap,before,after,cfg) {
  if(ari>=cfg$regression_review$similar_ari && overlap>=cfg$regression_review$similar_subspace)
    return("SIMILAR_RESULTS")
  if(is.finite(before) && is.finite(after) && after>=before-cfg$regression_review$association_tolerance)
    return("BASELINE_PREFERRED")
  "USER_REVIEW_REQUIRED"
}

analysis_run_regression_review <- function(input,out,cfg,provenance) {
  baseline <- readRDS(input)
  analysis_validate_object(baseline,"cluster-review")
  info <- baseline@misc$workflow_analysis
  if(!identical(cfg$assay,info$assay)) stop("Comparison assay must equal baseline assay")
  cfg$seed <- info$seed
  cfg$max_pcs <- info$pc_selection$max_pc
  diagnostics <- analysis_regression_diagnostics(baseline,cfg)
  review <- info$cluster_review
  associations <- analysis_cluster_associations(baseline@meta.data,review$columns,review$resolutions)
  analysis_write_csv(diagnostics,file.path(out,"baseline_pc_associations.csv"))
  analysis_write_csv(associations,file.path(out,"baseline_cluster_associations.csv"))
  resolution <- cfg$regression_review$comparison_resolution
  index <- which(abs(review$resolutions-resolution)<1e-10)
  # A comparison setting must be an actually computed baseline resolution.
  if(length(index)!=1L) stop("comparison_resolution must be present in the baseline resolution scan")
  labels <- as.character(baseline@meta.data[[review$columns[index]]])
  analysis_write_csv(analysis_cluster_metric_profiles(baseline@meta.data,labels),file.path(out,"baseline_cluster_profiles.csv"))
  baseline_markers <- analysis_marker_structure(baseline,labels,cfg,"baseline")
  if(!is.null(baseline_markers)) analysis_write_csv(baseline_markers,file.path(out,"baseline_marker_review.csv"))
  cycle_risk <- any(diagnostics$risk & diagnostics$variable %in% c("S.Score","G2M.Score")) ||
    any(associations$variable %in% c("S.Score","G2M.Score") & associations$eta_squared>=cfg$regression_review$cluster_eta_warning,na.rm=TRUE)
  cycle_review <- cfg$regression_review$cell_cycle_review
  cycle_approved <- is.list(cycle_review) && identical(cycle_review$confirmed_by,"user") &&
    is.character(cycle_review$user_statement) && nzchar(trimws(cycle_review$user_statement)) &&
    is.character(cycle_review$marker_interpretation) && nzchar(trimws(cycle_review$marker_interpretation))
  rows <- list(); warnings <- analysis_warning_table()
  schemes <- list(feature="nFeature_RNA",count="nCount_RNA",cell_cycle=c("S.Score","G2M.Score"))
  for(target in names(schemes)) {
    extra <- schemes[[target]]
    pc_risk <- any(diagnostics$risk & diagnostics$variable %in% extra)
    cluster_risk <- any(associations$variable %in% extra &
      associations$eta_squared>=cfg$regression_review$cluster_eta_warning,na.rm=TRUE)
    risk <- pc_risk || cluster_risk
    requested <- target %in% cfg$sensitivity || (target!="cell_cycle" && isTRUE(cfg$regression_review$auto_numeric) && risk)
    status <- if(all(extra %in% info$regress)) "SKIPPED_ALREADY_REGRESSED" else
      if(!all(extra %in% names(baseline@meta.data))) "SKIPPED_MISSING_METADATA" else
      if(!risk) "SKIPPED_NO_CLEAR_RISK" else if(target=="cell_cycle" && !cycle_approved) "MARKER_REVIEW_REQUIRED" else
      if(!requested && !(target=="cell_cycle" && cycle_approved)) "SKIPPED_NOT_REQUESTED" else "COMPLETE"
    if(status=="MARKER_REVIEW_REQUIRED") warnings <- analysis_warning(warnings,"CELL_CYCLE_BIOLOGICAL_REVIEW",
      "Cycle scores distinguish PCs/clusters; inspect proliferation markers before a cycle regression comparison",
      "Read baseline_marker_review.csv; record marker interpretation and the user's choice if a comparison is needed",severity="ACTION_REQUIRED")
    row <- data.frame(target=target,status=status,variables=paste(c("pMT",extra),collapse=";"),
      pc_risk=pc_risk,cluster_risk=cluster_risk,resolution=resolution,ari=NA_real_,subspace_overlap=NA_real_,
      recommendation=if(status=="MARKER_REVIEW_REQUIRED") "USER_REVIEW_REQUIRED" else "BASELINE_PREFERRED")
    if(status!="COMPLETE") { rows[[target]] <- row; next }
    directory <- file.path(out,target); dir.create(directory,recursive=TRUE,showWarnings=FALSE)
    log <- file.path(directory,"comparison.log")
    cat("START targeted comparison from baseline QC input; no formal output replacement\n",file=log)
    tryCatch({
      qc_input <- info$provenance$input
      if(is.null(qc_input) || !file.exists(qc_input)) stop("Original QC input is unavailable")
      if(!identical(digest::digest(file=qc_input,algo="sha256"),info$provenance$input_sha256)) stop("Original QC input changed since baseline PCA")
      comparison <- analysis_pca(readRDS(qc_input),cfg,c("pMT",extra),features=info$features %||%
        SeuratObject::VariableFeatures(baseline[[cfg$assay]]))
      if(!identical(colnames(comparison),colnames(baseline)) ||
         !identical(analysis_counts(comparison,cfg$assay),analysis_counts(baseline,cfg$assay)))
        stop("Comparison cell set or expression differs from baseline")
      dims <- seq_len(info$pc_selection$pc_used)
      if(ncol(SeuratObject::Embeddings(comparison[["pca"]]))<length(dims)) stop("Candidate cannot support fixed baseline PC count")
      if(review$reduction=="harmony") {
        design <- info$harmony_design
        if(!isTRUE(design$eligible)) stop("Baseline Harmony design unavailable")
        comparison <- analysis_apply_harmony(comparison,design$selected,dims,cfg$seed)
      }
      comparison <- Seurat::FindNeighbors(comparison,reduction=review$reduction,dims=dims,
        graph.name=c("comparison_nn","comparison_snn"),verbose=FALSE)
      comparison <- Seurat::FindClusters(comparison,graph.name="comparison_snn",resolution=resolution,
        random.seed=cfg$seed,verbose=FALSE)
      candidate <- as.character(SeuratObject::Idents(comparison))
      comparison$comparison_clusters <- candidate
      comparison@misc$workflow_analysis$comparison <- list(baseline_input=input,resolution=resolution,
        fixed_pc_used=length(dims),reduction=review$reduction,harmony_design=info$harmony_design,
        formal_use="comparison only; adoption requires a new isolated run and PC/resolution review")
      a <- SeuratObject::Embeddings(baseline[["pca"]])[,dims,drop=FALSE]
      b <- SeuratObject::Embeddings(comparison[["pca"]])[,dims,drop=FALSE]
      row$subspace_overlap <- sum(svd(crossprod(qr.Q(qr(a)),qr.Q(qr(b))),nu=0,nv=0)$d^2)/length(dims)
      row$ari <- analysis_ari(labels,candidate)
      before <- analysis_cluster_associations(data.frame(baseline@meta.data,comparison_clusters=labels),"comparison_clusters",resolution)
      after <- analysis_cluster_associations(comparison@meta.data,"comparison_clusters",resolution)
      max_association <- function(x) { y <- x$eta_squared[x$variable %in% extra]; if(any(is.finite(y))) max(y,na.rm=TRUE) else NA_real_ }
      row$recommendation <- analysis_comparison_recommendation(row$ari,row$subspace_overlap,max_association(before),max_association(after),cfg)
      before$scheme <- "baseline"; after$scheme <- target
      analysis_write_csv(rbind(before,after),file.path(directory,"cluster_associations.csv"))
      profiles <- rbind(data.frame(scheme="baseline",analysis_cluster_metric_profiles(baseline@meta.data,labels)),
        data.frame(scheme=target,analysis_cluster_metric_profiles(comparison@meta.data,candidate)))
      analysis_write_csv(profiles,file.path(directory,"cluster_metric_profiles.csv"))
      transitions <- as.data.frame(table(baseline=labels,candidate=candidate)); names(transitions)[3] <- "cells"
      analysis_write_csv(transitions,file.path(directory,"cluster_transitions.csv"))
      sizes <- rbind(data.frame(scheme="baseline",cluster=names(table(labels)),cells=as.integer(table(labels))),
        data.frame(scheme=target,cluster=names(table(candidate)),cells=as.integer(table(candidate))))
      analysis_write_csv(sizes,file.path(directory,"cluster_sizes.csv"))
      analysis_write_csv(data.frame(cell=colnames(baseline),baseline=labels,candidate=candidate),file.path(directory,"cell_cluster_changes.csv"))
      composition <- rbind(analysis_composition(baseline,labels,"baseline"),analysis_composition(comparison,candidate,target))
      analysis_write_csv(composition,file.path(directory,"sample_group_cycle_composition.csv"))
      marker <- rbind(baseline_markers,analysis_marker_structure(comparison,candidate,cfg,target))
      if(!is.null(marker)) analysis_write_csv(marker,file.path(directory,"marker_structure.csv"))
      before_neighbors <- analysis_neighbor_metrics(baseline,review$reduction,dims,"baseline",cfg$markers)
      after_neighbors <- analysis_neighbor_metrics(comparison,review$reduction,dims,target,cfg$markers)
      analysis_write_csv(rbind(before_neighbors$markers,after_neighbors$markers),file.path(directory,"marker_neighborhood.csv"))
      analysis_write_csv(rbind(before_neighbors$metrics,after_neighbors$metrics),file.path(directory,"sample_group_neighborhood.csv"))
      analysis_write_csv(analysis_regression_diagnostics(comparison,cfg),file.path(directory,"regression_diagnostics.csv"))
      analysis_write_csv(row,file.path(directory,"sensitivity_comparison.csv"))
      jsonlite::write_json(list(configuration=cfg,regress=c("pMT",extra),input=qc_input,input_sha256=info$provenance$input_sha256,
        baseline=input,features=info$features,pc_used=length(dims),resolution=resolution,
        seed=cfg$seed,harmony_design=info$harmony_design,metric_source=info$metric_source),
        file.path(directory,"comparison_config.json"),auto_unbox=TRUE,pretty=TRUE)
      analysis_plot(patchwork::wrap_plots(analysis_embedding_panel(baseline,"sample","pca","Baseline PCA"),
        analysis_embedding_panel(comparison,"sample","pca",paste(target,"PCA")),ncol=2),"pca_comparison","umap_2panel",directory)
      p <- ggplot2::ggplot(transitions,ggplot2::aes(x=baseline,y=candidate,fill=cells))+ggplot2::geom_tile()+theme_sc()+
        ggplot2::labs(title=paste(target,"cluster transitions; ARI",round(row$ari,3)))
      analysis_plot(p,"cluster_comparison","diagnostics",directory)
      if(!is.null(marker) && any(marker$status=="MEASURED")) {
        p <- ggplot2::ggplot(marker[marker$status=="MEASURED",],ggplot2::aes(x=cluster,y=marker,size=fraction_detected,color=mean_expression))+
          ggplot2::geom_point()+ggplot2::facet_wrap(~scheme,scales="free_x")+theme_sc()+ggplot2::labs(title="Marker structure; cluster labels are independent")
        analysis_plot(p,"marker_comparison","diagnostics",directory)
      }
      analysis_save(comparison,file.path(directory,"seurat_sensitivity.rds"))
      analysis_validate_object(readRDS(file.path(directory,"seurat_sensitivity.rds")),"sensitivity")
      cat("COMPLETE; diagnostic recommendation ",row$recommendation,"; baseline preserved\n",file=log,append=TRUE)
    },error=function(e) { cat("FAILED ",conditionMessage(e),"\n",file=log,append=TRUE); stop(e) })
    rows[[target]] <- row
  }
  summary <- do.call(rbind,rows)
  analysis_write_csv(summary,file.path(out,"sensitivity_summary.csv"))
  required <- any(summary$recommendation=="USER_REVIEW_REQUIRED")
  if(any(summary$status=="COMPLETE" & summary$recommendation=="USER_REVIEW_REQUIRED"))
    warnings <- analysis_warning(warnings,"REGRESSION_BIOLOGICAL_TRADEOFF","Targeted regression changes clustering/structure",
      "Review transitions, markers, phase and sample/group evidence; candidates never replace baseline",severity="ACTION_REQUIRED")
  record <- list(status=if(required) "USER_REVIEW_REQUIRED" else "COMPLETE",baseline_input=input,
    summary=summary,interpretation="Rule-based diagnostic advice; no biological correctness is established")
  saveRDS(record,file.path(out,"regression_review.rds"))
  analysis_finish("regression-review",record$status,out,summary,warnings,cfg,provenance)
}
