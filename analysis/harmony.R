analysis_design_check <- function(meta,cfg) {
  candidates <- unique(c("sample","batch","library_id","sequencing_run",cfg$harmony$technical_fields,
                         cfg$harmony$batch_field))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  selected <- cfg$harmony$batch_field
  trusted <- cfg$harmony$technical_fields
  if(is.null(selected)) {
    usable <- trusted[vapply(trusted,function(field) field %in% names(meta) &&
      !any(qc_blank(meta[[field]])) && length(unique(meta[[field]]))>1L,logical(1))]
    if(length(usable)) selected <- usable[1]
    else if(isTRUE(cfg$harmony$allow_sample)) selected <- "sample"
  }
  if(length(selected) && (length(selected)!=1L || identical(selected,"group")))
    stop("One nonbiological batch field is required")
  group_ok <- "group" %in% names(meta) && !any(qc_blank(meta$group))
  confounding <- function(field) {
    if(!group_ok || !field %in% names(meta) || any(qc_blank(meta[[field]])) ||
       length(unique(meta[[field]]))<2L || length(unique(meta$group))<2L) return(NA)
    frame <- data.frame(group=factor(meta$group),batch=factor(meta[[field]]))
    design <- stats::model.matrix(~ group + batch,data=frame)
    qr(design)$rank<ncol(design)
  }
  rows <- lapply(candidates,function(field) {
    available <- field %in% names(meta)
    data.frame(field=field,available=available,
      n_levels=if(available) length(unique(meta[[field]][!qc_blank(meta[[field]])])) else 0L,
      n_missing=if(available) sum(qc_blank(meta[[field]])) else nrow(meta),
      role=if(field=="sample") "BIOLOGICAL_SAMPLE" else
        if(field %in% trusted) "CONFIRMED_TECHNICAL" else "UNCONFIRMED_TECHNICAL_MEANING",
      confounded_with_group=confounding(field),selected=identical(field,selected))
  })
  table <- do.call(rbind,rows)
  acknowledgement <- cfg$harmony$confounding_acknowledgement
  acknowledged <- is.list(acknowledgement) && identical(acknowledgement$confirmed_by,"user") &&
    is.character(acknowledgement$user_statement) && length(acknowledgement$user_statement)==1L &&
    nzchar(trimws(acknowledgement$user_statement))
  reason <- if(!isTRUE(cfg$harmony$enabled)) "DISABLED_BY_CONFIGURATION" else
    if(is.null(selected)) "NO_CONFIRMED_TECHNICAL_BATCH" else
    if(!selected %in% names(meta)) "BATCH_FIELD_MISSING" else
    if(any(qc_blank(meta[[selected]]))) "BATCH_VALUES_MISSING" else
    if(length(unique(meta[[selected]]))<2L) "ONE_BATCH_LEVEL" else
    if(selected=="sample" && !isTRUE(cfg$harmony$allow_sample)) "SAMPLE_INTEGRATION_NOT_APPROVED" else
    if(selected!="sample" && !selected %in% trusted) "TECHNICAL_MEANING_UNCONFIRMED" else
    if(!group_ok && !acknowledged) "BIOLOGICAL_GROUP_UNAVAILABLE" else
    if(isTRUE(confounding(selected)) && !acknowledged) "BATCH_GROUP_CONFOUNDED" else "ELIGIBLE"
  table$decision <- reason
  list(table=table,selected=selected,eligible=identical(reason,"ELIGIBLE"),reason=reason,
       confounded=if(is.null(selected)) NA else confounding(selected),group_available=group_ok,
       acknowledgement=if(acknowledged) acknowledgement else NULL)
}

analysis_neighbors <- function(embedding,k=20L) {
  analysis_require("RANN")
  k <- min(as.integer(k),nrow(embedding)-1L)
  if(k<1L) stop("Neighbor diagnostics need at least two cells")
  # RANN includes the query cell; remove by ID rather than position (ties).
  found <- RANN::nn2(embedding,k=min(nrow(embedding),k+1L))$nn.idx
  t(vapply(seq_len(nrow(found)),function(i) head(found[i,found[i,]!=i],k),integer(k)))
}

analysis_neighbor_metrics <- function(object,reduction,dims,stage,markers,batch_field=NULL) {
  embedding <- SeuratObject::Embeddings(object[[reduction]])[,dims,drop=FALSE]
  neighbors <- analysis_neighbors(embedding)
  meta <- object@meta.data[rownames(embedding),,drop=FALSE]
  fields <- unique(intersect(c("sample","group","batch","library_id","sequencing_run",batch_field),names(meta)))
  rows <- lapply(fields,function(field) {
    values <- as.character(meta[[field]])
    valid <- !qc_blank(values)
    same <- rowMeans(matrix(values[neighbors],nrow=nrow(neighbors))==values)
    data.frame(representation=stage,metric="neighbor_same_label_fraction",field=field,
      value=mean(same[valid],na.rm=TRUE),cells=sum(valid),
      note="Descriptive kNN diagnostic; mixing alone does not establish integration quality")
  })
  marker_rows <- lapply(markers,function(marker) {
    available <- marker %in% rownames(object[[SeuratObject::DefaultAssay(object)]])
    if(!available) return(data.frame(representation=stage,marker=marker,status="MARKER_ABSENT",
      neighbor_spearman=NA_real_,mean_expression=NA_real_))
    expression <- as.numeric(SeuratObject::LayerData(object,assay=SeuratObject::DefaultAssay(object),layer="data")[marker,])
    smooth <- rowMeans(matrix(expression[neighbors],nrow=nrow(neighbors)))
    rho <- if(stats::sd(expression)>0 && stats::sd(smooth)>0) stats::cor(expression,smooth,method="spearman") else NA_real_
    data.frame(representation=stage,marker=marker,status=if(is.finite(rho)) "MEASURED" else "NOT_ESTIMABLE",
      neighbor_spearman=rho,mean_expression=mean(expression))
  })
  list(metrics=do.call(rbind,rows),markers=if(length(marker_rows)) do.call(rbind,marker_rows) else
    data.frame(representation=character(),marker=character(),status=character(),neighbor_spearman=double(),mean_expression=double()))
}

analysis_run_harmony <- function(input,out,cfg,provenance) {
  analysis_require(c("Seurat","SeuratObject","jsonlite","yaml","ggplot2","patchwork","RANN"))
  object <- readRDS(input)
  analysis_validate_object(object,"pca")
  original_rna <- analysis_counts(object,"RNA")
  warnings <- analysis_warning_table()
  design <- analysis_design_check(object@meta.data,cfg)
  analysis_write_csv(design$table,file.path(out,"batch_design_check.csv"))
  sample_fields <- intersect(c("sample","group","batch","library_id","sequencing_run"),names(object@meta.data))
  sample_table <- stats::aggregate(rep(1L,ncol(object)),object@meta.data[,sample_fields,drop=FALSE],sum)
  names(sample_table)[ncol(sample_table)] <- "cells"
  analysis_write_csv(sample_table,file.path(out,"sample_design.csv"))
  if(!is.null(design$selected) && design$selected %in% names(object@meta.data) && design$group_available)
    analysis_write_csv(as.data.frame(table(batch=object@meta.data[[design$selected]],group=object$group)),
      file.path(out,"batch_group_contingency.csv"))
  pc_used <- object@misc$workflow_analysis$pc_selection$pc_used
  dims <- seq_len(pc_used)
  object <- analysis_umap(object,"pca",dims,"umap_unintegrated",cfg$seed)
  object@misc$workflow_analysis$harmony_design <- design
  object@misc$workflow_analysis$harmony_provenance <- provenance
  analysis_save(object,file.path(out,"seurat_unintegrated.rds"))
  fields <- intersect(c("sample","group",design$selected),names(object@meta.data))
  fields <- unique(fields)
  before <- lapply(fields,function(field) analysis_embedding_panel(object,field,"umap_unintegrated",paste("Before:",field)))
  analysis_plot(patchwork::wrap_plots(before,ncol=min(3L,length(before))),"umap_before_integration",
    if(length(before)>2L) "umap_3panel" else if(length(before)==2L) "umap_2panel" else "umap_single",out)
  diagnostics <- analysis_neighbor_metrics(object,"pca",dims,"unintegrated",cfg$markers,design$selected)
  metrics <- diagnostics$metrics
  markers <- diagnostics$markers
  status <- "SKIPPED_DESIGN"
  if(isTRUE(design$confounded)) warnings <- analysis_warning(warnings,"BIOLOGICAL_TRADEOFF",
    "Batch/sample and biological group are rank-confounded; technical and group effects cannot be separated reliably",
    "Keep unintegrated results and review marker structure; any explicit Harmony comparison retains this limitation")
  if(!design$group_available) warnings <- analysis_warning(warnings,"BIOLOGICAL_TRADEOFF",
    "Biological group metadata are unavailable or incomplete; confounding cannot be checked",
    "Provide biological design metadata before interpreting corrected structure")
  if(design$eligible) {
    analysis_require("harmony")
    # Installed Harmony S3 methods differ between reduction.use and reduction.
    method <- getS3method("RunHarmony","Seurat",envir=asNamespace("harmony"))
    arguments <- list(object=object,group.by.vars=design$selected,dims.use=dims,
      reduction.save="harmony",project.dim=FALSE,verbose=FALSE)
    if("reduction.use" %in% names(formals(method))) arguments$reduction.use <- "pca"
    else if("reduction" %in% names(formals(method))) arguments$reduction <- "pca"
    else stop("Unsupported installed Harmony Seurat method: missing reduction argument")
    set.seed(cfg$seed)
    object <- do.call(harmony::RunHarmony,arguments)
    object <- analysis_umap(object,"harmony",dims,"umap_harmony",cfg$seed)
    object@misc$workflow_analysis$harmony_completed <- TRUE
    after <- lapply(fields,function(field) analysis_embedding_panel(object,field,"umap_harmony",paste("After:",field)))
    analysis_plot(patchwork::wrap_plots(after,ncol=min(3L,length(after))),"umap_after_integration",
      if(length(after)>2L) "umap_3panel" else if(length(after)==2L) "umap_2panel" else "umap_single",out)
    if(!identical(analysis_counts(object,"RNA"),original_rna)) stop("Harmony changed original RNA counts")
    post <- analysis_neighbor_metrics(object,"harmony",dims,"harmony",cfg$markers,design$selected)
    metrics <- rbind(metrics,post$metrics)
    markers <- rbind(markers,post$markers)
    present <- intersect(cfg$markers,rownames(object[[SeuratObject::DefaultAssay(object)]]))
    if(length(present)) {
      comparison <- unlist(lapply(present,function(marker) {
        lapply(c("umap_unintegrated","umap_harmony"),function(reduction) {
          plot <- Seurat::FeaturePlot(object,features=marker,reduction=reduction,
            pt.size=get_umap_point_size(ncol(object)),combine=FALSE)[[1]] + theme_sc_umap()
          plot + ggplot2::labs(title=paste(marker,reduction))
        })
      }),recursive=FALSE)
      # A real multi-page diagnostic instead of a blank marker placeholder.
      size <- sc_plot_dimensions("umap_2panel")
      grDevices::pdf(file.path(out,"marker_comparison.pdf"),width=size[1],height=size[2])
      tryCatch(for(i in seq(1L,length(comparison),by=2L))
        print(patchwork::wrap_plots(comparison[i:min(i+1L,length(comparison))],ncol=2)),
        finally=grDevices::dev.off())
    }
    analysis_save(object,file.path(out,"seurat_harmony.rds"))
    status <- "COMPLETE"
  } else warnings <- analysis_warning(warnings,"BIOLOGICAL_TRADEOFF",
    paste("Harmony was not run:",design$reason),
    "Review sample/technical metadata and biological design; use unintegrated PCA unless an explicit documented comparison is justified")
  if(length(cfg$markers) && any(markers$status=="MARKER_ABSENT"))
    warnings <- analysis_warning(warnings,"MISSING_DIAGNOSTIC",
      paste("Requested marker genes absent:",paste(unique(markers$marker[markers$status=="MARKER_ABSENT"]),collapse=",")),
      "Configure relevant marker symbols for the organism and tissue")
  analysis_write_csv(metrics,file.path(out,"integration_metrics.csv"))
  analysis_write_csv(markers,file.path(out,"marker_comparison.csv"))
  metric_plot <- ggplot2::ggplot(metrics,ggplot2::aes(x=field,y=value,fill=representation)) +
    ggplot2::geom_col(position="dodge") + ggplot2::ylim(0,1) + theme_sc() +
    ggplot2::labs(x=NULL,y="Mean same-label fraction among neighbors",title="kNN labels: descriptive, not an integration success criterion")
  marker_plot <- ggplot2::ggplot(markers[markers$status!="MARKER_ABSENT",,drop=FALSE],
    ggplot2::aes(x=marker,y=neighbor_spearman,color=representation,group=representation)) +
    ggplot2::geom_point(na.rm=TRUE) + ggplot2::labs(x=NULL,y="Marker-neighborhood Spearman") + theme_sc()
  if(any(markers$status=="MEASURED")) analysis_plot(patchwork::wrap_plots(metric_plot,marker_plot,ncol=2),
    "integration_diagnostics","umap_2panel",out)
  else analysis_plot(metric_plot,"integration_diagnostics","violin",out)
  yaml::write_yaml(list(configuration=cfg$harmony,selected_field=design$selected,decision=design$reason,
    group_confounded=design$confounded,expression_assay=object@misc$workflow_analysis$assay,
    pc_used=pc_used,seed=cfg$seed),file.path(out,"harmony_config.yaml"))
  summary <- data.frame(metric=c("cells","pc_used","batch_field","design_decision","group_confounded","harmony_status"),
    value=c(ncol(object),pc_used,design$selected %||% "NONE",design$reason,design$confounded,status))
  analysis_finish("harmony",status,out,summary,warnings)
}
