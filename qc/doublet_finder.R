qc_default_doublet_api <- function() list(
  paramSweep=DoubletFinder::paramSweep,
  summarizeSweep=DoubletFinder::summarizeSweep,
  find.pK=DoubletFinder::find.pK,
  modelHomotypic=DoubletFinder::modelHomotypic,
  doubletFinder=DoubletFinder::doubletFinder
)

qc_pc_selection <- function(stdev,cfg,max_pc) {
  if(length(max_pc)!=1L || !is.finite(max_pc) || max_pc<2L ||
     length(stdev)<2L || !all(is.finite(stdev)) || any(stdev<0) || sum(stdev)<=0)
    stop("PCA standard deviations unavailable or fewer than two usable PCs")
  max_pc <- min(as.integer(max_pc),length(stdev))
  stdev <- stdev[seq_len(max_pc)]
  # Reference-compatible selection uses standard-deviation share, NOT variance.
  pct <- stdev/sum(stdev)*100
  cumulative <- cumsum(pct)
  variance_pick <- which(cumulative > cfg$pc_cumulative_variance & pct < cfg$pc_individual_variance)[1]
  drop_candidates <- which(head(pct,-1L)-tail(pct,-1L) > cfg$pc_elbow_drop)
  drop_pick <- if(length(drop_candidates)) max(drop_candidates)+1L else NA_integer_
  candidate <- c(variance_pick,drop_pick)
  candidate <- candidate[is.finite(candidate)]
  method <- if(length(candidate)==2L) "reference_sd_min_cumulative_last_drop" else
    if(is.finite(variance_pick)) "reference_sd_cumulative_only" else
    if(is.finite(drop_pick)) "reference_sd_last_drop_only" else "fallback_min_20_max_pc"
  pc_used <- as.integer(max(2L,min(if(length(candidate)) min(candidate) else 20L,max_pc)))
  variance_explained <- stdev^2/sum(stdev^2)*100
  list(max_pc=max_pc,pc_used=pc_used,variance_pick=variance_pick,drop_pick=drop_pick,
    selection_method=method,standard_deviation_share=pct,variance_explained=variance_explained,
    selected_pc_cumulative_sd_share=cumulative[pc_used],
    selected_pc_cumulative_variance=sum(variance_explained[seq_len(pc_used)]))
}

qc_choose_pcs <- function(stdev,cfg,max_pc) qc_pc_selection(stdev,cfg,max_pc)$pc_used

qc_report_pcs <- function(report,sample,selection,merged=FALSE) {
  fields <- c("max_pc","pc_used","variance_pick","drop_pick","selection_method",
              "selected_pc_cumulative_sd_share","selected_pc_cumulative_variance")
  for(field in fields) qc_report_add(report,"PCA",sample,
    metric=paste0(if(merged) "merged_" else "",field),value=selection[[field]],
    note=if(field=="selected_pc_cumulative_variance") "Diagnostic only: cumulative stdev^2 / sum(stdev^2) * 100 across computed PCs" else
      "Selection uses stdev / sum(stdev) * 100; variance_pick is the legacy name for the SD-share threshold")
  for(field in c("standard_deviation_share","variance_explained")) qc_report_add(report,"PCA",sample,
    metric=paste0(if(merged) "merged_" else "",field),value=paste(signif(selection[[field]],10),collapse=";"),
    note=if(field=="variance_explained") "Diagnostic only; squared standard deviations, per computed PC (%)" else "Reference selection input; unsquared standard deviations, per computed PC (%)")
  if(selection$max_pc>=10L && selection$pc_used<=3L) qc_report_add(report,"PCA",sample,
    metric=if(merged) "merged_low_dim_warning" else "low_dim_warning",value=selection$pc_used,
    status="PCA_LOW_DIM_WARNING",note="At least 10 PCs available but only 2-3 selected; inspect diagnostics, do not force 20 PCs")
}

qc_preprocess_for_doublet <- function(object,cfg,sample,cluster=TRUE) {
  result <- tryCatch({
    object <- Seurat::NormalizeData(object,verbose=FALSE)
    object <- Seurat::FindVariableFeatures(object,nfeatures=min(2000L,nrow(object)-1L),verbose=FALSE)
    max_pc <- min(50L,as.integer(cfg$doublet$max_pcs),ncol(object)-1L,
                  length(Seurat::VariableFeatures(object))-1L,nrow(object)-1L)
    if (max_pc<3L) stop("fewer than three stable PCs")
    object <- Seurat::ScaleData(object,vars.to.regress="pMT",verbose=FALSE)
    object <- Seurat::RunPCA(object,npcs=max_pc,verbose=FALSE)
    selection <- qc_pc_selection(object[["pca"]]@stdev,cfg$doublet,max_pc)
    pc_used <- selection$pc_used
    if (pc_used<2L) stop("fewer than two usable PCs")
    if (cluster) {
      object <- Seurat::FindNeighbors(object,reduction="pca",dims=seq_len(pc_used),verbose=FALSE)
      object <- Seurat::FindClusters(object,resolution=cfg$doublet$clustering_resolution,verbose=FALSE)
    }
    list(object=object,pc_used=pc_used,max_pc=selection$max_pc,selection=selection)
  },error=function(e) stop("DoubletFinder preprocessing stopped for sample ",sample,": ",conditionMessage(e),call.=FALSE))
  result
}

qc_select_pk <- function(object,pc_used,api,sample) {
  result <- tryCatch({
    sweep <- api$paramSweep(object,PCs=seq_len(pc_used),sct=FALSE)
    stats <- api$summarizeSweep(sweep,GT=FALSE)
    # find.pK plots as a side effect; keep the QC directory limited to the
    # requested diagnostic PDFs instead of emitting Rplots.pdf in cwd.
    grDevices::pdf(file=NULL)
    on.exit(grDevices::dev.off(),add=TRUE)
    table <- api$find.pK(stats)
    values <- suppressWarnings(as.numeric(as.character(table$BCmetric)))
    pk <- suppressWarnings(as.numeric(as.character(table$pK)))
    candidates <- which(is.finite(values) & is.finite(pk) & pk>0 & pk<1)
    if(!length(candidates)) stop("find.pK returned no valid BCmetric/pK pair")
    pk[candidates[which.max(values[candidates])]]
  },error=function(e) stop("pK search stopped for sample ",sample,": ",conditionMessage(e),call.=FALSE))
  result
}

qc_new_df_column <- function(before,after,prefix,sample) {
  candidates <- setdiff(names(after@meta.data),before)
  candidates <- candidates[startsWith(candidates,prefix)]
  if(length(candidates)!=1L) stop("DoubletFinder expected one new ",prefix," column for sample ",sample)
  candidates
}

qc_run_doublet_sample <- function(object,sample,cfg,report,api) {
  original_metadata <- names(object@meta.data)
  processed <- qc_preprocess_for_doublet(object,cfg,sample)
  qc_report_pcs(report,sample,processed$selection)
  object <- processed$object
  pc_used <- processed$pc_used
  pk <- qc_select_pk(object,pc_used,api,sample)
  expected <- round(ncol(object)*cfg$doublet$expected_rate)
  if(expected<1L || expected>=ncol(object)) stop("Invalid expected doublet count for sample ",sample)
  homotypic <- api$modelHomotypic(object@meta.data[["seurat_clusters"]])
  if(length(homotypic)!=1L || !is.finite(homotypic) || homotypic<0 || homotypic>1) stop("Invalid homotypic proportion for sample ",sample)
  adjusted <- round(expected*(1-homotypic))
  if(adjusted<1L) stop("Adjusted expected doublet count is zero for sample ",sample,"; user review required")
  before <- names(object@meta.data)
  object <- api$doubletFinder(object,PCs=seq_len(pc_used),pN=cfg$doublet$pN,pK=pk,nExp=expected,reuse.pANN=NULL,sct=FALSE)
  pann <- qc_new_df_column(before,object,"pANN",sample)
  df <- qc_new_df_column(before,object,"DF.classifications",sample)
  object$pANN <- object@meta.data[[pann]]
  object$DF <- as.character(object@meta.data[[df]])
  before <- names(object@meta.data)
  object <- api$doubletFinder(object,PCs=seq_len(pc_used),pN=cfg$doublet$pN,pK=pk,nExp=adjusted,reuse.pANN="pANN",sct=FALSE)
  # DoubletFinder reuses the same classification column name when both nExp
  # values round to the same integer; the second call is still intentional.
  df_adj <- if(adjusted==expected) df else qc_new_df_column(before,object,"DF.classifications",sample)
  object$DF_adj <- as.character(object@meta.data[[df_adj]])
  if(!all(object$DF %in% c("Singlet","Doublet")) || !all(object$DF_adj %in% c("Singlet","Doublet"))) stop("Unexpected DoubletFinder classifications for sample ",sample)
  object$doublet_status <- "EVALUATED"
  metrics <- list(input_cells=ncol(object),expected_doublet_rate=cfg$doublet$expected_rate,
    nExp_poi=expected,homotypic_prop=homotypic,nExp_poi_adj=adjusted,
    DF_doublets=sum(object$DF=="Doublet"),DF_adj_doublets=sum(object$DF_adj=="Doublet"),
    DF_singlets=sum(object$DF=="Singlet"),DF_adj_singlets=sum(object$DF_adj=="Singlet"),
    selected_pK=pk,pc_used=pc_used)
  for (metric in names(metrics)) qc_report_add(report,"DOUBLET",sample,metric,metrics[[metric]],
    note=if(metric=="expected_doublet_rate") "Expected proportion, not a measured doublet rate" else "")
  # Discard temporary clustering and DoubletFinder implementation columns.
  public_metadata <- unique(c(original_metadata,"pANN","DF","DF_adj","doublet_status"))
  object@meta.data <- object@meta.data[,public_metadata,drop=FALSE]
  object
}

qc_run_doublets <- function(object,cfg,decision,paths,report,api=qc_default_doublet_api()) {
  object <- SeuratObject::JoinLayers(object)
  object <- qc_capture_libraries(object)
  for(library in unique(object$library_id)) {
    selected <- object$library_id==library
    source <- paste(unique(object$capture_mapping_source[selected]),collapse=";")
    qc_report_add(report,"CAPTURE_LIBRARY",library,"mapping_source",source,
      status=if(any(object$capture_mapping_source[selected]=="legacy_sample_mapping_unverified")) "CAPTURE_MAPPING_WARNING" else "PASS",
      note=paste("biological_samples=",paste(unique(object$sample[selected]),collapse=";"),
        "; review original capture manifest when independent library provenance is absent",sep=""))
  }
  parts <- Seurat::SplitObject(object,split.by="library_id")
  kept <- decision$kept_low_cell_samples %||% character()
  for (sample in names(parts)) {
    biological_samples <- unique(as.character(parts[[sample]]$sample))
    if(any(biological_samples %in% kept) && !all(biological_samples %in% kept))
      stop("AMBIGUOUS_CAPTURE_LIBRARY: partially skipped biological samples share library ",sample)
    if (all(biological_samples %in% kept)) {
      parts[[sample]]$pANN <- rep(NA_real_,ncol(parts[[sample]]))
      parts[[sample]]$DF <- "NotEvaluated"
      parts[[sample]]$DF_adj <- "NotEvaluated"
      parts[[sample]]$doublet_status <- "SKIPPED_LOW_CELL"
      qc_report_add(report,"DOUBLET",sample,"unevaluated_cells",ncol(parts[[sample]]),status="USER_CONFIRMED_SKIP",note=decision$user_statement)
    } else {
      if(ncol(parts[[sample]])<cfg$doublet$min_cells)
        stop("LOW_CELL_CAPTURE_LIBRARY: ",sample," has ",ncol(parts[[sample]])," cells; explicit user decision required")
      parts[[sample]] <- qc_run_doublet_sample(parts[[sample]],sample,cfg,report,api)
    }
  }
  merged <- if(length(parts)==1L) parts[[1]] else merge(parts[[1]],y=parts[-1],merge.data=FALSE)
  rm(parts); gc()
  merged <- SeuratObject::JoinLayers(merged)
  counts <- qc_counts(merged)
  metadata <- merged@meta.data[colnames(counts),,drop=FALSE]
  raw <- Seurat::CreateSeuratObject(counts=counts,min.cells=0,min.features=0,meta.data=metadata)
  rm(merged,counts); gc()
  processed <- qc_preprocess_for_doublet(raw,cfg,"ALL",cluster=FALSE)
  qc_report_pcs(report,"ALL",processed$selection,merged=TRUE)
  map <- Seurat::RunUMAP(processed$object,dims=seq_len(processed$pc_used),verbose=FALSE)
  qc_plot_doublets(map,paths,report)
  qc_report_add(report,"RUN",metric="doublet_umap_before_filter",value="doublet_umap.pdf; doublet_umap.png",status="PASS")
  metadata <- map@meta.data
  calls <- data.frame(cell_id=rownames(metadata),metadata[,intersect(c("sample","group","library_id","pANN","DF","DF_adj","doublet_status"),names(metadata)),drop=FALSE],row.names=NULL)
  utils::write.csv(calls,file.path(paths$qc,"doublet_cells_all.csv"),row.names=FALSE,na="")
  summary <- do.call(rbind,lapply(unique(as.character(metadata$sample)),function(sample) {
    selected <- metadata$sample==sample
    data.frame(sample=sample,library_ids=paste(unique(metadata$library_id[selected]),collapse=";"),
      input_cells=sum(selected),DF_doublets=sum(metadata$DF[selected]=="Doublet"),
      DF_adj_doublets=sum(metadata$DF_adj[selected]=="Doublet"),
      DF_fraction=mean(metadata$DF[selected]=="Doublet"),DF_adj_fraction=mean(metadata$DF_adj[selected]=="Doublet"),
      DF_adjustment_disagreements=sum(metadata$DF[selected]!=metadata$DF_adj[selected]),
      unevaluated_cells=sum(metadata$DF_adj[selected]=="NotEvaluated"))
  }))
  utils::write.csv(summary,file.path(paths$qc,"doublet_summary.csv"),row.names=FALSE)
  accepted <- metadata$DF_adj=="Singlet" | (metadata$DF_adj=="NotEvaluated" & metadata$sample %in% kept)
  if(!any(accepted)) stop("No cells remain after DF_adj filtering")
  original <- qc_counts(raw)
  retained <- rownames(metadata)[accepted]
  qc_plot_singlets(map,retained,paths,report)
  result <- Seurat::CreateSeuratObject(counts=original[,retained,drop=FALSE],min.cells=0,min.features=0,
                                       meta.data=metadata[retained,,drop=FALSE])
  list(object=result,pre_filter_cells=ncol(raw),unevaluated=sum(result$DF_adj=="NotEvaluated"))
}
