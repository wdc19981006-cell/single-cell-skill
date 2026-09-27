qc_default_doublet_api <- function() list(
  paramSweep=DoubletFinder::paramSweep,
  summarizeSweep=DoubletFinder::summarizeSweep,
  find.pK=DoubletFinder::find.pK,
  modelHomotypic=DoubletFinder::modelHomotypic,
  doubletFinder=DoubletFinder::doubletFinder
)

qc_choose_pcs <- function(stdev,cfg,max_pc) {
  if(length(stdev)<2L || !all(is.finite(stdev)) || sum(stdev)<=0) stop("PCA variance unavailable")
  pct <- stdev^2/sum(stdev^2)*100
  cumulative <- cumsum(pct)
  variance_pick <- which(cumulative >= cfg$pc_cumulative_variance & pct < cfg$pc_individual_variance)[1]
  drop_pick <- which(head(pct,-1L)-tail(pct,-1L) > cfg$pc_elbow_drop)[1]
  candidate <- c(variance_pick,if(length(drop_pick) && !is.na(drop_pick)) drop_pick+1L else NA_integer_)
  candidate <- candidate[is.finite(candidate) & candidate>=2L]
  if (!length(candidate)) return(min(20L,max_pc))
  min(min(candidate),max_pc)
}

qc_preprocess_for_doublet <- function(object,cfg,sample,cluster=TRUE) {
  result <- tryCatch({
    object <- Seurat::NormalizeData(object,verbose=FALSE)
    object <- Seurat::FindVariableFeatures(object,nfeatures=min(2000L,nrow(object)-1L),verbose=FALSE)
    max_pc <- min(as.integer(cfg$doublet$max_pcs),ncol(object)-1L,
                  length(Seurat::VariableFeatures(object))-1L,nrow(object)-1L)
    if (max_pc<3L) stop("fewer than three stable PCs")
    object <- Seurat::ScaleData(object,vars.to.regress="pMT",verbose=FALSE)
    object <- Seurat::RunPCA(object,npcs=max_pc,verbose=FALSE)
    pc_used <- qc_choose_pcs(object[["pca"]]@stdev,cfg$doublet,max_pc)
    if (pc_used<2L) stop("fewer than two usable PCs")
    if (cluster) {
      object <- Seurat::FindNeighbors(object,reduction="pca",dims=seq_len(pc_used),verbose=FALSE)
      object <- Seurat::FindClusters(object,resolution=cfg$doublet$clustering_resolution,verbose=FALSE)
    }
    list(object=object,pc_used=pc_used,max_pc=max_pc)
  },error=function(e) stop("DoubletFinder preprocessing stopped for sample ",sample,": ",conditionMessage(e),call.=FALSE))
  result
}

qc_select_pk <- function(object,pc_used,api,sample) {
  result <- tryCatch({
    sweep <- api$paramSweep(object,PCs=seq_len(pc_used),sct=FALSE)
    stats <- api$summarizeSweep(sweep,GT=FALSE)
    # find.pK plots as a side effect; keep the QC directory limited to the
    # six requested diagnostic PDFs instead of emitting Rplots.pdf in cwd.
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
  parts <- Seurat::SplitObject(object,split.by="sample")
  kept <- decision$kept_low_cell_samples %||% character()
  for (sample in names(parts)) {
    if (sample %in% kept) {
      parts[[sample]]$pANN <- rep(NA_real_,ncol(parts[[sample]]))
      parts[[sample]]$DF <- "NotEvaluated"
      parts[[sample]]$DF_adj <- "NotEvaluated"
      parts[[sample]]$doublet_status <- "SKIPPED_LOW_CELL"
      qc_report_add(report,"DOUBLET",sample,"unevaluated_cells",ncol(parts[[sample]]),status="USER_CONFIRMED_SKIP",note=decision$user_statement)
    } else {
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
  map <- Seurat::RunUMAP(processed$object,dims=seq_len(processed$pc_used),verbose=FALSE)
  qc_plot_doublets(map,paths,report)
  qc_report_add(report,"RUN",metric="doublet_umap_before_filter",value="doublet_umap.pdf; doublet_umap.png",status="PASS")
  metadata <- map@meta.data
  accepted <- metadata$DF_adj=="Singlet" | (metadata$DF_adj=="NotEvaluated" & metadata$sample %in% kept)
  if(!any(accepted)) stop("No cells remain after DF_adj filtering")
  original <- qc_counts(raw)
  retained <- rownames(metadata)[accepted]
  result <- Seurat::CreateSeuratObject(counts=original[,retained,drop=FALSE],min.cells=0,min.features=0,
                                       meta.data=metadata[retained,,drop=FALSE])
  list(object=result,pre_filter_cells=ncol(raw),unevaluated=sum(result$DF_adj=="NotEvaluated"))
}
