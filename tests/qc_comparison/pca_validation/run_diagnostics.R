# Test-only observations around the current fixed QC pipeline. Production code is unchanged.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==5L,isTRUE(getOption("r45.runtime.active")))
repo <- normalizePath(args[1],winslash="/",mustWork=TRUE)
root <- normalizePath(args[2],winslash="/",mustWork=TRUE)
gse <- args[3]; diagnostic_dir <- args[4]; mode <- args[5]
stopifnot(mode %in% c("latest","legacy_pc_control"))
dir.create(diagnostic_dir,recursive=TRUE,showWarnings=FALSE)
source(file.path(repo,"qc/run_qc.R"))
set.seed(155698L)

# The legacy control reproduces only the old selector on the identical modern
# pipeline. It is NOT labelled as a historical full old-Skill result.
current_selection <- qc_pc_selection
if(mode=="legacy_pc_control") qc_pc_selection <- function(stdev,cfg,max_pc) {
  selected <- current_selection(stdev,cfg,max_pc)
  pct <- stdev/sum(stdev)*100
  vp <- which(cumsum(pct)>=cfg$pc_cumulative_variance & pct<cfg$pc_individual_variance)[1]
  dp <- which(head(pct,-1L)-tail(pct,-1L)>cfg$pc_elbow_drop)[1]
  dp <- if(is.finite(dp)) dp+1L else NA_integer_
  candidates <- c(vp,dp); candidates <- candidates[is.finite(candidates) & candidates>=2L]
  selected$pc_used <- as.integer(if(length(candidates)) min(min(candidates),max_pc) else min(20L,max_pc))
  selected$variance_pick <- vp; selected$drop_pick <- dp
  selected$selection_method <- "CONTROL_reproduced_pre_cec9f45_first_drop"
  selected$selected_pc_cumulative_sd_share <- sum(selected$standard_deviation_share[seq_len(selected$pc_used)])
  selected$selected_pc_cumulative_variance <- sum(selected$variance_explained[seq_len(selected$pc_used)])
  selected
}

observations <- new.env(parent=emptyenv()); observations$rows <- list()
fixed_preprocess <- qc_preprocess_for_doublet
qc_preprocess_for_doublet <- function(object,cfg,sample,cluster=TRUE) {
  result <- fixed_preprocess(object,cfg,sample,cluster)
  stopifnot(result$pc_used>=2L,result$pc_used<=result$max_pc,result$max_pc<=50L,
    result$max_pc<=ncol(object)-1L,result$max_pc<=nrow(object)-1L,
    result$max_pc<=length(SeuratObject::VariableFeatures(result$object))-1L,
    !"umap" %in% names(result$object@reductions))
  s <- result$selection; stdev <- result$object[["pca"]]@stdev
  frame <- data.frame(sample=sample,PC=seq_along(stdev),stdev=stdev,
    sd_share=s$standard_deviation_share,variance_explained=s$variance_explained,
    cumulative_sd_share=cumsum(s$standard_deviation_share),
    cumulative_variance=cumsum(s$variance_explained),selected=seq_along(stdev)<=s$pc_used)
  write.csv(frame,file.path(diagnostic_dir,paste0(sample,"_pca_spectrum.csv")),row.names=FALSE)
  # ElbowPlot's y-axis is standard deviation, not variance explained.
  elbow <- Seurat::ElbowPlot(result$object,ndims=result$max_pc,reduction="pca") +
    ggplot2::geom_vline(xintercept=s$pc_used,color="#ea5455",linetype="dashed") +
    ggplot2::labs(title=paste0(sample," | selected PC ",s$pc_used),
      subtitle=paste0(mode,"\nSD cumulative ",round(s$selected_pc_cumulative_sd_share,2),
        "%; squared-SD cumulative ",round(s$selected_pc_cumulative_variance,2),"%"),
      y="PCA standard deviation")
  save_sc_plot(elbow,paste0("elbow_",sample),"umap_single",diagnostic_dir)
  observations$rows[[sample]] <- data.frame(sample=sample,cells=ncol(object),genes=nrow(object),
    HVGs=length(SeuratObject::VariableFeatures(result$object)),max_pc=s$max_pc,
    pc_used=s$pc_used,variance_pick=s$variance_pick,drop_pick=s$drop_pick,
    cumulative_sd_share=s$selected_pc_cumulative_sd_share,
    cumulative_variance_explained=s$selected_pc_cumulative_variance,
    selection_method=s$selection_method,mode=mode)
  write.csv(do.call(rbind,observations$rows),file.path(diagnostic_dir,"pca_metrics.csv"),row.names=FALSE)
  cat("PCA_DIAGNOSTIC",sample,"cells",ncol(object),"max_pc",s$max_pc,"selected",s$pc_used,"\n")
  flush.console()
  result
}
fixed_doublet_plot <- qc_plot_doublets
qc_plot_doublets <- function(object,paths,report) {
  embedding <- SeuratObject::Embeddings(object[["umap"]])
  meta <- object@meta.data[rownames(embedding),,drop=FALSE]
  fields <- intersect(c("sample","group","DF","DF_adj","doublet_status","pANN","nCount_RNA","nFeature_RNA"),names(meta))
  frame <- data.frame(cell=rownames(embedding),x=embedding[,1],y=embedding[,2],meta[,fields,drop=FALSE])
  write.csv(frame,file.path(diagnostic_dir,"pre_filter_umap_cells.csv"),row.names=FALSE)
  fixed_doublet_plot(object,paths,report)
}
fixed_singlet_plot <- qc_plot_singlets
qc_plot_singlets <- function(object,retained,paths,report) {
  before <- read.csv(file.path(diagnostic_dir,"pre_filter_umap_cells.csv"),stringsAsFactors=FALSE)
  subset <- before[match(retained,before$cell),,drop=FALSE]
  embedding <- SeuratObject::Embeddings(object[["umap"]])
  stopifnot(!anyNA(subset$cell),max(abs(subset$x-embedding[retained,1]))<1e-10,
    max(abs(subset$y-embedding[retained,2]))<1e-10,
    all(subset$DF_adj %in% c("Singlet","NotEvaluated")))
  write.csv(subset,file.path(diagnostic_dir,"post_filter_umap_cells.csv"),row.names=FALSE)
  fixed_singlet_plot(object,retained,paths,report)
}
api <- qc_default_doublet_api()
installed_sweep <- api$paramSweep
api$paramSweep <- function(seu,PCs,sct) {
  stopifnot(length(unique(seu$sample))==1L,!"umap" %in% names(seu@reductions))
  installed_sweep(seu,PCs=PCs,sct=sct)
}
status <- qc_run(root,gse,df_api=api)
report <- read.csv(qc_paths(root,gse)$report,stringsAsFactors=FALSE)
metrics <- do.call(rbind,observations$rows)
read_metric <- function(sample,key) {
  section <- if(sample=="ALL") "RUN" else "DOUBLET"
  values <- report$value[report$sample==sample & report$section==section & report$metric==key]
  if(length(values)) as.numeric(tail(values,1)) else NA_real_
}
for(field in c("selected_pK","DF_doublets","DF_singlets","DF_adj_doublets","DF_adj_singlets"))
  metrics[[field]] <- vapply(metrics$sample,function(s) read_metric(s,field),numeric(1))
before <- read.csv(file.path(diagnostic_dir,"pre_filter_umap_cells.csv"),stringsAsFactors=FALSE)
for(field in c("DF","DF_adj")) for(value in c("Doublet","Singlet"))
  metrics[metrics$sample=="ALL",paste(field,tolower(paste0(value,"s")),sep="_")] <- sum(before[[field]]==value)
write.csv(metrics,file.path(diagnostic_dir,"pca_metrics.csv"),row.names=FALSE)
jsonlite::write_json(list(status=status,mode=mode,seed=155698L,
  per_sample_sweep_independent=TRUE,no_per_sample_umap=TRUE,
  exact_post_filter_embedding_subset=TRUE,merged_pK="NOT_APPLICABLE: no merged DoubletFinder"),
  file.path(diagnostic_dir,"diagnostic_run.json"),pretty=TRUE,auto_unbox=TRUE)
cat("DIAGNOSTIC_QC",mode,status,"\n")
