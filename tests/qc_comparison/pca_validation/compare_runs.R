# Compare complete diagnostic runs on identical pre-QC input and matched cell sets.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==5L,isTRUE(getOption("r45.runtime.active")))
source(file.path(args[1],"qc/run_qc.R"))
source_dir <- args[2]; latest_dir <- args[3]; control_dir <- args[4]; output <- args[5]
dir.create(output,recursive=TRUE,showWarnings=FALSE)
latest <- read.csv(file.path(latest_dir,"diagnostics/pre_filter_umap_cells.csv"),stringsAsFactors=FALSE)
control <- read.csv(file.path(control_dir,"diagnostics/pre_filter_umap_cells.csv"),stringsAsFactors=FALSE)
stopifnot(!anyDuplicated(latest$cell),!anyDuplicated(control$cell),setequal(latest$cell,control$cell))
control <- control[match(latest$cell,control$cell),,drop=FALSE]
stopifnot(identical(latest$cell,control$cell),identical(latest$sample,control$sample),
          identical(latest$group,control$group))
raw <- readRDS(file.path(source_dir,"seurat_raw.rds")); rc <- qc_counts(raw)
cfg <- qc_config(list()); raw <- qc_add_percentages(raw,cfg$species)
expected_cell_qc <- colnames(raw)[qc_cell_keep(raw@meta.data,cfg$cell_qc)]
stopifnot(setequal(latest$cell,expected_cell_qc))
summary <- list(pre_filter_same_cells=TRUE,pre_filter_cells=nrow(latest),
  raw_cells=ncol(raw),raw_genes=nrow(raw),filter_matches_raw_thresholds=TRUE)
new_metrics <- read.csv(file.path(latest_dir,"diagnostics/pca_metrics.csv"),stringsAsFactors=FALSE)
old_metrics <- read.csv(file.path(control_dir,"diagnostics/pca_metrics.csv"),stringsAsFactors=FALSE)
stopifnot(setequal(new_metrics$sample,old_metrics$sample))
spectrum_delta <- vapply(new_metrics$sample,function(sample) {
  a <- read.csv(file.path(latest_dir,"diagnostics",paste0(sample,"_pca_spectrum.csv")))
  b <- read.csv(file.path(control_dir,"diagnostics",paste0(sample,"_pca_spectrum.csv")))
  stopifnot(identical(a$PC,b$PC))
  max(abs(a$stdev-b$stdev))
},numeric(1))
summary$max_control_current_stdev_delta <- max(spectrum_delta)
summary$PCA_spectra_identical_within_1e_8 <- all(spectrum_delta<1e-8)
object_metrics <- list(); cell_sets <- list(); gene_sets <- list()
versions <- list(original=file.path(source_dir,"seurat.all_qc.rds"),
  prior_fixed=file.path(source_dir,"qc/seurat_qc.rds"),
  current=file.path(latest_dir,"qc/seurat_qc.rds"),
  legacy_pc_control=file.path(control_dir,"qc/seurat_qc.rds"))
for(version in names(versions)) {
  object <- readRDS(versions[[version]]); counts <- qc_counts(object)
  genes <- intersect(rownames(raw),rownames(object)); cells <- intersect(colnames(raw),colnames(object))
  stopifnot(length(genes)>0,length(cells)==ncol(object))
  mismatch <- 0
  for(start in seq.int(1L,length(cells),by=256L)) {
    block <- cells[start:min(start+255L,length(cells))]
    delta <- counts[genes,block,drop=FALSE]-rc[genes,block,drop=FALSE]
    mismatch <- mismatch+sum(delta@x!=0)
  }
  stopifnot(mismatch==0)
  clean <- identical(names(object@assays),"RNA") && identical(SeuratObject::Layers(object[["RNA"]]),"counts") &&
    !length(object@reductions) && !length(object@graphs) && !length(object@neighbors) &&
    !length(SeuratObject::VariableFeatures(object)) &&
    !any(grepl("^seurat_clusters$|^RNA_snn_res\\.",names(object@meta.data)))
  if(version %in% c("current","legacy_pc_control")) {
    stopifnot(clean,identical(rownames(raw),rownames(object)),all(object$DF_adj=="Singlet"),
      all(is.finite(object$S.Score)),all(is.finite(object$G2M.Score)),
      max(abs(object$CC.Difference-object$S.Score+object$G2M.Score))<1e-12,
      all(object$Phase %in% c("G1","S","G2M")))
    pre <- if(version=="current") latest else control
    retained <- pre$cell[pre$DF_adj=="Singlet"]
    stopifnot(setequal(colnames(object),retained))
    for(field in c("sample","group","database"))
      stopifnot(identical(as.character(object@meta.data[cells,field]),as.character(raw@meta.data[cells,field])))
    for(field in c("pMT","pRP","pHB"))
      stopifnot(max(abs(object@meta.data[cells,field]-raw@meta.data[cells,field]))<1e-10)
  }
  object_metrics[[version]] <- data.frame(version=version,cells=ncol(object),genes=nrow(object),
    counts_mismatches_to_raw=mismatch,clean=clean,
    DF_doublets=if("DF" %in% names(object@meta.data)) sum(object$DF=="Doublet") else NA_integer_,
    DF_adj_doublets=if("DF_adj" %in% names(object@meta.data)) sum(object$DF_adj=="Doublet") else NA_integer_)
  cell_sets[[version]] <- colnames(object);gene_sets[[version]] <- rownames(object)
  rm(object,counts);gc()
}
write.csv(do.call(rbind,object_metrics),file.path(output,"object_validation.csv"),row.names=FALSE)
set_comparisons <- list()
for(version in setdiff(names(cell_sets),"current")) {
  ac <- cell_sets[[version]];bc <- cell_sets$current
  ag <- gene_sets[[version]];bg <- gene_sets$current
  set_comparisons[[version]] <- data.frame(reference=version,current="current",
    common_cells=length(intersect(ac,bc)),reference_only_cells=length(setdiff(ac,bc)),current_only_cells=length(setdiff(bc,ac)),
    common_genes=length(intersect(ag,bg)),reference_only_genes=length(setdiff(ag,bg)),current_only_genes=length(setdiff(bg,ag)))
  for(kind in c("cells","genes")) {
    a <- if(kind=="cells") ac else ag;b <- if(kind=="cells") bc else bg
    only_a <- setdiff(a,b);only_b <- setdiff(b,a)
    write.csv(data.frame(version=c(rep(version,length(only_a)),rep("current",length(only_b))),ID=c(only_a,only_b)),
      file.path(output,paste0(version,"_vs_current_",kind,"_differences.csv")),row.names=FALSE)
  }
}
write.csv(do.call(rbind,set_comparisons),file.path(output,"cell_gene_comparisons.csv"),row.names=FALSE)
for(field in c("DF","DF_adj")) {
  write.csv(as.data.frame(table(legacy_pc_control=control[[field]],current=latest[[field]])),
    file.path(output,paste0(field,"_all_pre_filter_cross.csv")),row.names=FALSE)
  summary[[paste0(field,"_changed_cells")]] <- sum(control[[field]]!=latest[[field]])
}
summary$common_post_filter_cells <- length(intersect(cell_sets$current,cell_sets$legacy_pc_control))
summary$current_retained_cells <- length(cell_sets$current)
summary$legacy_pc_control_retained_cells <- length(cell_sets$legacy_pc_control)
summary$all_raw_count_comparisons_exact <- TRUE
summary$counts_comparison_method <- "sparse 256-cell blocks; no whole expression densification"
summary$cell_cycle_no_additional_cell_removal <- TRUE
summary$actual_historical_pre_fix_skill_available <- FALSE
summary$regression <- "pMT only, unchanged in current and legacy-PC control"

# Each version uses its own coordinate system, fixed aspect and limits. Between
# versions the cell sets are identical; their rotations are not compared as errors.
palette <- sc_sample_colors(latest$sample)
make_plot <- function(frame,title,field="sample") {
  limits <- sc_umap_limits(as.matrix(frame[,c("x","y")]))
  sc_umap_panel(frame,field,limits,length(palette),palette) +
    ggplot2::labs(title=title,subtitle=paste0(nrow(frame)," identical cell IDs")) +
    ggplot2::theme(legend.position="none")
}
for(stage in c("pre_filter","post_filter_matched")) {
  ids <- if(stage=="pre_filter") latest$cell else intersect(cell_sets$current,cell_sets$legacy_pc_control)
  a <- control[match(ids,control$cell),,drop=FALSE];b <- latest[match(ids,latest$cell),,drop=FALSE]
  stopifnot(identical(a$cell,b$cell))
  for(field in c("sample","DF_adj")) {
    p <- patchwork::wrap_plots(list(
      make_plot(a,paste0("Legacy PC rule control | ",stage),field),
      make_plot(b,paste0("Current main | ",stage),field)),ncol=2)
    stem <- if(stage=="pre_filter") "pre_filter_matched" else "post_filter_common"
    save_sc_plot(p,paste0(stem,"_",field),"umap_2panel",output,nrow(a),length(palette))
  }
}
jsonlite::write_json(summary,file.path(output,"validation.json"),pretty=TRUE,auto_unbox=TRUE)
cat("PASS: equal pre-filter inputs, matched UMAP subsets, exact sparse counts, DF_adj, clean objects, metadata and cell-cycle\n")
