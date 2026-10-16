# Full independent validation of an isolated real QC run. No expression densification.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==5L,isTRUE(getOption("r45.runtime.active")))
source(file.path(args[1],"qc/run_qc.R"))
paths <- qc_paths(args[2],args[3])
out <- args[5]; dir.create(out,showWarnings=FALSE,recursive=TRUE)
raw <- readRDS(paths$raw); final <- readRDS(paths$final)
old <- readRDS(file.path(args[4],args[3],"qc/seurat_qc.rds"))
rc <- qc_validate_input(raw,args[3]); fc <- qc_validate_input(final,args[3]); oc <- qc_counts(old)
stopifnot(identical(rownames(final),rownames(raw)),all(colnames(final) %in% colnames(raw)),
          all(final$DF_adj=="Singlet"),identical(names(final@assays),"RNA"),
          identical(SeuratObject::Layers(final[["RNA"]]),"counts"),
          !length(final@reductions),!length(final@graphs),!length(final@neighbors),
          !length(SeuratObject::VariableFeatures(final)),
          !any(grepl("^seurat_clusters$|^RNA_snn_res\\.",names(final@meta.data))))
report <- read.csv(paths$report,stringsAsFactors=FALSE)
old_report <- read.csv(file.path(args[4],args[3],"qc/qc_report.csv"),stringsAsFactors=FALSE)
metric <- function(table,field,sample="ALL") tail(table$value[table$metric==field & table$sample==sample],1)
decision <- qc_decision(paths,unique(raw$sample)); cfg <- qc_config(paths,decision)
percent <- qc_add_percentages(raw,cfg$species)
retained_cell_qc <- colnames(raw)[qc_cell_keep(percent@meta.data,cfg$cell_qc)]
stopifnot(as.integer(metric(report,"after_cell_qc_cells"))==length(retained_cell_qc),
          identical(metric(report,"after_cell_qc_cells"),metric(old_report,"after_cell_qc_cells")),
          all(colnames(final) %in% retained_cell_qc))
for(field in c("sample","group","database"))
  stopifnot(identical(as.character(raw@meta.data[colnames(final),field]),as.character(final@meta.data[[field]])))
for(field in c("pMT","pRP","pHB"))
  stopifnot(isTRUE(all.equal(percent@meta.data[colnames(final),field],final@meta.data[[field]],tolerance=0)))
stopifnot(all(is.finite(final$S.Score)),all(is.finite(final$G2M.Score)),
          all(final$Phase %in% c("G1","S","G2M")),
          max(abs(final$CC.Difference-(final$S.Score-final$G2M.Score)))<1e-12)
common <- intersect(colnames(old),colnames(final))
compare_sparse <- function(a,b,cells) {
  stopifnot(identical(rownames(a),rownames(b)))
  mismatch <- 0
  for(start in seq.int(1L,length(cells),by=256L)) {
    block <- cells[start:min(start+255L,length(cells))]
    d <- a[,block,drop=FALSE]-b[,block,drop=FALSE]
    mismatch <- mismatch+sum(d@x!=0)
  }
  mismatch
}
raw_mismatch <- compare_sparse(rc,fc,colnames(final))
old_mismatch <- compare_sparse(oc,fc,common)
stopifnot(raw_mismatch==0,old_mismatch==0)
stopifnot(all(final$nCount_RNA==Matrix::colSums(fc)),all(final$nFeature_RNA==Matrix::colSums(fc>0)))
pcs <- do.call(rbind,lapply(unique(as.character(raw$sample)),function(sample) data.frame(
  sample=sample,old_pc=as.integer(metric(old_report,"pc_used",sample)),
  new_pc=as.integer(metric(report,"pc_used",sample)),max_pc=as.integer(metric(report,"max_pc",sample)),
  variance_pick=as.integer(metric(report,"variance_pick",sample)),drop_pick=as.integer(metric(report,"drop_pick",sample)),
  selected_pc_cumulative_variance=as.numeric(metric(report,"selected_pc_cumulative_variance",sample)),
  old_DF=as.integer(metric(old_report,"DF_doublets",sample)),new_DF=as.integer(metric(report,"DF_doublets",sample)),
  old_DF_adj=as.integer(metric(old_report,"DF_adj_doublets",sample)),new_DF_adj=as.integer(metric(report,"DF_adj_doublets",sample)))))
stopifnot(all(pcs$new_pc>=2 & pcs$new_pc<=pcs$max_pc),all(pcs$max_pc<=50),
          all(pcs$old_DF==pcs$new_DF),nrow(pcs)==14L)
write.csv(pcs,file.path(out,"GSE167297_pc_doublet_comparison.csv"),row.names=FALSE)
for(field in c("DF","DF_adj","Phase")) write.csv(
  as.data.frame(table(old=old@meta.data[common,field],fixed=final@meta.data[common,field],useNA="ifany")),
  file.path(out,paste0("GSE167297_old_fixed_",field,"_cross.csv")),row.names=FALSE)
exclusive_old <- setdiff(colnames(old),colnames(final)); exclusive_new <- setdiff(colnames(final),colnames(old))
write.csv(data.frame(version=c(rep("old_skill",length(exclusive_old)),rep("fixed_skill",length(exclusive_new))),
  cell=c(exclusive_old,exclusive_new)),file.path(out,"GSE167297_old_fixed_cell_differences.csv"),row.names=FALSE)
stats <- list(status="PASS",raw_cells=ncol(raw),cell_qc_cells=length(retained_cell_qc),
  old_final_cells=ncol(old),fixed_final_cells=ncol(final),genes=nrow(final),
  common_cells=length(common),old_only_cells=length(exclusive_old),fixed_only_cells=length(exclusive_new),
  raw_counts_mismatches=raw_mismatch,old_fixed_counts_mismatches=old_mismatch,
  merged_max_pc=as.integer(metric(report,"merged_max_pc")),merged_pc_used=as.integer(metric(report,"merged_pc_used")),
  merged_variance_pick=as.integer(metric(report,"merged_variance_pick")),merged_drop_pick=as.integer(metric(report,"merged_drop_pick")),
  merged_cumulative_variance=as.numeric(metric(report,"merged_selected_pc_cumulative_variance")),
  clean_object=TRUE,cell_cycle_fields_valid=TRUE,group_sample_database_unchanged=TRUE)
jsonlite::write_json(stats,file.path(out,"GSE167297_real_validation.json"),pretty=TRUE,auto_unbox=TRUE)
cat("PASS: real QC counts, cell QC stability, metadata, PCA reports, DF_adj filtering, cell cycle and clean object\n")
