# Validate actual existing-cohort checkpoints without rerunning computation.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3L) stop("Usage: validate_staged_real.R REPO TEST_ROOT GSE")
repo <- normalizePath(args[1],winslash="/",mustWork=TRUE)
project <- file.path(args[2],"results",args[3])
source(file.path(repo,"qc/run_qc.R"))
raw <- readRDS(file.path(project,"01_download/seurat_raw.rds"))
qc1 <- readRDS(file.path(project,"02_QC/01_basic_QC/seurat_QC-1.rds"))
qc2 <- readRDS(file.path(project,"02_QC/02_doublet_cellcycle/seurat_QC-2.rds"))
qc3 <- readRDS(file.path(project,"02_QC/03_ambient_RNA/seurat_QC-3.rds"))
rows <- list()
check <- function(test,condition,note="") {
  if(!isTRUE(condition)) stop("VALIDATION_FAILED: ",test)
  rows[[length(rows)+1L]] <<- data.frame(test=test,status="PASS",note=note)
}
original <- qc_counts(raw)
one <- qc_counts(qc1); two <- qc_counts(qc2); three <- qc_counts(qc3)
check("QC1 original RNA subset",identical(one,original[rownames(one),colnames(one),drop=FALSE]))
expressed <- Matrix::rowSums(original[,colnames(one),drop=FALSE]>0)
check("QC1 exact gene >3 rule",identical(rownames(one),rownames(original)[expressed>3L]))
check("QC2 original RNA subset",identical(two,one[rownames(two),colnames(two),drop=FALSE]))
check("QC2 retained DF_adj only",all(qc2$DF_adj %in% c("Singlet","NotEvaluated")))
all_calls <- read.csv(file.path(project,"02_QC/02_doublet_cellcycle/doublet_cells_all.csv"),stringsAsFactors=FALSE)
capture_verified <- "library_id" %in% names(raw@meta.data)
capture_review <- data.frame(sample=unique(as.character(raw$sample)),
  raw_library_id_present=capture_verified,
  status=if(capture_verified) "METADATA_PRESENT_REQUIRES_SOURCE_EVIDENCE" else "CAPTURE_MAPPING_WARNING",
  note=if(capture_verified) "Review source provenance for the explicit capture IDs" else
    "Historical input lacks capture IDs; this regression exercises the legacy per-sample route, not independently verified capture mapping")
write.csv(capture_review,file.path(project,"02_QC/02_doublet_cellcycle/capture_provenance_review.csv"),row.names=FALSE)
cell_field <- intersect(c("cell_id","cell"),names(all_calls))[1]
check("QC2 exact DF_adj retained cells",setequal(colnames(qc2),all_calls[[cell_field]][all_calls$DF_adj %in% c("Singlet","NotEvaluated")] ))
check("QC2 original biological metadata",all(vapply(c("sample","database","group"),function(field)
  identical(as.character(qc2@meta.data[[field]]),as.character(raw@meta.data[colnames(qc2),field])),logical(1))))
check("QC2 finite cell-cycle scores",all(is.finite(qc2$S.Score)) && all(is.finite(qc2$G2M.Score)))
for(name in c("QC1","QC2","QC3")) {
  object <- switch(name,QC1=qc1,QC2=qc2,QC3=qc3)
  counts <- qc_counts(object)
  check(paste(name,"current nCount/nFeature"),identical(as.numeric(object$nCount_RNA),as.numeric(Matrix::colSums(counts))) &&
    identical(as.numeric(object$nFeature_RNA),as.numeric(Matrix::colSums(counts>0))))
}
check("QC3 RNA and cells unchanged",identical(two,three))
corrected <- SeuratObject::LayerData(qc3,assay="decontX",layer="counts")
check("QC3 separate corrected counts readable",identical(dimnames(corrected),dimnames(three)) &&
  inherits(corrected,"sparseMatrix") && all(is.finite(corrected@x)) && all(corrected@x>=0))
check("QC3 finite contamination",all(is.finite(qc3$decontX_contamination)) && all(qc3$decontX_contamination>=0 & qc3$decontX_contamination<=1))
pca_path <- file.path(project,"03_analysis/seurat_pca.rds")
if(file.exists(pca_path)) {
  pca <- readRDS(pca_path)
  check("Formal PCA preserves RNA counts",identical(qc_counts(pca),three))
  check("Formal PCA fresh and baseline pMT",isTRUE(pca@misc$workflow_analysis$formal_pca) && identical(pca@misc$workflow_analysis$regress,"pMT"))
}
review_path <- file.path(project,"03_analysis/seurat_cluster_review.rds")
if(file.exists(review_path)) {
  review <- readRDS(review_path)
  check("Cluster review preserves RNA counts",identical(qc_counts(review),three))
  check("No unchosen formal final",!file.exists(file.path(project,"03_analysis/seurat_analysis.rds")),
    "Real-cohort resolution awaits user selection; regression does not invent approval")
}
summary <- data.frame(stage=c("raw","qc1","qc2","qc3"),
  cells=c(ncol(raw),ncol(qc1),ncol(qc2),ncol(qc3)),genes=c(nrow(raw),nrow(qc1),nrow(qc2),nrow(qc3)))
write.csv(summary,file.path(args[2],"stage_dimensions.csv"),row.names=FALSE)
write.csv(do.call(rbind,rows),file.path(args[2],"object_validation.csv"),row.names=FALSE)
cat("PASS real-cohort object invariants:",length(rows),"\n")
