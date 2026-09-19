args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if(length(args)) args[1] else ".",winslash="/",mustWork=TRUE)
scripts <- file.path(root,".agents/skills/geo-single-cell-loader/scripts")
source(file.path(scripts,"seurat_common.R"))
need(c("Seurat","Matrix"))

features <- c(paste0("base_",seq_len(200)),"gene_shared_low","gene_two_cells")
make_counts <- function(sample) {
  cells <- paste0(sample,"_cell",seq_len(3))
  x <- matrix(0,nrow=length(features),ncol=length(cells),dimnames=list(features,cells))
  x[seq_len(200),] <- 1
  x["gene_shared_low",1] <- 1
  if (sample %in% c("Sample_A","Sample_B")) x["gene_two_cells",1] <- 1
  Matrix::Matrix(x,sparse=TRUE)
}
counts_list <- setNames(lapply(c("Sample_A","Sample_B","Sample_C"),make_counts),c("Sample_A","Sample_B","Sample_C"))
# Exercise safe reordering when sets match but order differs.
counts_list[["Sample_B"]] <- counts_list[["Sample_B"]][rev(seq_len(nrow(counts_list[["Sample_B"]]))),,drop=FALSE]
cell_sample_map <- unlist(lapply(names(counts_list),function(sample) setNames(rep(sample,ncol(counts_list[[sample]])),colnames(counts_list[[sample]]))),use.names=TRUE)
reference_merged <- combine_aligned_inputs(align_count_inputs(counts_list))
object <- combine_counts_and_create(counts_list,cell_sample_map,"GSE999999999")
stopifnot(ncol(object)==9L,nrow(object)==201L)
stopifnot("gene-shared-low" %in% rownames(object),!"gene-two-cells" %in% rownames(object))
stopifnot(identical(unname(as.character(object$sample)),unname(as.character(object$orig.ident))))
stopifnot(all(table(object$sample)==3L),anyDuplicated(colnames(object))==0L)
stopifnot(length(SeuratObject::Layers(object[["RNA"]]))==1L)
reference_kept <- reference_merged[Matrix::rowSums(reference_merged != 0) >= 3,
                                   Matrix::colSums(reference_merged != 0) >= 200,drop=FALSE]
rownames(reference_kept) <- gsub("_","-",rownames(reference_kept),fixed=TRUE)
observed <- SeuratObject::LayerData(object,assay="RNA",layer="counts")
stopifnot(identical(dim(observed),dim(reference_kept)),
          identical(rownames(observed),rownames(reference_kept)),
          identical(colnames(observed),colnames(reference_kept)),
          identical(as.matrix(observed),as.matrix(reference_kept)))

# The previous per-sample filtering loses the feature because it is expressed in only one cell per sample.
old_objects <- lapply(counts_list,function(x) Seurat::CreateSeuratObject(x,min.cells=3,min.features=200))
stopifnot(all(!vapply(old_objects,function(x) "gene-shared-low" %in% rownames(x),logical(1))))

mismatch <- counts_list
rownames(mismatch[["Sample_C"]])[1] <- "different_feature"
message <- tryCatch({combine_counts_and_create(mismatch,cell_sample_map,"GSE999999999"); ""},error=function(e) conditionMessage(e))
stopifnot(identical(message,"Feature sets differ between samples; explicit reconciliation is required."))
cat("PASS: global min.cells keeps 3-cell cross-sample gene, removes 2-cell gene, reorders identical feature sets, and rejects differing sets.\n")
