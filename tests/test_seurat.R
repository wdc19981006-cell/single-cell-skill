args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if(length(args)) args[1] else ".",winslash="/",mustWork=TRUE)
scripts <- file.path(root,".agents/skills/geo-single-cell-loader/scripts")
source(file.path(scripts,"seurat_common.R"))
need(c("Seurat","Matrix"))
m <- read_manifest(file.path(root,"data/GSE999999999/.workflow/sample_manifest.csv"),root)
expect_error <- function(code) {result <- tryCatch({force(code); FALSE},error=function(e) {cat("Expected rejection:",conditionMessage(e),"\n"); TRUE}); stopifnot(result)}
counts <- lapply(seq_len(nrow(m)),function(i) read_counts(m[i,,drop=FALSE],root))
stopifnot(all(vapply(counts,function(x) identical(dim(x),c(250L,4L)),logical(1))))
stopifnot(all(vapply(counts,function(x) identical(as.numeric(x[1,]),c(2,3,4,5)),logical(1))))
bad <- m[3,,drop=FALSE]; bad$count_source <- "X"; expect_error(read_counts(bad,root))
raw <- m[3,,drop=FALSE]; raw$count_source <- "raw:X"; stopifnot(all(as.matrix(read_counts(raw,root))==as.matrix(counts[[3]])))
plain <- m[3,,drop=FALSE]; plain$local_path <- "data/GSE999999999/raw/counts.h5ad"; stopifnot(all(as.matrix(read_counts(plain,root))==as.matrix(counts[[3]])))
wrong <- m[3,,drop=FALSE]; wrong$file_type <- "10x_h5"; expect_error(read_counts(wrong,root))
transposed <- m[4,,drop=FALSE]; transposed$local_path <- "data/GSE999999999/raw/transposed.tsv"; transposed$orientation <- "cells_by_genes"; transposed$delimiter <- "tab"; transposed$feature_column <- "cell"; transposed$drop_columns <- ""
stopifnot(all(as.matrix(read_counts(transposed,root))==as.matrix(counts[[4]])))
missing_header <- file.path(root,"data/GSE999999999/raw/missing_feature_header.tsv")
writeLines(c("cellA\tcellB","gene1\t1\t0","gene2\t0\t2"),missing_header,useBytes=TRUE)
headerless <- m[4,,drop=FALSE]; headerless$local_path <- "data/GSE999999999/raw/missing_feature_header.tsv"; headerless$orientation <- "genes_by_cells"; headerless$delimiter <- "tab"; headerless$feature_column <- "__row_names__"; headerless$drop_columns <- ""; headerless$count_source <- "counts"; headerless$file_type <- "text"
parsed_headerless <- read_counts(headerless,root)
stopifnot(identical(dim(parsed_headerless),c(2L,2L)),identical(rownames(parsed_headerless),c("gene1","gene2")),identical(colnames(parsed_headerless),c("cellA","cellB")))
labeled_row_names <- file.path(root,"data/GSE999999999/raw/labeled_row_names.tsv")
writeLines(c("row.names\tcellA\tcellB","gene1\t1\t0","gene2\t0\t2"),labeled_row_names,useBytes=TRUE)
headerless$local_path <- "data/GSE999999999/raw/labeled_row_names.tsv"
parsed_labeled <- read_counts(headerless,root)
stopifnot(identical(dim(parsed_labeled),c(2L,2L)),identical(rownames(parsed_labeled),c("gene1","gene2")))
cache <- new.env(parent=emptyenv()); reader_calls <- 0L
shared <- headerless[rep(1,2),,drop=FALSE]; shared$sample <- c("SharedA","SharedB"); shared$local_path <- "data/GSE999999999/raw/shared.tsv"
fake_reader <- function(row,root) {reader_calls <<- reader_calls + 1L; Matrix::Matrix(matrix(c(1,0,0,1),2,2),sparse=TRUE)}
cached_a <- read_counts_cached(shared[1,,drop=FALSE],shared,root,cache,fake_reader)
cached_b <- read_counts_cached(shared[2,,drop=FALSE],shared,root,cache,fake_reader)
stopifnot(reader_calls==1L,identical(cached_a,cached_b),length(ls(cache,all.names=TRUE))==1L)
bad_shared <- shared; bad_shared$delimiter[2] <- "comma"
expect_error(read_counts_cached(bad_shared[1,,drop=FALSE],bad_shared,root,new.env(parent=emptyenv()),fake_reader))
object <- Seurat::CreateSeuratObject(counts[[1]],min.cells=3,min.features=200)
object <- Seurat::RenameCells(object,new.names=paste0("FixtureZ_",colnames(object))); object$sample <- rep("FixtureZ",ncol(object))
object <- map_metadata(object,m[nrow(m):1,,drop=FALSE],rep("FixtureZ",ncol(object)))
stopifnot(all(object$group=="Fixture_B")); validate_object(object,m[m$sample=="FixtureZ",,drop=FALSE])
expect_error(map_metadata(object,m,rep("FixtureA",ncol(object))))
for(field in required_fields) {badmeta <- object; badmeta[[field]] <- rep(NA_character_,ncol(object)); expect_error(validate_object(badmeta))}
badmeta <- object; badmeta$group <- rep("wrong",ncol(object)); expect_error(validate_object(badmeta,m[m$sample=="FixtureZ",,drop=FALSE]))
badmeta <- object; badmeta$orig.ident <- rep("wrong",ncol(object)); expect_error(validate_object(badmeta,m[m$sample=="FixtureZ",,drop=FALSE]))
badmeta <- Seurat::RenameCells(object,new.names=sub("^FixtureZ_","",colnames(object))); expect_error(validate_object(badmeta))
badmeta <- object; rownames(badmeta@meta.data) <- rev(rownames(badmeta@meta.data)); expect_error(validate_object(badmeta))
localrds <- file.path(root,"data/GSE999999999/raw/input.rds"); saveRDS(object,localrds)
row <- m[1,,drop=FALSE]; row$local_path <- "data/GSE999999999/raw/input.rds"; row$file_type <- "rds"
stopifnot(ncol(read_counts(row,root))==4L)
saveRDS(list(counts=counts[[1]]),file.path(root,"data/GSE999999999/raw/not_seurat.rds")); row$local_path <- "data/GSE999999999/raw/not_seurat.rds"; expect_error(read_counts(row,root))
cat("PASS: real readers for trio, H5, H5AD layer/raw/gzip, text orientations/missing feature headers, shared-input caching, RDS; count rejection; stable metadata; cell prefix/alignment.\n")
