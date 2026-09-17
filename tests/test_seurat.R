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
contracts <- lapply(seq_len(nrow(m)),function(i) read_expression(m[i,,drop=FALSE],root))
contract_fields <- c("counts","reader","input_signature","source_path","input_cells","input_features","original_cell_ids")
stopifnot(all(vapply(contracts,function(x) all(contract_fields %in% names(x)) && inherits(x$counts,"dgCMatrix"),logical(1))))
stopifnot(all(vapply(seq_along(contracts),function(i) identical(contracts[[i]]$input_signature,read_signature(m[i,,drop=FALSE])),logical(1))))
bad <- m[3,,drop=FALSE]; bad$count_source <- "X"; expect_error(read_counts(bad,root))
raw <- m[3,,drop=FALSE]; raw$count_source <- "raw:X"; stopifnot(all(as.matrix(read_counts(raw,root))==as.matrix(counts[[3]])))
plain <- m[3,,drop=FALSE]; plain$local_path <- "data/GSE999999999/raw/counts.h5ad"; stopifnot(all(as.matrix(read_counts(plain,root))==as.matrix(counts[[3]])))
wrong <- m[3,,drop=FALSE]; wrong$file_type <- "10x_h5"; expect_error(read_counts(wrong,root))
transposed <- m[4,,drop=FALSE]; transposed$local_path <- "data/GSE999999999/raw/transposed.tsv"; transposed$orientation <- "cells_by_genes"; transposed$delimiter <- "tab"; transposed$feature_column <- "cell"; transposed$drop_columns <- ""
stopifnot(all(as.matrix(read_counts(transposed,root))==as.matrix(counts[[4]])))
for (path in c("counts.tsv","counts.txt","counts.txt.gz")) {
  text <- m[4,,drop=FALSE]; text$local_path <- paste0("data/GSE999999999/raw/",path); text$orientation <- "genes_by_cells"; text$delimiter <- "tab"; text$feature_column <- "gene"; text$drop_columns <- ""
  stopifnot(all(as.matrix(read_counts(text,root))==as.matrix(counts[[4]])))
}
wrong_schema <- m[2,,drop=FALSE]; wrong_schema$local_path <- "data/GSE999999999/raw/not_10x.h5"; expect_error(read_counts(wrong_schema,root))
raw_reads <- m[1,,drop=FALSE]; raw_reads$file_type <- "raw_reads"; raw_reads$local_path <- "data/GSE999999999/raw/absent.fastq.gz"; expect_error(route_input(raw_reads,root))
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
author <- Seurat::CreateSeuratObject(counts[[1]],min.cells=0,min.features=0)
author <- Seurat::NormalizeData(author,verbose=FALSE)
embedding <- matrix(seq_len(8),nrow=4,dimnames=list(colnames(author),c("PC_1","PC_2")))
author[["pca"]] <- SeuratObject::CreateDimReducObject(embeddings=embedding,key="PC_",assay="RNA")
colnames(embedding) <- c("UMAP_1","UMAP_2")
author[["umap"]] <- SeuratObject::CreateDimReducObject(embeddings=embedding,key="UMAP_",assay="RNA")
author$author_cluster <- c("a","a","b","b")
expect_error(validate_object(author))
localrds <- file.path(root,"data/GSE999999999/raw/input.rds"); saveRDS(author,localrds)
row <- m[1,,drop=FALSE]; row$local_path <- "data/GSE999999999/raw/input.rds"; row$file_type <- "rds"
stopifnot(ncol(read_counts(row,root))==4L)
sparse <- counts[[1]]
saveRDS(sparse,file.path(root,"data/GSE999999999/raw/sparse_counts.rds"))
sparse_row <- row; sparse_row$local_path <- "data/GSE999999999/raw/sparse_counts.rds"
sparse_result <- read_expression(sparse_row,root)
stopifnot(identical(as.matrix(sparse_result$counts),as.matrix(sparse)),
          identical(sparse_result$reader,"readRDS/sparseMatrix raw counts"),
          identical(sparse_result$read_timings$gunzip_seconds,0),
          is.numeric(sparse_result$read_timings$read_rds_seconds))
stopifnot(!dense_rds_memory_safe(24747485896,32000000000),
          dense_rds_memory_safe(2000000000,32000000000))
dense_frame <- as.data.frame(as.matrix(sparse))
saveRDS(dense_frame,file.path(root,"data/GSE999999999/raw/dense_counts.rds"))
dense_row <- sparse_row; dense_row$local_path <- "data/GSE999999999/raw/dense_counts.rds"
dense_result <- read_expression(dense_row,root)
stopifnot(identical(as.matrix(dense_result$counts),as.matrix(sparse)),
          identical(dense_result$reader,"readRDS/dense raw counts"),
          identical(dense_result$read_timings$object_class,"data.frame"),
          dense_result$read_timings$estimated_memory_bytes > 0,
          dense_result$read_timings$total_ram_bytes > 0)
saveRDS(as.matrix(sparse),file.path(root,"data/GSE999999999/raw/dense_counts.rds"))
stopifnot(identical(as.matrix(read_counts(dense_row,root)),as.matrix(sparse)))
original_ram_reader <- total_physical_ram_bytes
total_physical_ram_bytes <- function() 100
expect_error(read_counts(dense_row,root))
total_physical_ram_bytes <- original_ram_reader
dense_frame[1,1] <- -1
saveRDS(dense_frame,file.path(root,"data/GSE999999999/raw/dense_counts.rds"))
expect_error(read_counts(dense_row,root))
logical_sparse <- sparse != 0
saveRDS(logical_sparse,file.path(root,"data/GSE999999999/raw/sparse_counts.rds"))
stopifnot(identical(as.matrix(read_counts(sparse_row,root)),1 * as.matrix(logical_sparse)))
saveRDS(sparse,file.path(root,"data/GSE999999999/raw/sparse_counts.rds"))
gz_path <- file.path(root,"data/GSE999999999/raw/sparse_counts.rds.gz")
con <- gzfile(gz_path,"wb"); saveRDS(sparse,con); close(con)
gz_row <- sparse_row; gz_row$local_path <- "data/GSE999999999/raw/sparse_counts.rds.gz"
gz_result <- read_expression(gz_row,root)
stopifnot(identical(as.matrix(gz_result$counts),as.matrix(sparse)),
          is.numeric(gz_result$read_timings$gunzip_seconds),
          is.numeric(gz_result$read_timings$read_rds_seconds))
bad_sparse <- sparse; bad_sparse@x[1] <- -1
saveRDS(bad_sparse,file.path(root,"data/GSE999999999/raw/bad_sparse.rds"))
bad_row <- sparse_row; bad_row$local_path <- "data/GSE999999999/raw/bad_sparse.rds"
expect_error(read_counts(bad_row,root))
bad_sparse <- sparse; bad_sparse@x[1] <- 1.5
saveRDS(bad_sparse,file.path(root,"data/GSE999999999/raw/bad_sparse.rds"))
expect_error(read_counts(bad_row,root))
bad_sparse <- sparse; rownames(bad_sparse)[2] <- rownames(bad_sparse)[1]
saveRDS(bad_sparse,file.path(root,"data/GSE999999999/raw/bad_sparse.rds"))
expect_error(read_counts(bad_row,root))
bad_sparse <- sparse; colnames(bad_sparse)[2] <- colnames(bad_sparse)[1]
saveRDS(bad_sparse,file.path(root,"data/GSE999999999/raw/bad_sparse.rds"))
expect_error(read_counts(bad_row,root))
prepared <- prepare_expression_inputs(row,root)
clean <- combine_counts_and_create(prepared$counts_list,prepared$cell_sample_map,row$database)
clean <- map_metadata(clean,row,unname(prepared$cell_sample_map[SeuratObject::Cells(clean)]))
validate_object(clean,row)
stopifnot(!"author_cluster" %in% names(clean@meta.data),length(clean@reductions)==0L,identical(SeuratObject::Layers(clean[["RNA"]]),"counts"))
split_author <- author
split_author[["RNA"]] <- split(split_author[["RNA"]],f=c("A","A","B","B"))
saveRDS(split_author,file.path(root,"data/GSE999999999/raw/split_input.rds"))
split_row <- row; split_row$local_path <- "data/GSE999999999/raw/split_input.rds"; expect_error(read_counts(split_row,root))
saveRDS(list(counts=counts[[1]]),file.path(root,"data/GSE999999999/raw/not_seurat.rds")); row$local_path <- "data/GSE999999999/raw/not_seurat.rds"; expect_error(read_counts(row,root))
saveRDS(data.frame(gene=c("a","b")),file.path(root,"data/GSE999999999/raw/not_seurat.rds")); expect_error(read_counts(row,root))
cat("PASS: reader contracts; trio, H5 schema, H5AD layer/raw/gzip, CSV/TSV/TXT orientations, Seurat/sparse RDS and rds.gz, count validation, and split-layer rejection.\n")
