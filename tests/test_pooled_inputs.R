args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if(length(args)) args[1] else ".",winslash="/",mustWork=TRUE)
scripts <- file.path(root,".agents/skills/geo-single-cell-loader/scripts")
source(file.path(scripts,"seurat_common.R"))
need("Matrix")

expect_error <- function(pattern,code) {
  message <- tryCatch({force(code); ""},error=function(e) conditionMessage(e))
  if (!grepl(pattern,message,fixed=TRUE)) stop("Expected error containing '",pattern,"', got: ",message)
}

fixture_root <- tempfile("pooled-root-")
dir.create(file.path(fixture_root,"data/GSE999999999/raw"),recursive=TRUE)
workflow <- file.path(fixture_root,"data/GSE999999999/.workflow")
dir.create(workflow)
file.create(file.path(fixture_root,"data/GSE999999999/raw/pooled_counts.tsv"))
fixture_root <- normalizePath(fixture_root,winslash="/",mustWork=TRUE)
workflow <- file.path(fixture_root,"data/GSE999999999/.workflow")
map_path <- file.path(workflow,"cell_map.csv")
write.csv(data.frame(cell=paste0("cell",1:6),sample=rep(c("SampleA","SampleB","SampleC"),each=2)),map_path,row.names=FALSE,quote=FALSE)
file.copy(map_path,file.path(workflow,"other_map.csv"))
on.exit(unlink(fixture_root,recursive=TRUE),add=TRUE)

manifest <- data.frame(database="GSE999999999",sample=c("SampleA","SampleB","SampleC"),
  local_path="data/GSE999999999/raw/pooled_counts.tsv",file_type="text",count_source="counts",
  delimiter="tab",orientation="genes_by_cells",feature_column="gene",drop_columns="",assay="RNA",
  cell_map_path="data/GSE999999999/.workflow/cell_map.csv",stringsAsFactors=FALSE,check.names=FALSE)
counts <- Matrix::Matrix(matrix(1:12,nrow=2,dimnames=list(c("gene1","gene2"),paste0("cell",1:6))),sparse=TRUE)
reader_calls <- 0L
reader <- function(...) {reader_calls <<- reader_calls+1L; x <- counts; attr(x,"geo_reader") <- "stub reader"; x}
cell_map_calls <- 0L
map_reader <- function(path) {cell_map_calls <<- cell_map_calls+1L; read.csv(path,stringsAsFactors=FALSE)}
prepared <- prepare_expression_inputs(manifest,fixture_root,reader,map_reader)
expected_cells <- c("SampleA_cell1","SampleA_cell2","SampleB_cell3","SampleB_cell4","SampleC_cell5","SampleC_cell6")
stopifnot(reader_calls==1L,cell_map_calls==1L,prepared$reader_calls==1L,prepared$cell_map_reads==1L)
stopifnot(length(prepared$counts_list)==1L,identical(colnames(prepared$counts_list[[1]]),expected_cells))
stopifnot(identical(unname(prepared$cell_sample_map[expected_cells]),rep(c("SampleA","SampleB","SampleC"),each=2)))

# A user-confirmed subset retains complete source mapping but only selected columns.
subset_manifest <- manifest[1:2,,drop=FALSE]
attr(subset_manifest,"subset_confirmed") <- TRUE
attr(subset_manifest,"source_samples") <- manifest$sample
reader_calls <- 0L; cell_map_calls <- 0L
selected <- prepare_expression_inputs(subset_manifest,fixture_root,reader,map_reader)
stopifnot(reader_calls==1L,cell_map_calls==1L)
stopifnot(identical(colnames(selected$counts_list[[1]]),expected_cells[1:4]))
stopifnot(identical(unname(selected$cell_sample_map),rep(c("SampleA","SampleB"),each=2)))
unconfirmed_subset <- manifest[1:2,,drop=FALSE]
expect_error("Cell map sample set must exactly match",prepare_expression_inputs(unconfirmed_subset,fixture_root,reader,map_reader))

# A separate author metadata table is a mapping input, not a new expression reader.
separate_path <- file.path(workflow,"separate_metadata.csv")
write.csv(data.frame(barcode=paste0("cell",1:6),author_sample=rep(c("SampleA","SampleB","SampleC"),each=2),author_cluster=letters[1:6]),separate_path,row.names=FALSE,quote=FALSE)
separate <- manifest
separate$cell_map_path <- "data/GSE999999999/.workflow/separate_metadata.csv"
separate$cell_map_cell_column <- "barcode"
separate$cell_map_sample_column <- "author_sample"
reader_calls <- 0L; cell_map_calls <- 0L
prepared_separate <- prepare_expression_inputs(separate,fixture_root,reader,map_reader)
stopifnot(reader_calls==1L,cell_map_calls==1L,prepared_separate$reader_calls==1L,prepared_separate$cell_map_reads==1L)
stopifnot(identical(unname(prepared_separate$cell_sample_map[expected_cells]),rep(c("SampleA","SampleB","SampleC"),each=2)))

bad <- manifest; bad$cell_map_path <- ""
expect_error("Shared input requires one consistent cell_map_path",prepare_expression_inputs(bad,fixture_root,reader,map_reader))
bad <- manifest; bad$cell_map_path[3] <- "data/GSE999999999/.workflow/other_map.csv"
expect_error("Shared input requires one consistent cell_map_path",prepare_expression_inputs(bad,fixture_root,reader,map_reader))
bad <- manifest; bad$delimiter[3] <- "comma"
expect_error("Shared input has inconsistent reader configuration",prepare_expression_inputs(bad,fixture_root,reader,map_reader))
bad <- manifest; bad$layer <- ""; bad$layer[3] <- "counts"
expect_error("Shared input has inconsistent reader configuration",prepare_expression_inputs(bad,fixture_root,reader,map_reader))

check_map <- function(mapping,pattern) {
  custom <- function(path) mapping
  expect_error(pattern,prepare_expression_inputs(manifest,fixture_root,reader,custom))
}
valid_map <- data.frame(cell=paste0("cell",1:6),sample=rep(c("SampleA","SampleB","SampleC"),each=2),stringsAsFactors=FALSE)
check_map(valid_map[-1,],"Cell map cell set must exactly match")
check_map(rbind(valid_map,data.frame(cell="cell7",sample="SampleC")),"Cell map cell set must exactly match")
check_map(rbind(valid_map,valid_map[1,]),"Invalid cell map")
unknown <- valid_map; unknown$sample[6] <- "UnknownSample"
check_map(unknown,"Cell map sample set must exactly match")
missing <- valid_map; missing$sample[5:6] <- "SampleB"
check_map(missing,"Cell map sample set must exactly match")
cat("PASS: pooled input and separate metadata are each read once, custom mapping keys work, cells stay unsplit, and shared-input errors are rejected.\n")
