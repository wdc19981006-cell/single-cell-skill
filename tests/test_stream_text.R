args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2L) stop("Usage: test_stream_text.R REPO_ROOT PYTHON_EXECUTABLE")
root <- normalizePath(args[1],winslash="/",mustWork=TRUE)
source(file.path(root,".agents/skills/geo-single-cell-loader/scripts/seurat_common.R"))
fixture <- tempfile("stream-text-test-")
dir.create(file.path(fixture,"data/GSE999999999/raw"),recursive=TRUE)
on.exit(unlink(fixture,recursive=TRUE),add=TRUE)
fixture <- normalizePath(fixture,winslash="/",mustWork=TRUE)
path <- file.path(fixture,"data/GSE999999999/raw/counts.txt.gz")
con <- gzfile(path,"wt")
writeLines(c("c1\tc2\tc3", "g1\t0\t2\t0", "g2\t4\t0\t6", "g3\t0\t0\t7"),con)
close(con)
old <- Sys.getenv(c("GEO_SINGLE_CELL_STREAM_TEXT","GEO_SINGLE_CELL_PYTHON"),unset=NA_character_)
on.exit({
  for (name in names(old)) if (is.na(old[[name]])) Sys.unsetenv(name) else do.call(Sys.setenv,setNames(list(old[[name]]),name))
},add=TRUE)
Sys.setenv(GEO_SINGLE_CELL_STREAM_TEXT="1",GEO_SINGLE_CELL_PYTHON=args[2])
row <- data.frame(database="GSE999999999",sample="GSM1",local_path="data/GSE999999999/raw/counts.txt.gz",
  file_type="text",count_source="counts",delimiter="tab",orientation="genes_by_cells",
  feature_column="__row_names__",drop_columns="",stringsAsFactors=FALSE)
Sys.unsetenv("GEO_SINGLE_CELL_STREAM_TEXT")
stopifnot(stream_text_candidate(row,256 * 1024^2),
          !stream_text_candidate(row,256 * 1024^2 - 1))
other_shape <- row
other_shape$orientation <- "cells_by_genes"
stopifnot(!stream_text_candidate(other_shape,256 * 1024^2))
Sys.setenv(GEO_SINGLE_CELL_STREAM_TEXT="1")
result <- read_expression(row,fixture)
expected <- matrix(c(0,4,0,2,0,0,0,6,7),nrow=3,
                   dimnames=list(c("g1","g2","g3"),c("c1","c2","c3")))
stopifnot(identical(as.matrix(result$counts),expected),
          inherits(result$counts,"dgCMatrix"),
          identical(result$reader,"Python NumPy streaming + Matrix::sparseMatrix"))
con <- gzfile(path,"wt")
writeLines(c("Index\tc1\tc2\tc3", "g1\t0\t2\t0", "g2\t4\t0\t6", "g3\t0\t0\t7"),con)
close(con)
labeled <- read_expression(row,fixture)
stopifnot(identical(as.matrix(labeled$counts),expected),
          identical(colnames(labeled$counts),c("c1","c2","c3")))
con <- gzfile(path,"wt")
writeLines(c("c1\tc2", "g1\t1\t-2"),con)
close(con)
error <- tryCatch({read_expression(row,fixture); ""},error=function(e) conditionMessage(e))
stopifnot(grepl("Streaming text conversion failed",error,fixed=TRUE))
cat("PASS: streaming text reader handles labeled and unlabeled gene columns, preserves sparse counts, and rejects negative input.\n")
