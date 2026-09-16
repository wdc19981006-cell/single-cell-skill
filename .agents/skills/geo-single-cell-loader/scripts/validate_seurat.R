#!/usr/bin/env Rscript
argv <- commandArgs(trailingOnly=TRUE)
script <- sub("^--file=", "", commandArgs()[grepl("^--file=",commandArgs())][1])
source(file.path(dirname(normalizePath(script)),"seurat_common.R"))
if (length(argv)!=3L) stop("Usage: Rscript validate_seurat.R REPOSITORY_ROOT MANIFEST_RELATIVE_PATH RDS_RELATIVE_PATH")
root <- normalizePath(argv[1],winslash="/",mustWork=TRUE)
m <- read_manifest(repo_path(root,argv[2],"data"),root)
expected <- paste0("data/",m$database[1],"/seurat_raw.rds")
if (!identical(argv[3],expected)) stop("Expected final data/<GSE>/seurat_raw.rds")
object <- readRDS(repo_path(root,argv[3],paste0("data/",m$database[1])))
# Shared validation enforces one sparse integer RNA counts layer, exact manifest metadata,
# and absence of normalized/downstream assays, reductions, graphs, neighbors and extra metadata.
validate_object(object,m)
cat("Standard raw Seurat counts, metadata, cell alignment and clean structure validated\n")
