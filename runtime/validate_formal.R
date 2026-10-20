# Separate R process, launched only through runtime/r45/run_r45.py.
script <- sub("^--file=","",commandArgs()[grepl("^--file=",commandArgs())][1])
source(file.path(dirname(normalizePath(script)),"seurat_storage.R"))
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2L) stop("Usage: validate_formal.R RDS qc|analysis")
object <- readRDS(args[1])
sc_verify_formal(object,operations=args[2]=="analysis")
cat("PASS: independent formal readback; exact expression/metadata hashes, valid Seurat, downstream APIs\n")
