#!/usr/bin/env Rscript
argv <- commandArgs(trailingOnly=TRUE)
script <- sub("^--file=", "", commandArgs()[grepl("^--file=",commandArgs())][1])
source(file.path(dirname(normalizePath(script)),"seurat_common.R"))
if (length(argv) != 2L) stop("Usage: Rscript build_seurat.R REPOSITORY_ROOT MANIFEST_RELATIVE_PATH")
root <- normalizePath(argv[1],winslash="/",mustWork=TRUE)
manifest_path <- repo_path(root,argv[2],"manifests")
m <- read_manifest(manifest_path,root)
need(c("Seurat","Matrix"))
gse <- m$database[1]
out <- file.path(root,"output",gse)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
final <- file.path(out,"seurat_raw.rds")
if (file.exists(final)) stop("Output exists; archive it explicitly before rebuilding")
warnings_seen <- character()
scelist <- list()
input_records <- character()
withCallingHandlers({
  for (i in seq_len(nrow(m))) {
    row <- m[i,,drop=FALSE]
    counts <- read_counts(row,root)
    reader_used <- attr(counts,"geo_reader")
    sharing <- m[m$local_path == row$local_path,,drop=FALSE]
    if (nzchar(value(row,"cell_map_path"))) {
      cellmap <- read.csv(repo_path(root,row$cell_map_path,paste0("datasets/",gse)), stringsAsFactors=FALSE,check.names=FALSE,colClasses="character")
      if (!all(c("cell","sample") %in% names(cellmap)) || anyDuplicated(cellmap$cell) || any(blank(cellmap$cell)) || any(blank(cellmap$sample))) stop("Invalid cell map")
      if (!setequal(cellmap$cell,colnames(counts)) || !setequal(cellmap$sample,sharing$sample)) stop("Cell map must match all matrix cells and all manifest samples for this input")
      idx <- match(colnames(counts),cellmap$cell)
      counts <- counts[,cellmap$sample[idx] == row$sample,drop=FALSE]
    } else if (nrow(sharing) > 1) stop("Shared matrix requires cell_map_path")
    assert_counts(counts)
    prefixed <- paste0(row$sample,"_",colnames(counts))
    if (anyDuplicated(prefixed)) stop("Duplicate prefixed barcode")
    colnames(counts) <- prefixed
    object <- Seurat::CreateSeuratObject(counts=counts,min.cells=3,min.features=200,project=gse)
    if (!ncol(object) || !nrow(object)) stop("Sample empty after requested construction thresholds: ",row$sample)
    object$sample <- rep(row$sample,ncol(object))
    object <- map_metadata(object,m,rep(row$sample,ncol(object)))
    scelist[[row$sample]] <- object
    input_records <- c(input_records,paste(row$sample,row$file_type,row$local_path,paste0("count_source=",row$count_source),paste0("reader=",reader_used),paste0("input_cells=",ncol(counts)),paste0("retained_cells=",ncol(object)),sep="\t"))
  }
  expected <- unlist(lapply(scelist,function(s) setNames(s$sample,colnames(s))),use.names=FALSE)
  # Explicit cell-key map protects against merge reordering.
  keys <- unlist(lapply(scelist,colnames),use.names=FALSE)
  sample_by_cell <- setNames(expected,keys)
  if (anyDuplicated(keys)) stop("Cell IDs collide across samples")
  seurat <- if (length(scelist)==1L) scelist[[1]] else merge(scelist[[1]],y=scelist[-1],merge.data=FALSE)
  if (inherits(seurat[["RNA"]],"Assay5")) seurat <- SeuratObject::JoinLayers(seurat,assay="RNA")
  seurat <- map_metadata(seurat,m,unname(sample_by_cell[colnames(seurat)]))
  validate_object(seurat,m)
  partial <- file.path(out,"seurat_raw.pending.rds")
  saveRDS(seurat,partial)
  validate_object(readRDS(partial),m)
  optional <- intersect(optional_fields,names(m))
  optional <- optional[vapply(m[optional],function(x) any(!blank(x)),logical(1))]
  partial_optional <- optional[vapply(m[optional],function(x) any(blank(x)),logical(1))]
  if (length(partial_optional)) warnings_seen <- c(warnings_seen,paste("Partial optional metadata (NA retained):",paste(partial_optional,collapse=", ")))
  summary <- c(paste("GSE:",gse),paste("Created UTC:",format(Sys.time(),tz="UTC",usetz=TRUE)),paste("R:",R.version.string),paste("Seurat:",packageVersion("Seurat")),paste("SeuratObject:",packageVersion("SeuratObject")),
    "Construction thresholds: min.cells=3; min.features=200. No additional QC or downstream analysis.",
    paste("Total cells:",ncol(seurat)),paste("Total genes:",nrow(seurat)),paste("Samples:",nrow(m)),
    "Cells per sample:",capture.output(table(seurat$sample)),"Cells per group:",capture.output(table(seurat$group)),
    paste("tissue:",paste(unique(m$tissue),collapse=", ")),paste("disease:",paste(unique(m$disease),collapse=", ")),paste("source_type:",paste(unique(m$source_type),collapse=", ")),
    paste("Optional metadata:",paste(optional,collapse=", ")),"Inputs and selected count matrices:",input_records,
    "Reader policy: H5AD prefers zellkonverter native R; explicit Python anndata fallback only when unavailable. Exact reader is recorded per input above.",
    paste("Manifest:",argv[2]),paste("Manifest MD5:",unname(tools::md5sum(manifest_path))),paste("Output:",file.path("output",gse,"seurat_raw.rds")),
    "Download provenance (timestamps and SHA256):")
  download_log <- file.path(root,"logs",paste0(gse,"_download.json"))
  summary <- c(summary,if(file.exists(download_log)) readLines(download_log,warn=FALSE) else "Local inputs: download date unavailable; not downloaded by this run.","Warnings:",if(length(warnings_seen)) unique(warnings_seen) else "None",capture.output(sessionInfo()))
  writeLines(summary,file.path(out,"run_summary.txt"))
  if (!file.rename(partial,final)) stop("Could not finalize validated RDS")
},warning=function(w) {warnings_seen <<- c(warnings_seen,conditionMessage(w))})
cat("Validated output:",final,"\n")
