#!/usr/bin/env Rscript
argv <- commandArgs(trailingOnly=TRUE)
script <- sub("^--file=", "", commandArgs()[grepl("^--file=",commandArgs())][1])
source(file.path(dirname(normalizePath(script)),"seurat_common.R"))
if (length(argv) != 2L) stop("Usage: Rscript build_seurat.R REPOSITORY_ROOT MANIFEST_RELATIVE_PATH")
root <- normalizePath(argv[1],winslash="/",mustWork=TRUE)
if (!grepl("^data/GSE[0-9]+/\\.workflow/[^/]+\\.csv$",argv[2])) stop("Manifest must be under data/<GSE>/.workflow/")
gse <- strsplit(argv[2],"/",fixed=TRUE)[[1]][2]
out <- file.path(root,"data",gse)
if (dir.exists(out) && !identical(tolower(normalizePath(out,winslash="/")),tolower(out))) stop("Dataset directory redirects outside its declared location")
workflow <- file.path(out,".workflow")
final <- file.path(out,"seurat_raw.rds")
if (file.exists(final)) stop("Output exists; archive it explicitly before rebuilding")
source(file.path(dirname(normalizePath(script)),"sample_info.R"))
tryCatch({
manifest_path <- repo_path(root,argv[2],paste0("data/",gse,"/.workflow"))
m <- read_manifest(manifest_path,root)
need(c("Seurat","Matrix"))
if (gse != m$database[1]) stop("Manifest belongs to a different GSE")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
warnings_seen <- character()
counts_list <- list()
cell_sample_map <- character()
input_details <- list()
counts_cache <- new.env(parent=emptyenv())
withCallingHandlers({
  for (i in seq_len(nrow(m))) {
    row <- m[i,,drop=FALSE]
    sharing <- m[m$local_path == row$local_path,,drop=FALSE]
    counts <- read_counts_cached(row,sharing,root,counts_cache)
    reader_used <- attr(counts,"geo_reader")
    if (nzchar(value(row,"cell_map_path"))) {
      cellmap <- read.csv(repo_path(root,row$cell_map_path,paste0("data/",gse,"/.workflow")), stringsAsFactors=FALSE,check.names=FALSE,colClasses="character")
      if (!all(c("cell","sample") %in% names(cellmap)) || anyDuplicated(cellmap$cell) || any(blank(cellmap$cell)) || any(blank(cellmap$sample))) stop("Invalid cell map")
      if (!setequal(cellmap$cell,colnames(counts)) || !setequal(cellmap$sample,sharing$sample)) stop("Cell map must match all matrix cells and all manifest samples for this input")
      idx <- match(colnames(counts),cellmap$cell)
      counts <- counts[,cellmap$sample[idx] == row$sample,drop=FALSE]
    } else if (nrow(sharing) > 1) stop("Shared matrix requires cell_map_path")
    assert_counts(counts)
    prefixed <- paste0(row$sample,"_",colnames(counts))
    if (anyDuplicated(prefixed)) stop("Duplicate prefixed barcode")
    colnames(counts) <- prefixed
    counts_list[[row$sample]] <- counts
    cell_sample_map <- c(cell_sample_map,setNames(rep(row$sample,ncol(counts)),colnames(counts)))
    input_details[[row$sample]] <- c(file_type=row$file_type,local_path=row$local_path,count_source=row$count_source,reader=reader_used,input_cells=ncol(counts))
  }
  if (anyDuplicated(names(cell_sample_map))) stop("Cell IDs collide across samples")
  seurat <- combine_counts_and_create(counts_list,cell_sample_map,gse)
  if (!setequal(unique(seurat$sample),m$sample)) stop("One or more samples are empty after requested construction thresholds")
  seurat <- map_metadata(seurat,m,unname(cell_sample_map[SeuratObject::Cells(seurat)]))
  retained <- table(seurat$sample)
  input_records <- vapply(names(input_details),function(sample) {
    detail <- input_details[[sample]]
    paste(sample,detail[["file_type"]],detail[["local_path"]],paste0("count_source=",detail[["count_source"]]),paste0("reader=",detail[["reader"]]),paste0("input_cells=",detail[["input_cells"]]),paste0("retained_cells=",unname(retained[sample])),sep="\t")
  },character(1))
  validate_object(seurat,m)
  partial <- file.path(workflow,"seurat_raw.pending.rds")
  saveRDS(seurat,partial)
  validate_object(readRDS(partial),m)
  optional <- intersect(optional_fields,names(m))
  optional <- optional[vapply(m[optional],function(x) any(!blank(x)),logical(1))]
  partial_optional <- optional[vapply(m[optional],function(x) any(blank(x)),logical(1))]
  if (length(partial_optional)) warnings_seen <- c(warnings_seen,paste("Partial optional metadata (NA retained):",paste(partial_optional,collapse=", ")))
  summary <- c(paste("GSE:",gse),paste("Created UTC:",format(Sys.time(),tz="UTC",usetz=TRUE)),paste("R:",R.version.string),paste("Seurat:",packageVersion("Seurat")),paste("SeuratObject:",packageVersion("SeuratObject")),
    "Construction strategy:","All sample count matrices were combined before CreateSeuratObject.",
    "CreateSeuratObject thresholds:","min.cells = 3","min.features = 200",
    "min.cells scope:","Entire merged GSE dataset",
    "No additional QC or downstream analysis.",
    paste("Total cells:",ncol(seurat)),paste("Total genes:",nrow(seurat)),paste("Samples:",nrow(m)),
    paste("Unique source matrices read:",length(ls(counts_cache,all.names=TRUE))),
    "Cells per sample:",capture.output(table(seurat$sample)),"Cells per group:",capture.output(table(seurat$group)),
    paste("tissue:",paste(unique(m$tissue),collapse=", ")),paste("disease:",paste(unique(m$disease),collapse=", ")),paste("source_type:",paste(unique(m$source_type),collapse=", ")),
    paste("Optional metadata:",paste(optional,collapse=", ")),"Inputs and selected count matrices:",input_records,
    "Reader policy: H5AD prefers zellkonverter native R; explicit Python anndata fallback only when unavailable. Exact reader is recorded per input above.",
    paste("Manifest:",argv[2]),paste("Manifest MD5:",unname(tools::md5sum(manifest_path))),paste("Output:",file.path("data",gse,"seurat_raw.rds")),
    "Download provenance (timestamps and SHA256):")
  download_log <- file.path(workflow,"download.json")
  summary <- c(summary,if(file.exists(download_log)) readLines(download_log,warn=FALSE) else "Local inputs: download date unavailable; not downloaded by this run.","Warnings:",if(length(warnings_seen)) unique(warnings_seen) else "None",capture.output(sessionInfo()))
  writeLines(enc2utf8(summary),file.path(workflow,"run_summary.txt"),useBytes=TRUE)
  if (!file.rename(partial,final)) stop("Could not finalize validated RDS")
  update_sample_info(out,"COMPLETE",m,seurat,input_records,warnings_seen)
},warning=function(w) {warnings_seen <<- c(warnings_seen,conditionMessage(w))})
cat("Validated output:",final,"\n")
},error=function(e) {
  info_path <- file.path(out,"sample_info.txt")
  waiting <- file.exists(info_path) && any(readLines(info_path,n=1L,warn=FALSE) == "STATUS: WAITING_FOR_GROUP_CONFIRMATION")
  if (dir.exists(out) && !waiting) tryCatch(update_sample_info(out,"BUILD_FAILED",failure=conditionMessage(e)),error=function(report_error) message("Could not update sample_info: ",conditionMessage(report_error)))
  stop(e)
})
