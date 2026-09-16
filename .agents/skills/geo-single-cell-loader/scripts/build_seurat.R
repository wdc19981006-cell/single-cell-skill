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
seconds_since <- function(started) unname(proc.time()[["elapsed"]] - started)
fmt_seconds <- function(x) formatC(x,format="f",digits=3)
total_started <- proc.time()[["elapsed"]]
tryCatch({
manifest_started <- proc.time()[["elapsed"]]
manifest_path <- repo_path(root,argv[2],paste0("data/",gse,"/.workflow"))
m <- read_manifest(manifest_path,root)
manifest_read_seconds <- seconds_since(manifest_started)
need(c("Seurat","Matrix"))
if (gse != m$database[1]) stop("Manifest belongs to a different GSE")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
warnings_seen <- character()
withCallingHandlers({
  prepared <- prepare_expression_inputs(m,root)
  alignment_started <- proc.time()[["elapsed"]]
  aligned <- align_count_inputs(prepared$counts_list)
  feature_alignment_seconds <- seconds_since(alignment_started)
  combine_started <- proc.time()[["elapsed"]]
  merged_counts <- combine_aligned_inputs(aligned)
  matrix_combine_seconds <- seconds_since(combine_started)
  create_started <- proc.time()[["elapsed"]]
  seurat <- create_seurat_from_merged(merged_counts,prepared$cell_sample_map,gse)
  create_seurat_seconds <- seconds_since(create_started)
  metadata_started <- proc.time()[["elapsed"]]
  if (!setequal(unique(seurat$sample),m$sample)) stop("One or more samples are empty after requested construction thresholds")
  seurat <- map_metadata(seurat,m,unname(prepared$cell_sample_map[SeuratObject::Cells(seurat)]))
  retained <- table(seurat$sample)
  metric_paths <- vapply(prepared$metrics,`[[`,character(1),"local_path")
  input_records <- vapply(seq_len(nrow(m)),function(i) {
    detail <- prepared$metrics[[match(m$local_path[i],metric_paths)]]
    paste(m$sample[i],m$file_type[i],m$local_path[i],paste0("count_source=",m$count_source[i]),paste0("reader=",detail$reader),paste0("physical_input_cells=",detail$input_cells),paste0("retained_cells=",unname(retained[m$sample[i]])),sep="\t")
  },character(1))
  metadata_mapping_seconds <- seconds_since(metadata_started)
  validation_started <- proc.time()[["elapsed"]]
  validate_object(seurat,m)
  validation_seconds <- seconds_since(validation_started)
  partial <- file.path(workflow,"seurat_raw.pending.rds")
  serialization_started <- proc.time()[["elapsed"]]
  saveRDS(seurat,partial)
  serialized <- readRDS(partial)
  serialization_seconds <- seconds_since(serialization_started)
  validation_started <- proc.time()[["elapsed"]]
  validate_object(serialized,m)
  validation_seconds <- validation_seconds + seconds_since(validation_started)
  optional <- intersect(optional_fields,names(m))
  optional <- optional[vapply(m[optional],function(x) any(!blank(x)),logical(1))]
  partial_optional <- optional[vapply(m[optional],function(x) any(blank(x)),logical(1))]
  if (length(partial_optional)) warnings_seen <- c(warnings_seen,paste("Partial optional metadata (NA retained):",paste(partial_optional,collapse=", ")))
  input_performance <- unlist(lapply(seq_along(prepared$metrics),function(i) {
    detail <- prepared$metrics[[i]]
    c(paste0("Input ",i,":"),paste("local_path:",detail$local_path),paste("file_type:",detail$file_type),paste("reader:",detail$reader),
      paste("input_signature:",gsub("\034"," | ",detail$input_signature,fixed=TRUE)),paste("file_size_bytes:",detail$file_size_bytes),
      paste("shared_samples:",detail$shared_samples),paste("input_features:",detail$input_features),paste("input_cells:",detail$input_cells),
      paste("read_seconds:",fmt_seconds(detail$read_seconds)),paste("mapping_seconds:",fmt_seconds(detail$mapping_seconds)))
  }))
  total_build_seconds <- seconds_since(total_started)
  summary <- c(paste("GSE:",gse),paste("Created UTC:",format(Sys.time(),tz="UTC",usetz=TRUE)),paste("R:",R.version.string),paste("Seurat:",packageVersion("Seurat")),paste("SeuratObject:",packageVersion("SeuratObject")),
    "Construction strategy:","All unique physical expression inputs were combined before one CreateSeuratObject call.",
    "CreateSeuratObject thresholds:","min.cells = 3","min.features = 200",
    "min.cells scope:","Entire merged GSE dataset",
    "No additional QC or downstream analysis.",
    paste("Total cells:",ncol(seurat)),paste("Total genes:",nrow(seurat)),paste("Samples:",nrow(m)),
    "Cells per sample:",capture.output(table(seurat$sample)),"Cells per group:",capture.output(table(seurat$group)),
    paste("tissue:",paste(unique(m$tissue),collapse=", ")),paste("disease:",paste(unique(m$disease),collapse=", ")),paste("source_type:",paste(unique(m$source_type),collapse=", ")),
    paste("Optional metadata:",paste(optional,collapse=", ")),"Inputs and selected count matrices:",input_records,
    "Reader policy: text inputs use data.table::fread; other routes load only their own dependencies. Exact reader is recorded per physical input.",
    "Input routing:",paste("Manifest rows:",prepared$manifest_rows),paste("Unique input signatures:",prepared$unique_inputs),
    paste("Unique physical expression inputs:",prepared$unique_inputs),
    paste("Expression matrices actually read:",prepared$reader_calls),paste("Cell maps actually read:",prepared$cell_map_reads),
    "Performance:",
    paste("manifest_read_seconds:",fmt_seconds(manifest_read_seconds)),input_performance,
    paste("feature_alignment_seconds:",fmt_seconds(feature_alignment_seconds)),paste("matrix_combine_seconds:",fmt_seconds(matrix_combine_seconds)),
    paste("create_seurat_seconds:",fmt_seconds(create_seurat_seconds)),paste("metadata_mapping_seconds:",fmt_seconds(metadata_mapping_seconds)),
    paste("validation_seconds:",fmt_seconds(validation_seconds)),paste("serialization_seconds:",fmt_seconds(serialization_seconds)),
    paste("total_build_seconds:",fmt_seconds(total_build_seconds)),
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
