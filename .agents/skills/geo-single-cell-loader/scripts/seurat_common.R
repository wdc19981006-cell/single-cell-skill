if (.Platform$OS.type == "windows" && !nzchar(Sys.getenv("PROCESSOR_ARCHITECTURE"))) {
  # cli <=3.6.6 dereferences this normally-present Windows variable on unload.
  arch <- switch(R.version$arch, x86_64="AMD64", aarch64="ARM64", i386="x86", R.version$arch)
  Sys.setenv(PROCESSOR_ARCHITECTURE=arch)
}
geo_scripts_dir <- dirname(normalizePath(sys.frame(1)$ofile,winslash="/",mustWork=TRUE))
required_fields <- c("database", "sample", "tissue", "disease", "source_type", "group")
optional_fields <- c("patient", "specimen", "cohort", "treatment")
blank <- function(x) is.na(x) | tolower(trimws(as.character(x))) %in% c("", "na", "nan", "null", "none")
need <- function(packages) {
  absent <- packages[!vapply(packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(absent)) stop("Missing R packages: ", paste(absent, collapse = ", "), ". See README dependency commands.")
}
repo_path <- function(root, relative, area = "data") {
  if (blank(relative) || grepl("^[/\\\\]|^[A-Za-z]:|\\\\", relative) || ".." %in% strsplit(relative, "/", fixed=TRUE)[[1]]) stop("Unsafe relative path: ", relative)
  base <- normalizePath(file.path(root, area), winslash="/", mustWork=TRUE)
  if (!identical(tolower(base),tolower(paste0(root,"/",area)))) stop("Permitted directory redirects outside its declared location")
  path <- normalizePath(file.path(root, relative), winslash="/", mustWork=TRUE)
  if (!startsWith(tolower(path), paste0(tolower(base), "/"))) stop("Input outside permitted dataset directory")
  path
}
read_manifest <- function(path, root) {
  need(c("jsonlite"))
  m <- read.csv(path, stringsAsFactors=FALSE, check.names=FALSE, na.strings=c("", "NA"), colClasses="character", fileEncoding="UTF-8-BOM")
  fields <- c(required_fields, "local_path", "file_type", "count_source", "count_evidence", "metadata_evidence")
  if (!nrow(m) || !all(fields %in% names(m))) stop("Missing manifest fields")
  if (any(vapply(m[fields], function(x) any(blank(x)), logical(1)))) stop("Missing required manifest value")
  if (anyDuplicated(m$sample) || any(!grepl("^[A-Za-z0-9][A-Za-z0-9._-]*$", m$sample))) stop("Duplicate/unsafe sample")
  if (length(unique(m$database)) != 1L || !grepl("^GSE[0-9]+$", m$database[1])) stop("One GSE per manifest required")
  receipt <- sub("\\.csv$", ".confirmation.json", path)
  if (identical(receipt, path) || !file.exists(receipt)) stop("User confirmation missing")
  confirmation <- jsonlite::fromJSON(receipt)
  if (!identical(unname(tools::md5sum(path)), confirmation$manifest_md5) || blank(confirmation$user_statement) || !identical(confirmation$confirmed_by, "user")) stop("Manifest changed or unconfirmed")
  if (!setequal(names(confirmation$groups), m$sample) || !identical(unname(unlist(confirmation$groups[m$sample])), m$group)) stop("Confirmed groups differ from manifest")
  for (field in intersect(optional_fields, names(m))) {
    evidence <- paste0(field, "_evidence")
    if (any(!blank(m[[field]])) && (!evidence %in% names(m) || any(!blank(m[[field]]) & blank(m[[evidence]])))) stop("Missing public evidence for ", field)
  }
  repo_path(root, substring(normalizePath(path,winslash="/"), nchar(root)+2L), paste0("data/",m$database[1],"/.workflow"))
  for (i in seq_len(nrow(m))) {
    repo_path(root, m$local_path[i], paste0("data/", m$database[i],"/raw"))
  }
  cell_map_paths <- unique(m$cell_map_path[!blank(m$cell_map_path)])
  for (relative in cell_map_paths) {
    rows <- m[!blank(m$cell_map_path) & m$cell_map_path == relative,,drop=FALSE]
    expected <- unique(vapply(seq_len(nrow(rows)),function(i) value(rows[i,,drop=FALSE],"cell_map_md5"),character(1)))
    if (length(expected) != 1L || !nzchar(expected)) stop("Shared cell map has inconsistent or missing confirmation hash")
    cellmap <- repo_path(root,relative,paste0("data/",rows$database[1],"/.workflow"))
    if (!identical(unname(tools::md5sum(cellmap)),expected)) stop("Cell map changed; review and reconfirm")
  }
  m
}
value <- function(row, field, default = "") {
  if (!field %in% names(row) || blank(row[[field]][1])) default else row[[field]][1]
}
reader_signature_fields <- c("local_path", "file_type", "count_source", "delimiter", "orientation", "feature_column", "drop_columns", "assay")
reader_configuration <- function(row) {
  defaults <- c(local_path="",file_type="",count_source="",delimiter="",orientation="",feature_column="",drop_columns="",assay="RNA")
  setNames(vapply(reader_signature_fields,function(field) value(row,field,defaults[[field]]),character(1)),reader_signature_fields)
}
read_signature <- function(row) paste(paste(reader_signature_fields,reader_configuration(row),sep="="),collapse="\034")
assert_counts <- function(counts) {
  if (length(dim(counts)) != 2L || any(dim(counts) == 0L)) stop("Empty counts matrix")
  vals <- if (inherits(counts, "sparseMatrix")) counts@x else as.vector(counts)
  if (!is.numeric(vals) || any(!is.finite(vals)) || any(vals < 0) || any(abs(vals - round(vals)) > 1e-8)) stop("Counts must be finite, nonnegative integers; normalized data unsupported")
  if (is.null(rownames(counts)) || is.null(colnames(counts)) || any(blank(rownames(counts))) || any(blank(colnames(counts))) || anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts))) stop("Missing/duplicate gene or cell IDs")
  counts
}
select_rna <- function(x) {
  if (is.list(x)) {
    if (!"Gene Expression" %in% names(x)) stop("Ambiguous feature types; Gene Expression missing")
    x <- x[["Gene Expression"]]
  }
  x
}
read_counts <- function(row, root) {
  need(c("Seurat", "Matrix"))
  path <- repo_path(root, row$local_path, paste0("data/", row$database,"/raw"))
  type <- row$file_type
  chosen <- row$count_source
  if (type == "10x_mtx") {
    if (chosen != "counts") stop("10x count_source must be counts")
    hits <- function(names) unlist(lapply(names, function(n) c(file.path(path,n),file.path(path,paste0(n,".gz")))))
    for (opts in list("matrix.mtx", c("features.tsv", "genes.tsv"), "barcodes.tsv")) {
      if (sum(file.exists(hits(opts))) != 1L) stop("Missing or ambiguous 10x trio")
    }
    # Read10X requires matching compression and canonical names. Normalize only a temporary copy.
    tmp <- tempfile("trio-"); dir.create(tmp); on.exit(unlink(tmp, recursive=TRUE), add=TRUE)
    for (pair in list(c("matrix.mtx", "matrix.mtx"), c("features.tsv", "features.tsv"), c("barcodes.tsv", "barcodes.tsv"))) {
      opts <- if (pair[1] == "features.tsv") c("features.tsv", "genes.tsv") else pair[1]
      src <- hits(opts); src <- src[file.exists(src)]
      input <- if (grepl("\\.gz$", src)) gzfile(src, "rb") else file(src, "rb")
      output <- gzfile(file.path(tmp, paste0(pair[2], ".gz")), "wb")
      tryCatch(repeat { block <- readBin(input, "raw", 1048576); if (!length(block)) break; writeBin(block, output) }, finally={close(input); close(output)})
    }
    x <- select_rna(Seurat::Read10X(tmp)); reader_used <- "Seurat::Read10X"
  } else if (type == "10x_h5") {
    need("hdf5r")
    if (grepl("\\.h5ad(\\.gz)?$", path, ignore.case=TRUE)) stop("H5AD is not 10x H5")
    if (chosen != "counts") stop("10x count_source must be counts")
    h <- hdf5r::H5File$new(path, mode="r")
    valid <- tryCatch(h$exists("matrix/data") && h$exists("matrix/barcodes") && h$exists("matrix/features"), finally=h$close_all())
    if (!valid) stop("Not a supported 10x HDF5 matrix schema")
    x <- select_rna(Seurat::Read10X_h5(path)); reader_used <- "Seurat::Read10X_h5"
  } else if (type == "h5ad") {
    if (requireNamespace("zellkonverter",quietly=TRUE) && requireNamespace("SingleCellExperiment",quietly=TRUE) && requireNamespace("SummarizedExperiment",quietly=TRUE)) {
      if (grepl("\\.gz$", path)) {
        tmp <- tempfile(fileext=".h5ad"); on.exit(unlink(tmp), add=TRUE)
        input <- gzfile(path,"rb"); output <- file(tmp,"wb")
        tryCatch(repeat {block <- readBin(input,"raw",1048576); if (!length(block)) break; writeBin(block,output)}, finally={close(input);close(output)})
        path <- tmp
      }
      sce <- withCallingHandlers(zellkonverter::readH5AD(path, reader="R"), warning=function(w) stop("H5AD reader warning: ", conditionMessage(w)))
      if (startsWith(chosen, "raw:")) {
        if (!"raw" %in% SingleCellExperiment::altExpNames(sce)) stop("Requested raw matrix absent")
        sce <- SingleCellExperiment::altExp(sce,"raw"); chosen <- substring(chosen,5)
      }
      available <- SummarizedExperiment::assayNames(sce)
      if (!chosen %in% available) stop("Requested H5AD assay absent. Available: ", paste(available,collapse=", "))
      x <- SummarizedExperiment::assay(sce, chosen); reader_used <- "zellkonverter::readH5AD(reader=R)"
    } else {
      python <- Sys.getenv("GEO_SINGLE_CELL_PYTHON")
      if (!nzchar(python) || !file.exists(python)) stop("zellkonverter unavailable; set GEO_SINGLE_CELL_PYTHON to an explicit Python executable containing anndata, numpy, scipy")
      helper <- file.path(geo_scripts_dir,"h5ad_to_mtx.py")
      tmp <- tempfile("h5ad-export-")
      on.exit(unlink(tmp,recursive=TRUE),add=TRUE)
      status <- system2(python,shQuote(c(helper,path,chosen,tmp)),stdout=TRUE,stderr=TRUE)
      if (!is.null(attr(status,"status")) && attr(status,"status") != 0L) stop("anndata H5AD export failed: ",paste(status,collapse="\n"))
      x <- Matrix::readMM(file.path(tmp,"matrix.mtx"))
      rownames(x) <- readLines(file.path(tmp,"features.tsv"),encoding="UTF-8",warn=FALSE)
      colnames(x) <- readLines(file.path(tmp,"barcodes.tsv"),encoding="UTF-8",warn=FALSE)
      reader_used <- "Python anndata fallback"
    }
  } else if (type == "text") {
    need("data.table")
    if (chosen != "counts") stop("Text count_source must be counts")
    sep <- switch(value(row,"delimiter"), comma=",", tab="\t", space=" ", stop("Explicit inspected delimiter required"))
    key <- value(row,"feature_column")
    if (!nzchar(key)) stop("Explicit feature_column required")
    source <- path
    if (grepl("\\.gz$",path,ignore.case=TRUE) || identical(key,"__row_names__")) {
      source <- tempfile("geo-text-",fileext=".txt")
      input <- if (grepl("\\.gz$",path,ignore.case=TRUE)) gzfile(path,"rb") else file(path,"rb")
      output <- file(source,"wb")
      on.exit(unlink(source),add=TRUE)
      tryCatch({
        if (identical(key,"__row_names__")) writeBin(charToRaw(paste0("__feature_id__",sep)),output)
        repeat {
          block <- readBin(input,"raw",4L*1024L*1024L)
          if (!length(block)) break
          writeBin(block,output)
        }
      },finally={close(input);close(output)})
      if (identical(key,"__row_names__")) key <- "__feature_id__"
    }
    tab <- data.table::fread(source,header=TRUE,sep=sep,data.table=TRUE,check.names=FALSE,showProgress=FALSE)
    if (!key %in% names(tab) || anyDuplicated(names(tab))) stop("Missing ID column or duplicate headers")
    ids <- tab[[key]]
    dropped <- strsplit(value(row,"drop_columns"),";",fixed=TRUE)[[1]]
    dropped <- dropped[nzchar(dropped)]
    if (!all(dropped %in% names(tab)) || key %in% dropped) stop("Invalid excluded columns")
    expression_columns <- setdiff(names(tab),c(key,dropped))
    if (!length(expression_columns) || !all(vapply(tab[,expression_columns,with=FALSE],is.numeric,logical(1)))) stop("Nonexpression column found; inspect and specify drop_columns")
    data.table::set(tab,j=c(key,dropped),value=NULL)
    x <- as.matrix(tab); rownames(x) <- ids
    orientation <- value(row,"orientation")
    if (orientation == "cells_by_genes") x <- t(x) else if (orientation != "genes_by_cells") stop("Explicit inspected orientation required")
    x <- Matrix::Matrix(x, sparse=TRUE); reader_used <- "data.table::fread"
  } else if (type == "rds") {
    object <- readRDS(path)
    if (!inherits(object,"Seurat")) stop("Unsupported RDS class: ", paste(class(object),collapse=", "))
    assay <- value(row,"assay","RNA")
    if (!assay %in% names(object@assays)) stop("Requested RDS assay missing")
    if (inherits(object[[assay]],"Assay5")) {
      layers <- SeuratObject::Layers(object[[assay]])
      if (!chosen %in% layers || !startsWith(chosen,"counts")) stop("Select an explicit counts layer; split RDS layers must first be resolved")
      x <- SeuratObject::LayerData(object, assay=assay, layer=chosen)
    } else {
      if (chosen != "counts") stop("Legacy RDS requires counts")
      x <- SeuratObject::GetAssayData(object, assay=assay, slot="counts")
    }
    reader_used <- "readRDS/Seurat counts extraction"
  } else stop("Unsupported input format: ", type)
  x <- assert_counts(x); attr(x,"geo_reader") <- reader_used; x
}
map_metadata <- function(object, manifest, expected_samples) {
  if (!identical(unname(as.character(object$sample)), unname(expected_samples))) stop("Cell sample mapping changed")
  idx <- match(expected_samples, manifest$sample)
  if (anyNA(idx)) stop("Unmapped cell sample")
  for (field in c(required_fields, intersect(optional_fields,names(manifest)))) object[[field]] <- manifest[[field]][idx]
  object$orig.ident <- expected_samples
  if (!identical(unname(as.character(object$sample)), unname(expected_samples))) stop("Sample changed during mapping")
  object
}
default_cell_map_reader <- function(path) read.csv(path,stringsAsFactors=FALSE,check.names=FALSE,colClasses="character",fileEncoding="UTF-8-BOM")
prepare_expression_inputs <- function(manifest, root, reader=read_counts, cell_map_reader=default_cell_map_reader) {
  if (!nrow(manifest)) stop("Manifest has no rows")
  paths <- unique(manifest$local_path)
  counts_list <- list()
  cell_sample_map <- character()
  metrics <- vector("list",length(paths))
  reader_calls <- 0L
  cell_map_reads <- 0L
  cell_map_cache <- new.env(parent=emptyenv())
  for (input_index in seq_along(paths)) {
    local_path <- paths[[input_index]]
    rows <- manifest[manifest$local_path == local_path,,drop=FALSE]
    signatures <- unique(vapply(seq_len(nrow(rows)),function(i) read_signature(rows[i,,drop=FALSE]),character(1)))
    if (length(signatures) != 1L) stop("Shared input has inconsistent reader configuration across manifest rows.")
    samples <- rows$sample
    map_values <- vapply(seq_len(nrow(rows)),function(i) value(rows[i,,drop=FALSE],"cell_map_path"),character(1))
    if (nrow(rows) > 1L && (any(!nzchar(map_values)) || length(unique(map_values)) != 1L)) stop("Shared input requires one consistent cell_map_path across manifest rows.")
    if (any(nzchar(map_values)) && length(unique(map_values[nzchar(map_values)])) != 1L) stop("Shared input requires one consistent cell_map_path across manifest rows.")
    read_started <- proc.time()[["elapsed"]]
    counts <- reader(rows[1,,drop=FALSE],root)
    reader_calls <- reader_calls + 1L
    reader_used <- attr(counts,"geo_reader")
    if (is.null(reader_used) || !nzchar(reader_used)) reader_used <- "injected reader"
    counts <- assert_counts(counts)
    read_seconds <- proc.time()[["elapsed"]] - read_started
    mapping_started <- proc.time()[["elapsed"]]
    if (any(nzchar(map_values))) {
      map_path <- repo_path(root,unique(map_values[nzchar(map_values)]),paste0("data/",rows$database[1],"/.workflow"))
      if (exists(map_path,envir=cell_map_cache,inherits=FALSE)) {
        cellmap <- get(map_path,envir=cell_map_cache,inherits=FALSE)
      } else {
        cellmap <- cell_map_reader(map_path)
        assign(map_path,cellmap,envir=cell_map_cache)
        cell_map_reads <- cell_map_reads + 1L
      }
      if (!all(c("cell","sample") %in% names(cellmap)) || any(blank(cellmap$cell)) || any(blank(cellmap$sample)) || anyDuplicated(cellmap$cell)) stop("Invalid cell map: cell/sample must be nonblank and cell must be unique.")
      missing_cells <- setdiff(colnames(counts),cellmap$cell)
      extra_cells <- setdiff(cellmap$cell,colnames(counts))
      if (length(missing_cells) || length(extra_cells)) stop("Cell map cell set must exactly match the expression matrix columns.")
      unknown_samples <- setdiff(unique(cellmap$sample),samples)
      missing_samples <- setdiff(samples,unique(cellmap$sample))
      if (length(unknown_samples) || length(missing_samples)) stop("Cell map sample set must exactly match the manifest samples for this input.")
      sample_for_cell <- cellmap$sample[match(colnames(counts),cellmap$cell)]
    } else {
      if (nrow(rows) != 1L) stop("Shared input requires one consistent cell_map_path across manifest rows.")
      sample_for_cell <- rep(samples[[1]],ncol(counts))
    }
    new_cell_ids <- paste0(sample_for_cell,"_",colnames(counts))
    if (anyDuplicated(new_cell_ids)) stop("Duplicate prefixed barcode")
    colnames(counts) <- new_cell_ids
    counts_list[[sprintf("input_%03d",input_index)]] <- counts
    cell_sample_map <- c(cell_sample_map,setNames(sample_for_cell,new_cell_ids))
    metrics[[input_index]] <- list(local_path=local_path,file_type=rows$file_type[1],reader=reader_used,
      shared_samples=paste(samples,collapse=","),input_features=nrow(counts),input_cells=ncol(counts),
      read_seconds=read_seconds,mapping_seconds=proc.time()[["elapsed"]]-mapping_started)
  }
  if (anyDuplicated(names(cell_sample_map))) stop("Cell IDs collide across expression inputs")
  list(counts_list=counts_list,cell_sample_map=cell_sample_map,metrics=metrics,
    manifest_rows=nrow(manifest),unique_inputs=length(paths),reader_calls=reader_calls,cell_map_reads=cell_map_reads)
}
align_count_inputs <- function(counts_list) {
  need(c("Seurat", "Matrix"))
  if (!length(counts_list) || is.null(names(counts_list)) || any(blank(names(counts_list))) || anyDuplicated(names(counts_list))) stop("Counts list requires unique input names")
  reference_features <- rownames(counts_list[[1]])
  aligned <- vector("list", length(counts_list)); names(aligned) <- names(counts_list)
  for (i in seq_along(counts_list)) {
    counts <- assert_counts(counts_list[[i]])
    if (!identical(rownames(counts), reference_features)) {
      if (!setequal(rownames(counts), reference_features)) stop("Feature sets differ between samples; explicit reconciliation is required.")
      counts <- counts[match(reference_features, rownames(counts)),,drop=FALSE]
      if (!identical(rownames(counts), reference_features)) stop("Feature order reconciliation failed")
    }
    if (!inherits(counts,"sparseMatrix")) counts <- Matrix::Matrix(counts,sparse=TRUE)
    aligned[[i]] <- counts
  }
  aligned
}
combine_aligned_inputs <- function(aligned) {
  merged_counts <- if (length(aligned) == 1L) aligned[[1]] else do.call(cbind,unname(aligned))
  if (!inherits(merged_counts,"sparseMatrix")) stop("Merged counts must remain sparse")
  assert_counts(merged_counts)
}
create_seurat_from_merged <- function(merged_counts, cell_sample_map, gse) {
  need(c("Seurat","Matrix"))
  merged_counts <- assert_counts(merged_counts)
  cells <- colnames(merged_counts)
  if (anyDuplicated(cells)) stop("Cell IDs collide across samples")
  if (is.null(names(cell_sample_map)) || anyDuplicated(names(cell_sample_map)) || !setequal(names(cell_sample_map),cells)) stop("Cell-to-sample map must exactly match merged count matrix cells")
  expected_samples <- unname(cell_sample_map[cells])
  if (any(blank(expected_samples))) stop("Cell-to-sample map contains an unknown sample")
  seurat <- Seurat::CreateSeuratObject(counts=merged_counts,min.cells=3,min.features=200,project=gse)
  if (!ncol(seurat) || !nrow(seurat)) stop("Merged dataset empty after requested construction thresholds")
  retained_samples <- unname(cell_sample_map[SeuratObject::Cells(seurat)])
  if (any(blank(retained_samples))) stop("Retained cell is missing from cell-to-sample map")
  seurat$sample <- retained_samples
  seurat$orig.ident <- retained_samples
  seurat
}
combine_counts_and_create <- function(counts_list, cell_sample_map, gse) {
  aligned <- align_count_inputs(counts_list)
  merged_counts <- combine_aligned_inputs(aligned)
  create_seurat_from_merged(merged_counts,cell_sample_map,gse)
}
validate_object <- function(object, manifest=NULL) {
  need("Seurat")
  if (!inherits(object,"Seurat")) stop("Object is not Seurat")
  if (ncol(object) < 1 || nrow(object) < 1 || anyDuplicated(colnames(object))) stop("Empty object/duplicate cells")
  if (!identical(rownames(object@meta.data), SeuratObject::Cells(object)) || !identical(colnames(object), SeuratObject::Cells(object))) stop("Metadata/cell order mismatch")
  for (field in required_fields) if (!field %in% names(object@meta.data) || any(blank(object[[field]][,1]))) stop("Missing metadata: ", field)
  if (!"orig.ident" %in% names(object@meta.data) || !identical(unname(as.character(object$orig.ident)),unname(as.character(object$sample)))) stop("orig.ident/sample mismatch")
  if (any(!startsWith(colnames(object), paste0(object$sample,"_")))) stop("Sample prefix mismatch")
  if (!is.null(manifest)) {
    if (!setequal(unique(object$sample),manifest$sample)) stop("Matrix/manifest sample sets differ")
    idx <- match(object$sample,manifest$sample)
    for (field in c(required_fields,intersect(optional_fields,names(manifest)))) {
      if (!identical(unname(as.character(object[[field]][,1])), unname(as.character(manifest[[field]][idx])))) stop("Metadata mapping error: ", field)
    }
  }
  print(table(object$sample)); print(table(object$group)); print(table(object$database))
  invisible(TRUE)
}
