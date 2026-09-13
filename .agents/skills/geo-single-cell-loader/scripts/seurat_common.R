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
repo_path <- function(root, relative, area = "datasets") {
  if (blank(relative) || grepl("^[/\\\\]|^[A-Za-z]:|\\\\", relative) || ".." %in% strsplit(relative, "/", fixed=TRUE)[[1]]) stop("Unsafe relative path: ", relative)
  base <- normalizePath(file.path(root, area), winslash="/", mustWork=TRUE)
  path <- normalizePath(file.path(root, relative), winslash="/", mustWork=TRUE)
  if (!startsWith(tolower(path), paste0(tolower(base), "/"))) stop("Input outside datasets")
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
  for (i in seq_len(nrow(m))) repo_path(root, m$local_path[i], paste0("datasets/", m$database[i]))
  m
}
value <- function(row, field, default = "") {
  if (!field %in% names(row) || blank(row[[field]][1])) default else row[[field]][1]
}
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
  path <- repo_path(root, row$local_path, paste0("datasets/", row$database))
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
    if (chosen != "counts") stop("Text count_source must be counts")
    sep <- switch(value(row,"delimiter"), comma=",", tab="\t", space="", stop("Explicit inspected delimiter required"))
    tab <- read.table(path, header=TRUE, sep=sep, check.names=FALSE, stringsAsFactors=FALSE, comment.char="", quote="\"", row.names=NULL)
    key <- value(row,"feature_column")
    if (!key %in% names(tab) || anyDuplicated(names(tab))) stop("Missing ID column or duplicate headers")
    ids <- tab[[key]]
    dropped <- strsplit(value(row,"drop_columns"),";",fixed=TRUE)[[1]]
    dropped <- dropped[nzchar(dropped)]
    if (!all(dropped %in% names(tab)) || key %in% dropped) stop("Invalid excluded columns")
    numeric <- tab[setdiff(names(tab), c(key,dropped))]
    if (!ncol(numeric) || !all(vapply(numeric,is.numeric,logical(1)))) stop("Nonexpression column found; inspect and specify drop_columns")
    x <- as.matrix(numeric); rownames(x) <- ids
    orientation <- value(row,"orientation")
    if (orientation == "cells_by_genes") x <- t(x) else if (orientation != "genes_by_cells") stop("Explicit inspected orientation required")
    x <- Matrix::Matrix(x, sparse=TRUE); reader_used <- "base::read.table"
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
  if (!identical(unname(as.character(object$sample)), unname(expected_samples))) stop("Sample changed during mapping")
  object
}
validate_object <- function(object, manifest=NULL) {
  need("Seurat")
  if (!inherits(object,"Seurat")) stop("Object is not Seurat")
  if (ncol(object) < 1 || nrow(object) < 1 || anyDuplicated(colnames(object))) stop("Empty object/duplicate cells")
  if (!identical(rownames(object@meta.data), SeuratObject::Cells(object)) || !identical(colnames(object), SeuratObject::Cells(object))) stop("Metadata/cell order mismatch")
  for (field in required_fields) if (!field %in% names(object@meta.data) || any(blank(object[[field]][,1]))) stop("Missing metadata: ", field)
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
