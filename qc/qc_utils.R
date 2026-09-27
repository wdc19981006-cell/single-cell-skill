qc_script_dir <- dirname(normalizePath(sys.frame(1)$ofile, winslash="/", mustWork=TRUE))
qc_required_packages <- c("Seurat", "SeuratObject", "Matrix", "yaml", "jsonlite", "ggplot2", "patchwork", "DoubletFinder")
qc_blank <- function(x) is.na(x) | !nzchar(trimws(as.character(x)))

qc_paths <- function(root, gse) {
  if (length(gse) != 1L || !grepl("^GSE[0-9]+$", gse)) stop("Expected one GSE accession")
  root <- normalizePath(root, winslash="/", mustWork=TRUE)
  dataset <- file.path(root, "data", gse)
  if (!dir.exists(dataset) || !identical(normalizePath(dataset,winslash="/"), dataset)) stop("Dataset directory missing or redirected")
  qc <- file.path(dataset, "qc")
  dir.create(qc, showWarnings=FALSE)
  list(root=root, gse=gse, raw=file.path(dataset,"seurat_raw.rds"), qc=qc,
       report=file.path(qc,"qc_report.csv"), decision=file.path(qc,"qc_decision.json"),
       final=file.path(qc,"seurat_qc.rds"))
}

qc_report_new <- function() {
  report <- new.env(parent=emptyenv())
  report$rows <- list()
  report
}
qc_report_add <- function(report, section, sample="ALL", metric, value, status="PASS", note="") {
  report$rows[[length(report$rows)+1L]] <- data.frame(
    section=as.character(section), sample=as.character(sample), metric=as.character(metric),
    value=as.character(value), status=as.character(status), note=as.character(note),
    stringsAsFactors=FALSE)
  invisible(report)
}
qc_report_write <- function(report, path) {
  rows <- if (length(report$rows)) do.call(rbind,report$rows) else
    data.frame(section=character(),sample=character(),metric=character(),value=character(),status=character(),note=character())
  temporary <- paste0(path,".pending")
  utils::write.csv(rows,temporary,row.names=FALSE,na="",fileEncoding="UTF-8")
  if (file.exists(path) && !file.remove(path)) stop("Cannot replace qc_report.csv")
  if (!file.rename(temporary,path)) stop("Cannot finalize qc_report.csv")
  invisible(rows)
}

qc_config <- function(paths, decision=list()) {
  cfg <- yaml::read_yaml(file.path(qc_script_dir,"default_config.yaml"))
  overrides <- decision$parameter_overrides
  if (!is.null(overrides)) {
    allowed <- list(species=TRUE,cell_qc=names(cfg$cell_qc),doublet=c("min_cells","max_pcs","expected_rate"))
    if (any(!names(overrides) %in% names(allowed))) stop("Unsupported QC parameter override")
    for (section in names(overrides)) {
      if (section == "species") { cfg$species <- overrides$species; next }
      if (!is.list(overrides[[section]]) || any(!names(overrides[[section]]) %in% allowed[[section]])) stop("Unsupported QC parameter override")
      cfg[[section]] <- utils::modifyList(cfg[[section]],overrides[[section]])
    }
  }
  if (!cfg$species %in% c("human","mouse")) stop("QC species must be human or mouse")
  if (!identical(cfg$sample_column,"sample") || !identical(cfg$doublet$final_filter,"DF_adj") ||
      !identical(cfg$cell_cycle$regression,"none") || !isTRUE(cfg$cell_cycle$score_only) ||
      !isTRUE(cfg$doublet$enabled) || !isTRUE(cfg$cell_cycle$enabled)) stop("Required QC pipeline settings changed")
  for (entry in c("min_cells","max_pcs")) if (!is.numeric(cfg$doublet[[entry]]) || cfg$doublet[[entry]] < 2) stop("Invalid doublet configuration")
  if (!is.numeric(cfg$doublet$expected_rate) || cfg$doublet$expected_rate <= 0 || cfg$doublet$expected_rate >= 1) stop("Invalid expected doublet rate")
  cfg
}

qc_decision <- function(paths, samples) {
  if (!file.exists(paths$decision)) return(list())
  d <- jsonlite::fromJSON(paths$decision,simplifyVector=TRUE)
  if (!identical(d$confirmed_by,"user") || is.null(d$user_statement) || qc_blank(d$user_statement)) stop("QC decision requires the literal user statement and confirmed_by=user")
  if (is.null(d$decision_type) || !length(d$decision_type) || is.null(d$affected_samples)) stop("QC decision requires decision_type and affected_samples")
  removed <- as.character(d$removed_samples %||% character())
  kept <- as.character(d$kept_low_cell_samples %||% character())
  affected <- as.character(d$affected_samples)
  types <- as.character(d$decision_type)
  if (!length(affected) || anyDuplicated(affected) || anyDuplicated(c(removed,kept)) ||
      any(!c(removed,kept,affected) %in% samples) ||
      any(!types %in% c("remove_low_cell_samples","keep_without_doubletfinder","override_parameters","accept_default_conflicts")))
    stop("QC decision sample list or type invalid")
  if (length(c(removed,kept)) && !setequal(c(removed,kept), affected))
    stop("QC affected_samples must exactly match removed/kept samples")
  if (length(removed) && !"remove_low_cell_samples" %in% d$decision_type) stop("Removal decision_type missing")
  if (length(kept) && !"keep_without_doubletfinder" %in% d$decision_type) stop("Keep decision_type missing")
  if (!is.null(d$parameter_overrides) && !"override_parameters" %in% d$decision_type) stop("Override decision_type missing")
  d$removed_samples <- removed
  d$kept_low_cell_samples <- kept
  d
}
`%||%` <- function(x,y) if (is.null(x)) y else x

qc_counts <- function(object) {
  if (!inherits(object,"Seurat") || !"RNA" %in% names(object@assays)) stop("FAILED_INPUT: Seurat object with RNA assay required")
  layers <- SeuratObject::Layers(object[["RNA"]])
  if (!"counts" %in% layers || any(startsWith(layers,"counts.") & layers != "counts")) stop("FAILED_INPUT: one unambiguous RNA counts layer required")
  counts <- SeuratObject::LayerData(object,assay="RNA",layer="counts")
  if (!inherits(counts,"sparseMatrix") || length(dim(counts)) != 2L || any(dim(counts)==0L) ||
      is.null(rownames(counts)) || is.null(colnames(counts)) || any(qc_blank(rownames(counts))) ||
      any(qc_blank(colnames(counts))) || anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts)) ||
      any(!is.finite(counts@x)) || any(counts@x < 0) || any(counts@x != round(counts@x)))
    stop("FAILED_INPUT: RNA counts must be sparse, finite nonnegative integers with unique gene/cell IDs")
  counts
}
qc_validate_input <- function(object, gse) {
  counts <- qc_counts(object)
  metadata <- object@meta.data
  if (!identical(rownames(metadata),colnames(counts))) stop("FAILED_INPUT: metadata and count cells differ")
  for (field in c("sample","database","group"))
    if (!field %in% names(metadata) || any(qc_blank(metadata[[field]]))) stop("FAILED_INPUT: missing or blank metadata field: ",field)
  if (any(as.character(metadata$database) != gse)) stop("FAILED_INPUT: database metadata differs from GSE")
  counts
}

qc_patterns <- function(species) {
  if (species == "human") c(pMT="^MT-",pRP="^RP[SL]",pHB="^HB[^(P)]")
  else c(pMT="^mt-",pRP="^Rp[sl]",pHB="^Hb[^(p)]")
}
qc_add_percentages <- function(object, species) {
  for (metric in names(qc_patterns(species))) {
    if (!metric %in% names(object@meta.data)) {
      object <- Seurat::PercentageFeatureSet(object,pattern=qc_patterns(species)[[metric]],col.name=metric)
    } else if (any(!is.finite(object@meta.data[[metric]]) | object@meta.data[[metric]] < 0 | object@meta.data[[metric]] > 100)) {
      stop("FAILED_INPUT: invalid existing QC percentage: ",metric)
    }
  }
  object
}

qc_check_dependencies <- function(report, package_available=function(p) requireNamespace(p,quietly=TRUE)) {
  missing <- qc_required_packages[!vapply(qc_required_packages,package_available,logical(1))]
  for (p in missing) qc_report_add(report,"RUN",metric="missing_package",value=p,status="FAILED_INPUT",
    note=paste("R:",R.version.string,"; .libPaths():",paste(.libPaths(),collapse="; ")))
  missing
}

qc_report_config <- function(report,cfg) {
  qc_report_add(report,"CONFIG",metric="species",value=cfg$species,status="USED")
  for (field in names(cfg$cell_qc)) qc_report_add(report,"CONFIG",metric=field,value=cfg$cell_qc[[field]],status="USED")
  for (field in names(cfg$doublet)) qc_report_add(report,"CONFIG",metric=field,
    value=paste(cfg$doublet[[field]],collapse=";"),status="USED",
    note=if(field=="expected_rate") "Expected doublet proportion, not a direct measurement" else "")
  qc_report_add(report,"CELL_CYCLE",metric="regression",value="NONE",status="USED",note="Cell-cycle scoring only; no regression performed")
}
