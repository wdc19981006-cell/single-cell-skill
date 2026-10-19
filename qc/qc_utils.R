if (!isTRUE(getOption("r45.runtime.active"))) stop("Use runtime/r45/run_r45.py")
qc_source_files <- unlist(lapply(sys.frames(),function(frame) frame$ofile),use.names=FALSE)
qc_script_dir <- dirname(normalizePath(tail(qc_source_files,1L),winslash="/",mustWork=TRUE))
source(file.path(dirname(qc_script_dir),"runtime/project_paths.R"))
qc_required_packages <- c("Seurat", "SeuratObject", "Matrix", "yaml", "jsonlite", "ggplot2", "patchwork", "DoubletFinder")
qc_blank <- function(x) is.na(x) | !nzchar(trimws(as.character(x)))

qc_paths <- function(root, gse) {
  sc_dataset_name(gse)
  root <- normalizePath(root, winslash="/", mustWork=TRUE)
  dataset <- sc_project_path(root,gse)
  legacy <- file.path(root,"data",gse)
  canonical_raw <- file.path(dataset,"01_download","seurat_raw.rds")
  raw <- if(file.exists(canonical_raw)) canonical_raw else file.path(legacy,"seurat_raw.rds")
  qc <- file.path(dataset,"02_QC")
  dir.create(qc,recursive=TRUE,showWarnings=FALSE)
  if(!identical(normalizePath(qc,winslash="/"),qc)) stop("QC directory redirected")
  decision <- file.path(qc,"qc_decision.json")
  if(!file.exists(decision) && file.exists(file.path(legacy,"qc","qc_decision.json")))
    decision <- file.path(legacy,"qc","qc_decision.json")
  dirs <- c(qc1="01_basic_QC",qc2="02_doublet_cellcycle",qc3="03_ambient_RNA")
  list(root=root,gse=gse,dataset=dataset,raw=raw,qc=qc,
       report=file.path(qc,"qc_report.csv"),decision=decision,
       stage_dirs=stats::setNames(file.path(qc,dirs),names(dirs)),
       final=file.path(qc,dirs[["qc2"]],"seurat_QC-2.rds"),
       legacy_final=file.path(legacy,"qc","seurat_qc.rds"))
}

qc_stage_paths <- function(paths,step) {
  if(!step %in% names(paths$stage_dirs)) stop("Unknown QC step")
  paths$qc <- paths$stage_dirs[[step]]
  dir.create(paths$qc,recursive=TRUE,showWarnings=FALSE)
  if(!identical(normalizePath(paths$qc,winslash="/"),paths$qc)) stop("QC step directory redirected")
  paths$report <- file.path(paths$qc,"qc_report.csv")
  paths$final <- file.path(paths$qc,paste0("seurat_QC-",match(step,c("qc1","qc2","qc3")),".rds"))
  paths
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
  # This R 4.5.0 host can double-convert strings marked UTF-8 when writing
  # through the native locale. All workflow text arrives as UTF-8 bytes;
  # remove only the encoding flag so write.csv emits those bytes unchanged.
  for (field in names(rows)) Encoding(rows[[field]]) <- "unknown"
  utils::write.csv(rows,temporary,row.names=FALSE,na="")
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
  limit <- cfg$cell_qc$max_predicted_qc_removal_fraction
  if (!is.numeric(limit) || length(limit)!=1L || !is.finite(limit) || limit<0 || limit>=1)
    stop("Invalid predicted QC removal fraction limit")
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
      any(!types %in% c("remove_low_cell_samples","keep_without_doubletfinder","override_parameters","accept_default_conflicts","accept_qc_removal")))
    stop("QC decision sample list or type invalid")
  if (length(c(removed,kept)) &&
      !(if ("accept_qc_removal" %in% types) all(c(removed,kept) %in% affected)
        else setequal(c(removed,kept), affected)))
    stop("QC affected_samples must include removed/kept samples; match them exactly without a separate removal-acceptance decision")
  if (length(removed) && !"remove_low_cell_samples" %in% d$decision_type) stop("Removal decision_type missing")
  if (length(kept) && !"keep_without_doubletfinder" %in% d$decision_type) stop("Keep decision_type missing")
  if (!is.null(d$parameter_overrides) && !"override_parameters" %in% d$decision_type) stop("Override decision_type missing")
  if (!is.null(d$accept_qc_removal) && (!identical(d$accept_qc_removal, TRUE) ||
      !"accept_qc_removal" %in% types)) stop("accept_qc_removal requires an explicit user decision")
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
  for (field in c("nFeature_RNA","nCount_RNA"))
    if (!field %in% names(metadata) || !is.numeric(metadata[[field]]) ||
        any(!is.finite(metadata[[field]]) | metadata[[field]]<0))
      stop("FAILED_INPUT: missing or invalid QC metric: ",field)
  if (grepl("^GSE[0-9]+$",gse) && any(as.character(metadata$database) != gse)) stop("FAILED_INPUT: database metadata differs from GSE")
  counts
}

qc_refresh_metrics <- function(object,preserve=TRUE) {
  counts <- qc_counts(object)
  current <- list(nCount_RNA=as.numeric(Matrix::colSums(counts)),nFeature_RNA=as.numeric(Matrix::colSums(counts>0)))
  for(field in c("nCount_RNA","nFeature_RNA")) {
    original <- paste0(field,"_before_QC")
    metadata_field <- paste0(field,"_input_metadata")
    if(preserve && !metadata_field %in% names(object@meta.data)) object[[metadata_field]] <- object@meta.data[[field]]
    if(preserve && !original %in% names(object@meta.data)) object[[original]] <- current[[field]]
  }
  object$nCount_RNA <- as.numeric(Matrix::colSums(counts))
  object$nFeature_RNA <- as.numeric(Matrix::colSums(counts>0))
  object
}

qc_capture_libraries <- function(object) {
  meta <- object@meta.data
  if("library_id" %in% names(meta)) {
    if(any(qc_blank(meta$library_id))) stop("AMBIGUOUS_CAPTURE_LIBRARY: library_id is incomplete; provide an explicit cell-to-capture-library mapping")
    object$library_id <- as.character(meta$library_id)
    if(!"capture_mapping_source" %in% names(meta)) object$capture_mapping_source <- "library_id"
    else if(any(qc_blank(meta$capture_mapping_source))) stop("AMBIGUOUS_CAPTURE_LIBRARY: capture_mapping_source is incomplete")
    return(object)
  }
  if("capture_id" %in% names(meta) && !any(qc_blank(meta$capture_id))) {
    object$library_id <- as.character(meta$capture_id)
    object$capture_mapping_source <- "capture_id"
    return(object)
  }
  clues <- intersect(c("source_type","technology","capture_layout","pool_id"),names(meta))
  pooled <- length(clues) && any(grepl("pool|multiplex|demultiplex|multiple.capture",unlist(meta[clues]),ignore.case=TRUE))
  runs <- if("sequencing_run" %in% names(meta)) tapply(as.character(meta$sequencing_run),meta$sample,function(x) length(unique(x[!qc_blank(x)]))) else 1L
  if(pooled || any(runs>1L)) stop("AMBIGUOUS_CAPTURE_LIBRARY: pooled/multiple-run metadata cannot establish capture libraries; supply library_id before DoubletFinder")
  object$library_id <- as.character(meta$sample)
  object$capture_mapping_source <- "legacy_sample_mapping_unverified"
  object
}

qc_patterns <- function(species) {
  if (species == "human") c(pMT="^MT-",pRP="^RP[SL]",pHB="^HB[^(P)]")
  else c(pMT="^mt-",pRP="^Rp[sl]",pHB="^Hb[^(p)]")
}
qc_add_percentages <- function(object, species) {
  object <- qc_refresh_metrics(object)
  for (metric in names(qc_patterns(species))) {
    original <- paste0(metric,"_before_QC")
    metadata_field <- paste0(metric,"_input_metadata")
    if(metric %in% names(object@meta.data) && !metadata_field %in% names(object@meta.data)) object[[metadata_field]] <- object@meta.data[[metric]]
    object <- Seurat::PercentageFeatureSet(object,pattern=qc_patterns(species)[[metric]],col.name=metric)
    if(!original %in% names(object@meta.data)) object[[original]] <- object@meta.data[[metric]]
  }
  object
}

qc_cell_filter_rules <- function(meta, thresholds) {
  list(nFeature_min=meta$nFeature_RNA <= thresholds$nFeature_min,
       nFeature_max=meta$nFeature_RNA >= thresholds$nFeature_max,
       nCount_min=meta$nCount_RNA <= thresholds$nCount_min,
       nCount_max=meta$nCount_RNA >= thresholds$nCount_max,
       pMT_max=meta$pMT >= thresholds$pMT_max,
       pHB_max=meta$pHB >= thresholds$pHB_max,
       pRP_max=meta$pRP >= thresholds$pRP_max)
}
qc_cell_keep <- function(meta, thresholds) !Reduce(`|`,qc_cell_filter_rules(meta,thresholds))

# Rule hits overlap. Only total_removed_cells is the deduplicated union.
qc_filter_reasons <- function(meta,thresholds,removed_samples=character()) {
  rules <- as.data.frame(qc_cell_filter_rules(meta,thresholds))
  hits <- rowSums(rules)
  do.call(rbind,lapply(unique(as.character(meta$sample)),function(sample) {
    selected <- meta$sample==sample
    excluded <- sample %in% removed_samples
    row <- data.frame(sample=sample,input_cells=sum(selected),stringsAsFactors=FALSE)
    for(rule in names(rules)) {
      row[[paste0(rule,"_hits")]] <- sum(rules[[rule]][selected])
      row[[paste0(rule,"_only")]] <- sum(rules[[rule]][selected] & hits[selected]==1L)
    }
    row$multiple_conditions <- sum(hits[selected]>1L)
    row$qc_removed_cells <- sum(hits[selected]>0L)
    row$user_removed_sample_cells <- if(excluded) sum(selected) else 0L
    row$total_removed_cells <- if(excluded) sum(selected) else row$qc_removed_cells
    row$removed_fraction <- row$total_removed_cells/row$input_cells
    row
  }))
}

qc_cycle_features <- function(object) {
  reference <- Seurat::cc.genes
  available <- rownames(object)
  matched <- function(names) {
    indexes <- match(toupper(names),toupper(available))
    unique(available[stats::na.omit(indexes)])
  }
  list(s=matched(reference$s.genes),g2m=matched(reference$g2m.genes))
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
