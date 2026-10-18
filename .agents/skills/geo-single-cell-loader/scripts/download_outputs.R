# Stage-one user artifacts, derived only from a validated raw Seurat object.
immutable_download_csv <- function(value,path) {
  if (file.exists(path) && !identical(tolower(normalizePath(path,winslash="/")),tolower(path)))
    stop("Download artifact redirects outside its declared path: ",path)
  pending <- tempfile("download-csv-",tmpdir=dirname(path))
  on.exit(unlink(pending),add=TRUE)
  write.csv(value,pending,row.names=FALSE,na="",fileEncoding="UTF-8")
  if (file.exists(path)) {
    if (!identical(unname(tools::md5sum(path)),unname(tools::md5sum(pending))))
      stop("Existing download artifact differs; refusing overwrite: ",path)
  } else if (!file.rename(pending,path)) stop("Cannot publish download artifact: ",path)
}
write_download_outputs <- function(object,out,manifest=NULL,info=character()) {
  validate_object(object,manifest)
  counts <- if (inherits(object[["RNA"]],"Assay5"))
    SeuratObject::LayerData(object,assay="RNA",layer="counts") else
    SeuratObject::GetAssayData(object,assay="RNA",slot="counts")
  metadata <- data.frame(cell=colnames(object),object@meta.data,check.names=FALSE)
  immutable_download_csv(metadata,file.path(out,"metadata.csv"))
  samples <- unique(as.character(object$sample))
  summary <- do.call(rbind,lapply(samples,function(sample) {
    cells <- which(object$sample==sample)
    source <- if (!is.null(manifest)) manifest[match(sample,manifest$sample),,drop=FALSE] else NULL
    # Archived manifests remain in the private audit. The original TXT retains
    # actual reader records, which can be reused without changing any receipt.
    records <- strsplit(info[startsWith(info,paste0(sample,"\t"))],"\t",fixed=TRUE)
    record <- if (length(records)) tail(records,1)[[1]] else character()
    data.frame(database=unique(as.character(object$database[cells])),sample=sample,
      group=unique(as.character(object$group[cells])),cells=length(cells),
      genes_detected=sum(Matrix::rowSums(counts[,cells,drop=FALSE]>0)>0),
      object_features=nrow(counts),total_counts=sum(counts[,cells,drop=FALSE]),
      file_type=if (!is.null(source)) source$file_type else if (length(record)>=2L) record[2] else NA_character_,
      local_path=if (!is.null(source)) source$local_path else if (length(record)>=3L) record[3] else NA_character_,
      count_source=if (!is.null(source)) source$count_source else if (length(record)>=4L) sub("^count_source=","",record[4]) else NA_character_,
      structure=paste(class(counts),collapse="/"),
      integrity="validated sparse finite nonnegative integer RNA counts; confirmed metadata and cell alignment",
      stringsAsFactors=FALSE)
  }))
  immutable_download_csv(summary,file.path(out,"download_summary.csv"))
  report <- file.path(out,"download_report.pdf")
  if (file.exists(report) && !identical(tolower(normalizePath(report,winslash="/")),tolower(report)))
    stop("Download report redirects outside its declared path: ",report)
  if (!file.exists(report)) {
    pending <- tempfile("download-report-",tmpdir=out,fileext=".pdf")
    on.exit(unlink(pending),add=TRUE)
    grDevices::pdf(pending,width=8.3,height=11.7,onefile=TRUE)
    device <- grDevices::dev.cur()
    on.exit(if (device %in% grDevices::dev.list()) grDevices::dev.off(device),add=TRUE)
    graphics::par(mar=c(4,5,4,2),family="sans",cex=0.85)
    graphics::plot.new()
    lines <- c("Download and raw-object integrity report",paste("Dataset:",unique(object$database)),
      paste("Validated cells:",ncol(object)),paste("Validated features:",nrow(object)),
      paste("Samples:",length(samples)),paste("Sparse counts structure:",paste(class(counts),collapse="/")),
      "Counts: finite, nonnegative integers; one RNA counts layer.",
      "Cell keys, sample prefixes, user-confirmed groups and dimensions: PASS.",
      "nCount_RNA and nFeature_RNA match the current counts matrix.",
      "Construction thresholds: min.cells=3; min.features=200.",
      "No additional QC, normalization, integration or clustering was run.",
      "Genes detected per sample are measured within the constructed object.",
      "Small samples are reported and retained; no sample is removed here.",
      "Source expression inputs remain unchanged in raw/.",
      paste("R:",getRversion()),paste("Seurat:",packageVersion("Seurat")))
    graphics::text(0,seq(0.98,0.1,length.out=length(lines)),labels=lines,adj=c(0,1),cex=0.9)
    for (indices in split(seq_len(nrow(summary)),ceiling(seq_len(nrow(summary))/25))) {
      page <- summary[indices,,drop=FALSE]
      graphics::par(mfrow=c(2,1),mar=c(8,5,4,2),cex=0.75)
      graphics::barplot(page$cells,names.arg=page$sample,las=2,col="#3B82A0",ylab="Cells",main="Validated cells by sample")
      graphics::barplot(page$genes_detected,names.arg=page$sample,las=2,col="#659868",ylab="Detected genes",main="Detected genes by sample")
    }
    grDevices::dev.off(device)
    if (file.exists(report) || !file.rename(pending,report)) stop("Cannot publish download report without overwriting: ",report)
  }
  invisible(summary)
}
