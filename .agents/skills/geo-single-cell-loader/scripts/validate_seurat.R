# Invoke through runtime/r45/run_r45.py.
argv <- commandArgs(trailingOnly=TRUE)
script <- sub("^--file=", "", commandArgs()[grepl("^--file=",commandArgs())][1])
source(file.path(dirname(normalizePath(script)),"seurat_common.R"))
if (length(argv)==3L && identical(argv[1],"--existing")) {
  # Read-only regression after audit cleanup. Full manifest mode below remains
  # mandatory for a new build. Never infer a new grouping from the object.
  root <- normalizePath(argv[2],winslash="/",mustWork=TRUE)
  gse <- argv[3]
  if (!grepl("^GSE[0-9]+$",gse)) stop("Expected GSE accession")
  base <- paste0("data/",gse)
  need("jsonlite")
  approval <- jsonlite::fromJSON(repo_path(root,paste0(base,"/group_confirmation.json"),base))
  if (!identical(approval$confirmed_by,"user") || blank(approval$user_statement)) stop("User confirmation missing")
  info <- readLines(repo_path(root,paste0(base,"/sample_info.txt"),base),encoding="UTF-8",warn=FALSE)
  if (!identical(info[1],"STATUS: COMPLETE_RAW")) stop("Archived raw object is not recorded COMPLETE_RAW")
  object <- readRDS(repo_path(root,paste0(base,"/seurat_raw.rds"),base))
  validate_object(object)
  expected_groups <- unname(as.character(unlist(approval$groups[object$sample])))
  if (!setequal(unique(object$sample),names(approval$groups)) ||
      !identical(unname(as.character(object$group)),expected_groups))
    stop("Archived object differs from confirmed groups")
  if (!identical(unique(as.character(object$database)),gse)) stop("Dataset identity mismatch")
  expected <- function(label) {
    position <- tail(which(info==label),1)
    if (!length(position)) stop("Archived summary lacks ",label)
    as.integer(info[position+1L])
  }
  if (ncol(object)!=expected("Total cells:") || nrow(object)!=expected("Total genes:")) stop("Archived dimensions differ")
  cat("Existing raw artifact validated read-only:",ncol(object),"cells;",nrow(object),"genes\n")
  quit(save="no",status=0L)
}
record_completion <- length(argv)==4L && identical(argv[4],"--record-completion")
if (length(argv)!=3L && !record_completion) stop("Usage: run_r45.py validate_seurat.R ROOT MANIFEST RDS [--record-completion], or --existing ROOT GSE")
root <- normalizePath(argv[1],winslash="/",mustWork=TRUE)
m <- read_manifest(repo_path(root,argv[2],"data"),root)
expected <- paste0("data/",m$database[1],"/seurat_raw.rds")
if (!identical(argv[3],expected)) stop("Expected final data/<GSE>/seurat_raw.rds")
object <- readRDS(repo_path(root,argv[3],paste0("data/",m$database[1])))
# Shared validation enforces one sparse integer RNA counts layer, exact manifest metadata,
# and absence of normalized/downstream assays, reductions, graphs, neighbors and extra metadata.
validate_object(object,m)
if (record_completion) {
  out <- file.path(root,"data",m$database[1])
  info <- file.path(out,"sample_info.txt")
  if (!file.exists(info) || !identical(readLines(info,n=1L,warn=FALSE),"STATUS: COMPLETE_RAW")) {
    source(file.path(dirname(normalizePath(script)),"sample_info.R"))
    update_sample_info(out,"COMPLETE_RAW",m,object,
      inputs=paste(m$sample,m$file_type,m$local_path,paste0("count_source=",m$count_source),sep="\t"),
      warnings="Completion recovered by independent validation of the existing artifact.")
  }
  retained <- table(object$sample)
  write.csv(data.frame(sample=names(retained),cells=as.integer(retained),
    group=m$group[match(names(retained),m$sample)]),file.path(out,".workflow","sample_summary.csv"),row.names=FALSE)
}
cat("Standard raw Seurat counts, metadata, cell alignment and clean structure validated\n")
