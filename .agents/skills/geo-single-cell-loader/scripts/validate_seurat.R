# Invoke through runtime/r45/run_r45.py.
argv <- commandArgs(trailingOnly=TRUE)
script <- sub("^--file=", "", commandArgs()[grepl("^--file=",commandArgs())][1])
source(file.path(dirname(normalizePath(script)),"seurat_common.R"))
if (length(argv)==3L && identical(argv[1],"--existing")) {
  # Read-only regression after audit cleanup. Full manifest mode below remains
  # mandatory for a new build. Never infer a new grouping from the object.
  root <- normalizePath(argv[2],winslash="/",mustWork=TRUE)
  gse <- argv[3]
  object <- load_existing_raw(root,gse)$object
  cat("Existing raw artifact validated read-only:",ncol(object),"cells;",nrow(object),"genes\n")
  quit(save="no",status=0L)
}
record_completion <- length(argv)==4L && identical(argv[4],"--record-completion")
if (length(argv)!=3L && !record_completion) stop("Usage: run_r45.py validate_seurat.R ROOT MANIFEST RDS [--record-completion], or --existing ROOT GSE")
root <- normalizePath(argv[1],winslash="/",mustWork=TRUE)
base <- manifest_dataset(argv[2])
m <- read_manifest(repo_path(root,argv[2],paste0(base,"/.workflow")),root)
if (!identical(base,dataset_relative(root,m$database[1]))) stop("Manifest is outside the active dataset download root")
expected <- paste0(base,"/seurat_raw.rds")
if (!identical(argv[3],expected)) stop("Expected final seurat_raw.rds at the manifest's dataset root")
object <- readRDS(repo_path(root,argv[3],base))
# Shared validation enforces one sparse integer RNA counts layer, exact manifest metadata,
# and absence of normalized/downstream assays, reductions, graphs, neighbors and extra metadata.
validate_object(object,m)
if (record_completion) {
  out <- file.path(root,base)
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
  source(file.path(dirname(normalizePath(script)),"download_outputs.R"))
  write_download_outputs(object,out,m,readLines(info,encoding="UTF-8",warn=FALSE))
}
cat("Standard raw Seurat counts, metadata, cell alignment and clean structure validated\n")
