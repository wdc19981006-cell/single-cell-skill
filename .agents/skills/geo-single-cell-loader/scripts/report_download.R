# Independently report an existing validated raw artifact without discovery,
# downloads, manifest changes, or another Seurat build.
argv <- commandArgs(trailingOnly=TRUE)
script <- sub("^--file=", "", commandArgs()[grepl("^--file=",commandArgs())][1])
scripts <- dirname(normalizePath(script,winslash="/",mustWork=TRUE))
source(file.path(scripts,"seurat_common.R"))
source(file.path(scripts,"download_outputs.R"))
if (length(argv)!=3L || !identical(argv[1],"--existing")) stop("Usage: run_r45.py report_download.R --existing ROOT GSE")
root <- normalizePath(argv[2],winslash="/",mustWork=TRUE)
existing <- load_existing_raw(root,argv[3])
out <- file.path(root,existing$base)
log_path <- file.path(out,"download.log")
if (file.exists(log_path) && !identical(tolower(normalizePath(log_path,winslash="/")),tolower(log_path)))
  stop("Download log redirects outside its declared path")
write_download_outputs(existing$object,out,info=existing$info)
cat(paste0(format(Sys.time(),"%Y-%m-%d %H:%M:%S %Z"),
  " Existing raw artifact independently validated; stage-one CSV/PDF report generated or verified; no new download.\n"),
  file=log_path,append=TRUE)
cat("Stage-one artifacts generated or verified:",out,"\n")
