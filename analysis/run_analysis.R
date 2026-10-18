# All analysis computation must be invoked through the fixed R 4.5.0 launcher.
analysis_source <- if(sys.nframe()>=1L && !is.null(sys.frame(1)$ofile)) sys.frame(1)$ofile else
  sub("^--file=","",commandArgs()[grepl("^--file=",commandArgs())][1])
analysis_runner_dir <- dirname(normalizePath(analysis_source,winslash="/",mustWork=TRUE))
for(module in c("analysis_utils.R","validate_analysis.R","pca_regression.R","harmony.R","clustering.R"))
  source(file.path(analysis_runner_dir,module))

if(sys.nframe()==0L) {
  args <- commandArgs(trailingOnly=TRUE)
  if(length(args)!=4L) stop("Usage: run_analysis.R STEP INPUT_RDS OUTPUT_DIRECTORY REQUEST_JSON")
  request <- jsonlite::fromJSON(args[4],simplifyVector=FALSE)
  cfg <- analysis_config(request$configuration)
  step <- args[1]
  if(!dir.exists(args[3])) stop("Analysis output directory missing")
  switch(step,
    pca=analysis_run_pca(args[2],args[3],cfg,request$provenance),
    harmony=analysis_run_harmony(args[2],args[3],cfg,request$provenance),
    `cluster-review`=analysis_run_cluster_review(args[2],args[3],cfg,request$provenance),
    finalize=analysis_finalize(args[2],args[3],cfg,request$provenance,request$decision),
    stop("Unknown analysis step"))
  cat("PASS: analysis stage ",step,"\n",sep="")
}
