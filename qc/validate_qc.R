# Read-only validation: never rebuilds or replaces an RDS.
args <- commandArgs(trailingOnly=TRUE)
script <- sub("^--file=","",commandArgs()[grepl("^--file=",commandArgs())][1])
source(file.path(dirname(normalizePath(script)),"run_qc.R"))
if(!length(args) %in% 2:4) stop("Usage: run_r45.py qc/validate_qc.R ROOT GSE [qc1|qc2|qc3|legacy] [INPUT]")
base <- qc_paths(args[1],args[2]); step <- if(length(args)>2L) args[3] else "qc2"
if(step=="legacy") {
  object <- readRDS(base$legacy_final)
  qc_stage_validate(object,"qc2",args[2],legacy=TRUE)
  path <- file.path(dirname(base$legacy_final),"qc_report.csv")
} else {
  paths <- qc_stage_paths(base,step)
  object <- readRDS(paths$final)
  input <- if(length(args)>3L) args[4] else object@misc$qc_pipeline[[step]]$input
  qc_stage_validate(object,step,args[2],input,paths)
  path <- paths$report
}
report <- utils::read.csv(path,stringsAsFactors=FALSE)
status <- tail(report$value[report$section=="RUN" & report$metric=="status"],1)
if(!status %in% c("COMPLETE_QC1","COMPLETE_QC","QC_COMPLETE_WITH_UNEVALUATED_DOUBLETS","COMPLETE_QC3")) stop("QC report incomplete")
for(metric in c("final_cells","final_genes")) {
  value <- tail(report$value[report$section=="RUN" & report$metric==metric],1)
  dimension <- if(metric=="final_cells") ncol(object) else nrow(object)
  if(length(value)!=1L || as.integer(value)!=dimension) stop("QC report dimensions disagree")
}
cat("Read-only QC artifact validation:",step,status,"\n")
