# Diagnose a saved QC-3; input RDS and historical figures remain untouched.
file_arg <- sub("^--file=","",commandArgs()[grepl("^--file=",commandArgs())][1])
directory <- dirname(normalizePath(file_arg,winslash="/",mustWork=TRUE))
source(file.path(directory,"run_qc.R"))
args <- commandArgs(trailingOnly=TRUE)
if(!length(args) %in% c(2L,3L)) stop("Usage: review_ambient.R INPUT_RDS NEW_OUTPUT [CONFIG_YAML]")
cfg <- if(length(args)==3L) yaml::read_yaml(args[3]) else NULL
object <- readRDS(args[1])
report <- qc_report_new()
review <- qc_review_ambient(object,list(qc=args[2]),report,cfg)
qc_report_write(report,file.path(args[2],"qc_report.csv"))
jsonlite::write_json(list(stage="ambient-review",status="COMPLETE",input=args[1],
  outputs=list.files(args[2],full.names=TRUE),parameters=qc_ambient_config(cfg),statistics=review$summary,
  unresolved_issues=review$warnings,needs_user_decision=FALSE,
  next_command="Review warnings and markers; no automatic deletion or QC recomputation"),
  file.path(args[2],"ambient_review_result.json"),auto_unbox=TRUE,pretty=TRUE)
cat("PASS: read-only ambient diagnostic review\n")
