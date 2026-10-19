# Actual celda DecontX and non-overlapping QC reason accounting on a small fixture.
args <- commandArgs(trailingOnly=TRUE)
repo <- normalizePath(args[1],winslash="/",mustWork=TRUE)
source(file.path(repo,"qc/run_qc.R"))
source(file.path(repo,"analysis/run_analysis.R"))
out <- args[3]; dir.create(out,recursive=TRUE,showWarnings=FALSE)
object <- readRDS(args[2])
SeuratObject::DefaultAssay(object) <- "RNA"
if("decontX" %in% names(object@assays)) object[["decontX"]] <- NULL
for(field in intersect(c("nCount_decontX","nFeature_decontX"),names(object@meta.data))) object[[field]] <- NULL
before <- qc_counts(object); meta_before <- object@meta.data
thresholds <- yaml::read_yaml(file.path(repo,"qc/default_config.yaml"))$cell_qc
reasons_meta <- data.frame(sample=c("A","A","A","B"),nFeature_RNA=c(200,201,4000,201),
  nCount_RNA=c(500,501,30000,501),pMT=c(25,1,1,1),pHB=0,pRP=1)
reason <- qc_filter_reasons(reasons_meta,thresholds)
stopifnot(reason$input_cells[1]==3,reason$total_removed_cells[1]==2,reason$multiple_conditions[1]==2,
  reason$nFeature_min_hits[1]==1,reason$nFeature_min_only[1]==0)
qc_write_filter_reasons(reasons_meta,list(cell_qc=thresholds),list(qc=out))
report <- qc_report_new()
corrected <- qc_run_ambient(object,list(qc=out),report)
qc_report_write(report,file.path(out,"qc_report.csv"))
stopifnot(identical(qc_counts(corrected),before),identical(colnames(corrected),colnames(object)),
  identical(corrected@meta.data[,names(meta_before),drop=FALSE],meta_before))
matrix <- analysis_counts(corrected,"decontX")
stopifnot(any(matrix@x!=round(matrix@x)),all(is.finite(matrix@x)))
path <- file.path(out,"synthetic_actual_decontX.rds"); saveRDS(corrected,path)
stopifnot(identical(analysis_counts(readRDS(path),"decontX"),matrix))
# Diagnose an extreme estimate without filtering and without special accession rules.
extreme <- corrected; extreme$decontX_contamination <- rep(.99,ncol(extreme))
extreme_out <- file.path(out,"extreme_diagnostic_fixture"); dir.create(extreme_out)
review <- qc_review_ambient(extreme,list(qc=extreme_out),qc_report_new(),list(markers=c("g1","g50","absent")))
stopifnot(any(review$warnings$category=="EXTREME_CONTAMINATION"),ncol(extreme)==ncol(object))
jsonlite::write_json(list(status="PASS",method="actual installed celda::decontX",cells=ncol(corrected),
  fractional_counts_preserved=TRUE,original_RNA_and_metadata_preserved=TRUE,
  extreme_estimate_flagged_without_deletion=TRUE),file.path(out,"test_results.json"),auto_unbox=TRUE,pretty=TRUE)
cat("PASS: actual DecontX, fractional Seurat5 counts, unchanged cells/metadata, QC reason union and extreme diagnostic\n")
