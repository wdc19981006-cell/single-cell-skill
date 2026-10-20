# Fixture regression. Real installed DoubletFinder is exercised in test_qc_real_df.R.
args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if(length(args)) args[1] else ".",winslash="/",mustWork=TRUE)
source(file.path(root,"qc/run_qc.R"))
set.seed(716)
area <- tempfile("qc-stages-",tmpdir=file.path(root,".runtime")); dir.create(area,recursive=TRUE)
area <- normalizePath(area,winslash="/")
gse <- "GSE999999991"
raw_dir <- file.path(area,"results",gse,"01_download"); dir.create(raw_dir,recursive=TRUE)
genes <- unique(c("MT-TEST","RPS-TEST","HBA1",Seurat::cc.genes$s.genes,Seurat::cc.genes$g2m.genes,paste0("GENE",1:5000),"RARE3","RARE4"))
cells <- paste0("cell",1:160)
counts <- Matrix::Matrix(matrix(stats::rpois(length(genes)*length(cells),.5),nrow=length(genes),dimnames=list(genes,cells)),sparse=TRUE)
counts["HBA1",] <- 1; counts["RARE3",] <- 0; counts["RARE3",1:3] <- 1
counts["RARE4",] <- 0; counts["RARE4",1:4] <- 1
meta <- data.frame(sample=rep(c("A","B"),each=80),database=gse,group=rep(c("Tumor","Control"),each=80),library_id="CAPTURE1",row.names=cells)
object <- Seurat::CreateSeuratObject(counts=counts,min.cells=0,min.features=0,meta.data=meta)
# Stale input summaries must be preserved separately and recomputed from counts.
object$nCount_RNA <- 1
raw <- file.path(raw_dir,"seurat_raw.rds"); saveRDS(object,raw)
raw_md5 <- qc_md5(raw)
paths <- qc_paths(area,gse)
stopifnot(grepl("results/.*/02_QC$",paths$qc),identical(paths$raw,raw))
stopifnot(qc_precheck(area,gse,step="qc1")=="PASS",qc_precheck(area,gse,step="qc2")=="PASS")
boundary <- data.frame(nFeature_RNA=c(200,201,4000,3999,300),
  nCount_RNA=c(600,500,600,30000,501),pMT=c(1,1,25,1,24.9),pHB=c(0,0,0,1,.9),pRP=c(1,1,1,100,99.9))
stopifnot(identical(qc_cell_keep(boundary,qc_config(paths)$cell_qc),c(FALSE,FALSE,FALSE,FALSE,TRUE)))
# Consent schemas keep real sample IDs; an ALL summary cannot grant permission.
jsonlite::write_json(list(decision_type="accept_qc_removal",affected_samples="ALL",confirmed_by="user",
  user_statement="SIMULATED FIXTURE DECISION",accept_qc_removal=TRUE),paths$decision,auto_unbox=TRUE)
stopifnot(inherits(try(qc_decision(paths,c("A","B")),silent=TRUE),"try-error"))
unlink(paths$decision)
mock_api <- list(paramSweep=function(seu,PCs,sct) {
    stopifnot(length(unique(seu$library_id))==1L,length(unique(seu$sample))==2L)
    list(PCs=PCs)
  },summarizeSweep=function(sweep,GT) sweep,
  find.pK=function(stats) data.frame(pK=c("0.05","0.1"),BCmetric=c(1,2)),
  modelHomotypic=function(clusters) .2,
  doubletFinder=function(seu,PCs,pN,pK,nExp,reuse.pANN,sct) {
    if(is.null(reuse.pANN)) seu[[paste0("pANN_",nExp)]] <- seq_len(ncol(seu))/ncol(seu)
    calls <- rep("Singlet",ncol(seu)); calls[tail(seq_len(ncol(seu)),nExp)] <- "Doublet"
    seu[[paste0("DF.classifications_",nExp)]] <- calls
    seu
  })
stopifnot(qc_run(area,gse,step="qc1")=="COMPLETE_QC1")
qc1_path <- qc_stage_paths(paths,"qc1")$final
qc1 <- readRDS(qc1_path)
stopifnot(!"RARE3" %in% rownames(qc1),"RARE4" %in% rownames(qc1),
  all(qc1$nCount_RNA_input_metadata==1),all(qc1$nCount_RNA>500),
  identical(qc_counts(qc1),counts[rownames(qc1),colnames(qc1),drop=FALSE]))
before <- qc_md5(qc1_path)
stopifnot(qc_run(area,gse,step="qc1")=="COMPLETE_QC1",identical(before,qc_md5(qc1_path)))
stopifnot(qc_run(area,gse,step="qc2",df_api=mock_api)=="COMPLETE_QC")
qc2_path <- qc_stage_paths(paths,"qc2")$final
qc2 <- readRDS(qc2_path)
stopifnot(all(qc2$DF_adj=="Singlet"),any(qc2$DF=="Doublet"),ncol(qc2)==155L,
  identical(qc_counts(qc2),qc_counts(qc1)[,colnames(qc2),drop=FALSE]),
  all(c("S.Score","G2M.Score","Phase") %in% names(qc2@meta.data)))
# Exercise actual small-feature cycle scoring without repeating spatial rendering.
old_prepare <- qc_prepare_cycle_plot_object; old_plot <- qc_plot_cycle
qc_prepare_cycle_plot_object <- identity; qc_plot_cycle <- function(...) invisible(NULL)
small_paths <- qc_stage_paths(paths,"qc2"); small_paths$qc <- file.path(area,"small_cycle_test"); dir.create(small_paths$qc)
small_cycle <- qc_score_cycle(qc2[seq_len(800),],small_paths,qc_report_new())
qc_prepare_cycle_plot_object <- old_prepare; qc_plot_cycle <- old_plot
stopifnot(ncol(small_cycle)==ncol(qc2),all(is.finite(small_cycle$S.Score)),all(is.finite(small_cycle$G2M.Score)))
mock_ambient <- list(run=function(sce,batch,seed) {
    SummarizedExperiment::colData(sce)$decontX_contamination <- rep(.25,ncol(sce)); sce
  },counts=function(sce) SummarizedExperiment::assay(sce,"counts")*.75)
stopifnot(qc_run(area,gse,step="qc3",ambient_api=mock_ambient)=="COMPLETE_QC3")
qc3 <- readRDS(qc_stage_paths(paths,"qc3")$final)
if(Sys.getenv("SC_TEST_COMPACT_QC")=="1") {
  # Reuse the scientific fixture, but validate a separate formal compact copy.
  old_mode <- Sys.getenv("SINGLE_CELL_STORAGE_MODE")
  Sys.setenv(SINGLE_CELL_STORAGE_MODE="compact")
  compact_qc <- sc_compact_formal(qc3,file.path(area,"compact_qc.rds"),"qc")
  saveRDS(compact_qc,file.path(area,"compact_qc.rds"))
  sc_verify_formal(readRDS(file.path(area,"compact_qc.rds")))
  qc_stage_validate(compact_qc,"qc3",gse,qc2_path,qc_stage_paths(paths,"qc3"))
  stopifnot(identical(qc_counts(compact_qc),qc_counts(qc3)),
    identical(SeuratObject::LayerData(compact_qc,assay="decontX",layer="counts"),
              SeuratObject::LayerData(qc3,assay="decontX",layer="counts")))
  if(nzchar(old_mode)) Sys.setenv(SINGLE_CELL_STORAGE_MODE=old_mode) else Sys.unsetenv("SINGLE_CELL_STORAGE_MODE")
}
stopifnot(identical(qc_counts(qc3),qc_counts(qc2)),identical(colnames(qc3),colnames(qc2)),
  all(qc3$decontX_contamination==.25),"decontX" %in% names(qc3@assays),
  identical(SeuratObject::LayerData(qc3,assay="decontX",layer="counts"),qc_counts(qc2)*.75),
  identical(SeuratObject::DefaultAssay(qc3),"RNA"),identical(qc_md5(raw),raw_md5))
for(step in c("qc1","qc2","qc3")) {
  p <- qc_stage_paths(paths,step); x <- readRDS(p$final)
  qc_stage_validate(x,step,gse,x@misc$qc_pipeline[[step]]$input,p)
  stopifnot(file.exists(file.path(p$qc,paste0(step,"_config.yaml"))),file.info(file.path(p$qc,paste0(step,".log")))$size>0)
}
ambiguous <- object; ambiguous$library_id <- NULL; ambiguous$source_type <- "pooled_matrix"
stopifnot(inherits(try(qc_capture_libraries(ambiguous),silent=TRUE),"try-error"))
# Preserve precheck decisions and conflict handling from the historical suite.
scenario <- function(id,source=object) {
  source$database <- id
  directory <- file.path(area,"results",id,"01_download"); dir.create(directory,recursive=TRUE)
  saveRDS(source,file.path(directory,"seurat_raw.rds")); qc_paths(area,id)
}
write_decision <- function(p,value) jsonlite::write_json(value,p$decision,auto_unbox=TRUE)
low <- object[,1:40]; low$library_id <- "LOW_CAPTURE"
p_low <- scenario("GSE999999981",low)
stopifnot(qc_precheck(area,p_low$gse)=="NEEDS_USER_DECISION",
  any(utils::read.csv(p_low$report)$status=="LOW_CELL_SAMPLE"))
write_decision(p_low,list(decision_type="keep_without_doubletfinder",affected_samples="A",confirmed_by="user",
  user_statement="SIMULATED FIXTURE DECISION: keep A without DoubletFinder",kept_low_cell_samples="A"))
stopifnot(qc_precheck(area,p_low$gse)=="PASS")
never_df <- mock_api; never_df$paramSweep <- function(...) stop("Skipped capture must never run DoubletFinder")
stopifnot(qc_run(area,p_low$gse,step="qc2",df_api=never_df)=="QC_COMPLETE_WITH_UNEVALUATED_DOUBLETS")
kept_low <- readRDS(qc_stage_paths(p_low,"qc2")$final)
stopifnot(ncol(kept_low)==40L,all(kept_low$DF=="NotEvaluated"),all(kept_low$DF_adj=="NotEvaluated"),
  all(kept_low$doublet_status=="SKIPPED_LOW_CELL"),all(is.na(kept_low$pANN)))
p_remove <- scenario("GSE999999982")
write_decision(p_remove,list(decision_type=c("remove_low_cell_samples","override_parameters"),affected_samples="B",
  confirmed_by="user",user_statement="SIMULATED FIXTURE DECISION: remove B; minimum 75",removed_samples="B",
  parameter_overrides=list(doublet=list(min_cells=75))))
stopifnot(qc_precheck(area,p_remove$gse)=="PASS",qc_config(p_remove,qc_decision(p_remove,c("A","B")))$doublet$min_cells==75)
stopifnot(qc_run(area,p_remove$gse,step="qc1")=="COMPLETE_QC1")
removed <- readRDS(qc_stage_paths(p_remove,"qc1")$final)
removed_table <- utils::read.csv(file.path(qc_stage_paths(p_remove,"qc1")$qc,"cell_QC_all.csv"))
stopifnot(ncol(removed)==80L,all(removed$sample=="A"),nrow(removed_table)==160L,
  all(removed_table$removal_reason[removed_table$sample=="B"]=="user_confirmed_sample_removal"))
removal_counts <- counts; removal_counts[,1:30] <- 0; removal_counts[1:201,1:30] <- 1
removal <- Seurat::CreateSeuratObject(counts=removal_counts,min.cells=0,min.features=0,meta.data=meta)
p_removal <- scenario("GSE999999983",removal)
stopifnot(qc_precheck(area,p_removal$gse)=="NEEDS_USER_DECISION")
removal_report <- utils::read.csv(p_removal$report)
stopifnot(any(removal_report$sample=="ALL" & removal_report$metric=="predicted_removed_fraction" & removal_report$status=="NEEDS_USER_DECISION"))
write_decision(p_removal,list(decision_type="accept_qc_removal",affected_samples="A",confirmed_by="user",
  user_statement="SIMULATED FIXTURE DECISION: incomplete ALL acceptance",accept_qc_removal=TRUE))
stopifnot(qc_precheck(area,p_removal$gse)=="NEEDS_USER_DECISION")
write_decision(p_removal,list(decision_type="accept_qc_removal",affected_samples=c("A","B"),confirmed_by="user",
  user_statement="SIMULATED FIXTURE DECISION: accept current thresholds",accept_qc_removal=TRUE))
stopifnot(qc_precheck(area,p_removal$gse)=="PASS")
predicted <- utils::read.csv(p_removal$report)
stopifnot(as.integer(predicted$value[predicted$sample=="ALL" & predicted$metric=="predicted_after_qc_cells"][1])==130L)
p_species <- scenario("GSE999999984"); bad_species <- readRDS(p_species$raw); bad_species$species <- "mouse"; saveRDS(bad_species,p_species$raw)
stopifnot(qc_precheck(area,p_species$gse)=="NEEDS_USER_DECISION")
p_nuclei <- scenario("GSE999999985"); nuclei <- readRDS(p_nuclei$raw); nuclei$modality <- "snRNA-seq"; saveRDS(nuclei,p_nuclei$raw)
stopifnot(qc_precheck(area,p_nuclei$gse)=="NEEDS_USER_DECISION")
write_decision(p_nuclei,list(decision_type="accept_default_conflicts",affected_samples=c("A","B"),confirmed_by="user",
  user_statement="SIMULATED FIXTURE DECISION: retain defaults for nuclei",accept_default_conflicts=TRUE))
stopifnot(qc_precheck(area,p_nuclei$gse)=="PASS")
missing_df <- function(package) package!="DoubletFinder"
stopifnot(qc_precheck(area,gse,missing_df)=="FAILED_INPUT",qc_precheck(area,gse,missing_df,step="qc1")=="PASS")
coverage_counts <- counts; rownames(coverage_counts) <- paste0("ENSG",seq_len(nrow(coverage_counts)))
coverage <- Seurat::CreateSeuratObject(counts=coverage_counts,min.cells=0,min.features=0,meta.data=meta)
p_coverage <- scenario("GSE999999986",coverage)
stopifnot(qc_precheck(area,p_coverage$gse)=="NEEDS_USER_DECISION")
coverage_report <- utils::read.csv(p_coverage$report)
stopifnot(any(coverage_report$metric=="S_features_matched" & coverage_report$status=="NEEDS_USER_DECISION"),
  any(coverage_report$metric=="pMT_matched_features" & coverage_report$status=="NEEDS_USER_DECISION"))
ambiguous <- object; ambiguous$library_id <- c(NA,rep("x",159))
stopifnot(inherits(try(qc_capture_libraries(ambiguous),silent=TRUE),"try-error"))
# Existing checkpoints fail validation on changed inputs and remain byte-identical.
prior_md5 <- qc_md5(qc1_path); changed <- readRDS(raw); changed$group <- "Changed"; saveRDS(changed,raw)
stopifnot(inherits(try(qc_run(area,gse,step="qc1"),silent=TRUE),"try-error"),identical(qc_md5(qc1_path),prior_md5))
saveRDS(object,raw)
stopifnot(identical(qc_md5(raw),raw_md5))
# Missing optional DecontX must stop only QC-3: the earlier checkpoints stay intact.
stopifnot(file.exists(qc1_path),file.exists(qc2_path))
writeLines(area,file.path(root,".runtime","qc_stages_fixture_path.txt"))
cat("PASS: staged QC paths, retained-cell >3 gene rule, refreshed metrics, capture grouping, DF_adj, cycle scoring, separate DecontX assay, immutable inputs/checkpoints and figures.\n")
