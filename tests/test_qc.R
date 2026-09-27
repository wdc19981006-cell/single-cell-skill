args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if(length(args)) args[1] else ".",winslash="/",mustWork=TRUE)
source(file.path(root,"qc","run_qc.R"))
stopifnot(grepl("4.5.0",R.version.string,fixed=TRUE),
          "D:/R/R-4.5.0/library" %in% .libPaths(),
          requireNamespace("DoubletFinder",quietly=TRUE))
stopifnot(identical(unname(sc_plot_dimensions("umap_2panel")),c(11,5)),
          identical(unname(sc_plot_dimensions("umap_3panel")),c(16,5)),
          get_umap_point_size(9999)==0.5,get_umap_point_size(10000)==0.3,
          get_umap_point_size(50001)==0.15,get_umap_point_size(150001)==0.08)
mock_embedding <- matrix(c(0,1,2,4,4,8),ncol=2)
limits <- sc_umap_limits(mock_embedding)
mock_frame <- data.frame(x=mock_embedding[,1],y=mock_embedding[,2],sample=c("A","B","C"),
                         DF=c("Doublet","Singlet","NotEvaluated"),DF_adj=c("Doublet","Singlet","NotEvaluated"))
umap_panels <- lapply(c("sample","DF","DF_adj"),function(field)
  sc_umap_panel(mock_frame,field,limits,3L))
stopifnot(all(vapply(umap_panels,function(p) identical(p$coordinates$ratio,1) &&
                       !p$coordinates$is_free(),logical(1))),
          all(vapply(umap_panels,function(p) identical(p$coordinates$limits$x,limits$x) &&
                       identical(p$coordinates$limits$y,limits$y),logical(1))),
          tail(umap_panels[[2]]$data$DF,1)=="Doublet")
set.seed(17)
source_genes <- Seurat::cc.genes
s_genes <- unique(source_genes$s.genes)[seq_len(24)]
g2m_genes <- setdiff(unique(source_genes$g2m.genes),s_genes)[seq_len(24)]
features <- unique(c(paste0("MT-",seq_len(4)),paste0("RPS",seq_len(4)),"HBA1","HBB",
                     s_genes,g2m_genes,paste0("GENE",seq_len(950))))

make_fixture <- function(base,gse,low=TRUE) {
  dataset <- file.path(base,"data",gse)
  dir.create(dataset,recursive=TRUE)
  sample <- c(rep("SampleA",125),if(low) rep("SampleB",42) else character())
  cells <- paste0(sample,"_cell",seq_along(sample))
  dense <- matrix(stats::rpois(length(features)*length(cells),lambda=2),nrow=length(features),
                  dimnames=list(features,cells))
  counts <- Matrix::Matrix(dense,sparse=TRUE)
  meta <- data.frame(sample=sample,database=gse,group=ifelse(sample=="SampleA","Tumor","Control"),
                     row.names=cells)
  object <- Seurat::CreateSeuratObject(counts=counts,min.cells=3,min.features=200,meta.data=meta)
  saveRDS(object,file.path(dataset,"seurat_raw.rds"))
  invisible(object)
}

mock_api <- list(
  paramSweep=function(seu,PCs,sct) list(PCs=PCs),
  summarizeSweep=function(sweep,GT) sweep,
  find.pK=function(stats) data.frame(pK=c("0.05","0.1"),BCmetric=c(1,2)),
  modelHomotypic=function(clusters) 0.2,
  doubletFinder=function(seu,PCs,pN,pK,nExp,reuse.pANN,sct) {
    score <- seq_len(ncol(seu))/ncol(seu)
    if(is.null(reuse.pANN)) seu[[paste0("pANN_",nExp)]] <- score
    calls <- rep("Singlet",ncol(seu))
    calls[tail(seq_len(ncol(seu)),nExp)] <- "Doublet"
    seu[[paste0("DF.classifications_",nExp)]] <- calls
    seu
  }
)
bad_pk_api <- mock_api
bad_pk_api$find.pK <- function(stats) data.frame(pK=NA_real_,BCmetric=NA_real_)
pk_error <- tryCatch(qc_select_pk(NULL,3L,bad_pk_api,"SampleA"),error=conditionMessage)
stopifnot(grepl("pK search stopped for sample SampleA",pk_error,fixed=TRUE),
          qc_choose_pcs(c(10,9,8),list(pc_cumulative_variance=90,
            pc_individual_variance=5,pc_elbow_drop=100),3L)==3L)

temporary <- tempfile("qc-fixture-",tmpdir=file.path(root,"data"))
dir.create(temporary)
missing_gse <- "GSE999999995"
dir.create(file.path(temporary,"data",missing_gse),recursive=TRUE)
stopifnot(identical(qc_precheck(temporary,missing_gse),"FAILED_INPUT"),
          !dir.exists(file.path(temporary,"data",missing_gse,"raw")),
          grepl("seurat_raw.rds missing",paste(utils::read.csv(
            file.path(temporary,"data",missing_gse,"qc","qc_report.csv"))$note,collapse=";"),fixed=TRUE))
unicode_report <- qc_report_new()
qc_report_add(unicode_report,"RUN",metric="user_statement",value="用户确认保留样本")
unicode_path <- file.path(temporary,"unicode_report.csv")
qc_report_write(unicode_report,unicode_path)
stopifnot(identical(utils::read.csv(unicode_path)$value,
                    "用户确认保留样本"))
gse <- "GSE999999997"
original <- make_fixture(temporary,gse)
paths <- qc_paths(temporary,gse)
raw_checksum <- unname(tools::md5sum(paths$raw))
available <- function(p) TRUE
status <- qc_precheck(temporary,gse,available)
stopifnot(identical(status,"NEEDS_USER_DECISION"),
          any(utils::read.csv(paths$report)$status=="LOW_CELL_SAMPLE"),
          !file.exists(paths$final))

jsonlite::write_json(list(decision_type="remove_low_cell_samples",affected_samples="SampleB",
  user_statement="SIMULATED USER DECISION: remove SampleB",confirmed_by="user",
  removed_samples="SampleB",kept_low_cell_samples=character()),paths$decision,auto_unbox=TRUE)
stopifnot(identical(qc_precheck(temporary,gse,available),"PASS"))
failed_run <- tryCatch(qc_run(temporary,gse,available,bad_pk_api),error=conditionMessage)
failed_report <- utils::read.csv(paths$report,stringsAsFactors=FALSE)
stopifnot(grepl("pK search stopped",failed_run,fixed=TRUE),!file.exists(paths$final),
          sum(failed_report$section=="RUN" & failed_report$metric=="status")==1L,
          failed_report$value[failed_report$section=="RUN" & failed_report$metric=="status"]=="FAILED")
status <- qc_run(temporary,gse,available,mock_api)
stopifnot(identical(status,"COMPLETE_QC"))
final <- readRDS(paths$final)
report <- utils::read.csv(paths$report,stringsAsFactors=FALSE)
stopifnot(identical(raw_checksum,unname(tools::md5sum(paths$raw))),
          !dir.exists(file.path(temporary,"data",gse,"raw")),
          all(final$sample=="SampleA"),all(final$DF_adj=="Singlet"),
          all(c("pANN","DF","DF_adj","doublet_status","S.Score","G2M.Score","Phase","CC.Difference") %in% names(final@meta.data)),
          identical(SeuratObject::Layers(final[["RNA"]]),"counts"),
          length(final@reductions)==0L,length(final@graphs)==0L,length(final@neighbors)==0L,
          !any(grepl("^seurat_clusters$|^RNA_snn_res\\.",names(final@meta.data))),
          identical(qc_counts(final),qc_counts(original)[,colnames(final),drop=FALSE]),
          all(c("CONFIG","SAMPLE","CELL_QC","DOUBLET","CELL_CYCLE","PLOT","RUN") %in% report$section),
          any(report$status=="USER_CONFIRMED_REMOVAL"))
stopifnot(sum(report$section=="RUN" & report$metric=="status")==1L,
          report$value[report$section=="RUN" & report$metric=="status"]=="COMPLETE_QC")
expected_plots <- as.vector(outer(c("before_QC","after_QC","doublet_umap","doublet_vlnplot",
  "cell_cycle_phase","cell_cycle_score"),c(".pdf",".png"),paste0))
plot_files <- file.path(paths$qc,expected_plots)
png_header <- file(file.path(paths$qc,"before_QC.png"),open="rb")
invisible(seek(png_header,where=16L))
png_dimensions <- readBin(png_header,integer(),n=2L,size=4L,endian="big")
close(png_header)
stopifnot(all(file.exists(plot_files)),all(file.info(plot_files)$size>0),
          identical(png_dimensions,c(3000L,1800L)),
          setequal(report$metric[report$section=="PLOT"],expected_plots),
          any(report$metric=="doublet_umap_before_filter"),
          length(list.files(paths$qc,pattern="^qc_report\\.csv$"))==1L)
stopifnot(!length(list.files(paths$qc,pattern="summary\\.csv$")))

gse <- "GSE999999998"
original <- make_fixture(temporary,gse)
paths <- qc_paths(temporary,gse)
writeLines(paste0(
  '{"decision_type":"keep_without_doubletfinder","affected_samples":"SampleB",',
  '"user_statement":"模拟用户确认：保留SampleB但不运行DoubletFinder",',
  '"confirmed_by":"user","removed_samples":[],"kept_low_cell_samples":"SampleB"}'),
  paths$decision,useBytes=TRUE)
stopifnot(identical(qc_precheck(temporary,gse,available),"PASS"))
status <- qc_run(temporary,gse,available,mock_api)
final <- readRDS(paths$final)
report <- utils::read.csv(paths$report,stringsAsFactors=FALSE)
stopifnot(identical(status,"QC_COMPLETE_WITH_UNEVALUATED_DOUBLETS"),
          all(final$DF_adj[final$sample=="SampleB"]=="NotEvaluated"),
          all(final$DF[final$sample=="SampleB"]=="NotEvaluated"),
          all(final$doublet_status[final$sample=="SampleB"]=="SKIPPED_LOW_CELL"),
          all(final$DF_adj %in% c("Singlet","NotEvaluated")),
          any(report$value=="QC_COMPLETE_WITH_UNEVALUATED_DOUBLETS"),
          sum(report$section=="RUN" & report$metric=="status")==1L,
          any(report$metric=="unevaluated_cells" & report$sample=="SampleB"),
          any(grepl("模拟用户确认",report$note,fixed=TRUE)))
gse <- "GSE999999999"
make_fixture(temporary,gse)
paths <- qc_paths(temporary,gse)
jsonlite::write_json(list(decision_type="override_parameters",affected_samples="SampleB",
  user_statement="SIMULATED USER DECISION: use minimum 40 cells",confirmed_by="user",
  parameter_overrides=list(doublet=list(min_cells=40))),paths$decision,auto_unbox=TRUE)
stopifnot(identical(qc_precheck(temporary,gse,available),"PASS"))
stopifnot(identical(qc_run(temporary,gse,available,mock_api),"COMPLETE_QC"))
final <- readRDS(paths$final)
stopifnot(all(final$DF_adj=="Singlet"),all(final$doublet_status=="EVALUATED"),
          setequal(unique(final$sample),c("SampleA","SampleB")))
missing <- function(p) p!="DoubletFinder"
stopifnot(identical(qc_precheck(temporary,gse,missing),"FAILED_INPUT"))
report <- utils::read.csv(paths$report,stringsAsFactors=FALSE)
stopifnot(any(report$metric=="missing_package" & report$value=="DoubletFinder" &
  grepl("R version 4.5.0",report$note,fixed=TRUE) &
  grepl("D:/R/R-4.5.0/library",report$note,fixed=TRUE)))
unlink(temporary,recursive=TRUE)
cat("PASS: QC precheck, all three user decisions, missing-package report, per-sample doublet control flow, raw counts, reports, PDFs and clean final object.\n")
