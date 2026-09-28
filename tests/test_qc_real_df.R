args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if(length(args)) args[1] else ".",winslash="/",mustWork=TRUE)
source(file.path(root,"qc","run_qc.R"))
stopifnot(isTRUE(getOption("r45.runtime.active")),
          Sys.getenv("R45_EXPECTED_LIBRARY") %in% .libPaths(),
          requireNamespace("DoubletFinder",quietly=TRUE))
set.seed(29)
features <- c(paste0("MT-",seq_len(4)),paste0("RPS",seq_len(4)),
              "HBA1","HBB",paste0("GENE",seq_len(300)))
cells <- paste0("SampleA_cell",seq_len(130))
counts <- Matrix::Matrix(matrix(stats::rpois(length(features)*length(cells),2),
  nrow=length(features),dimnames=list(features,cells)),sparse=TRUE)
object <- Seurat::CreateSeuratObject(counts=counts,min.cells=0,min.features=0,
  meta.data=data.frame(sample=rep("SampleA",length(cells)),
    database=rep("GSE999999996",length(cells)),group=rep("Test",length(cells)),row.names=cells))
object <- qc_add_percentages(object,"human")
report <- qc_report_new()
cfg <- qc_config(list())
result <- qc_run_doublet_sample(object,"SampleA",cfg,report,qc_default_doublet_api())
metrics <- do.call(rbind,report$rows)
stopifnot(ncol(result)==130L,all(c("pANN","DF","DF_adj") %in% names(result@meta.data)),
          !any(grepl("^seurat_clusters$|^RNA_snn_res\\.",names(result@meta.data))),
          all(result$DF %in% c("Singlet","Doublet")),
          all(result$DF_adj %in% c("Singlet","Doublet")),
          all(c("selected_pK","pc_used","nExp_poi","nExp_poi_adj") %in% metrics$metric))
cat("PASS: real DoubletFinder paramSweep, pK selection, two calls and classifications.\n")
