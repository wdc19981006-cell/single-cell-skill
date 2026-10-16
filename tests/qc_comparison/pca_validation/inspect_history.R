# Read-only inventory and raw-count checks of historical inputs.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==3L,isTRUE(getOption("r45.runtime.active")))
source(file.path(args[1],"qc/run_qc.R"))
dataset <- args[2]; output <- args[3]
dir.create(output,recursive=TRUE,showWarnings=FALSE)
raw <- readRDS(file.path(dataset,"seurat_raw.rds"))
rc <- qc_validate_input(raw,basename(dataset))
cfg <- qc_config(list())
raw <- qc_add_percentages(raw,cfg$species)
cell_qc_ids <- colnames(raw)[qc_cell_keep(raw@meta.data,cfg$cell_qc)]
write.csv(data.frame(cell=cell_qc_ids),file.path(output,"cell_qc_ids.csv"),row.names=FALSE)
summaries <- list()
for(version in c("raw","original_pre_doublet","original_final","prior_fixed_skill")) {
  path <- switch(version,raw="seurat_raw.rds",original_pre_doublet="seurat.qc.rds",
    original_final="seurat.all_qc.rds",prior_fixed_skill="qc/seurat_qc.rds")
  object <- if(version=="raw") raw else readRDS(file.path(dataset,path))
  counts <- qc_counts(object)
  common_cells <- intersect(colnames(raw),colnames(object))
  common_genes <- intersect(rownames(raw),rownames(object))
  mismatch <- 0L
  if(length(common_cells) && length(common_genes)) for(start in seq.int(1L,length(common_cells),by=256L)) {
    block <- common_cells[start:min(start+255L,length(common_cells))]
    d <- counts[common_genes,block,drop=FALSE]-rc[common_genes,block,drop=FALSE]
    mismatch <- mismatch+sum(d@x!=0)
  }
  summaries[[version]] <- list(path=path,cells=ncol(object),genes=nrow(object),
    RNA_layers=SeuratObject::Layers(object[["RNA"]]),reductions=names(object@reductions),
    metadata=names(object@meta.data),common_raw_cells=length(common_cells),
    common_raw_genes=length(common_genes),counts_mismatches=mismatch,
    cells_match_current_cell_qc=setequal(colnames(object),cell_qc_ids),
    DF=if("DF" %in% names(object@meta.data)) as.list(table(object$DF)) else NULL,
    DF_adj=if("DF_adj" %in% names(object@meta.data)) as.list(table(object$DF_adj)) else NULL)
  write.csv(data.frame(cell=colnames(object),object@meta.data),
    file.path(output,paste0(version,"_metadata.csv")),row.names=FALSE)
  if("umap" %in% names(object@reductions)) {
    e <- SeuratObject::Embeddings(object[["umap"]])
    write.csv(data.frame(cell=rownames(e),x=e[,1],y=e[,2]),file.path(output,paste0(version,"_umap_cells.csv")),row.names=FALSE)
  }
  cat("HISTORY",version,"cells",ncol(object),"genes",nrow(object),"raw_mismatch",mismatch,"\n")
  if(version!="raw") {rm(object,counts);gc()}
}
jsonlite::write_json(summaries,file.path(output,"history_inventory.json"),pretty=TRUE,auto_unbox=TRUE)
