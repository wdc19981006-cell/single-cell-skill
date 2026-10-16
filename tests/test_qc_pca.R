# Algorithm and plotting invariants independent of the historical GSE accession.
args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if(length(args)) args[1] else ".",winslash="/",mustWork=TRUE)
source(file.path(root,"qc/run_qc.R"))
cfg <- qc_config(list())
selection <- qc_pc_selection(c(20,15,14,10,9,8,7,6,5,4,1,1),cfg$doublet,12L)
pct <- c(20,15,14,10,9,8,7,6,5,4,1,1)
pct <- pct/sum(pct)*100
expected_drop <- max(which(head(pct,-1)-tail(pct,-1)>0.1))+1L
expected_variance <- which(cumsum(pct)>90 & pct<5)[1]
stopifnot(selection$drop_pick==expected_drop,expected_drop>2L,
          selection$pc_used==min(expected_drop,expected_variance),
          abs(sum(selection$variance_explained)-100)<1e-10,
          !isTRUE(all.equal(selection$variance_explained,selection$standard_deviation_share)),
          selection$selected_pc_cumulative_variance==sum(selection$variance_explained[seq_len(selection$pc_used)]))
only_drop <- cfg$doublet; only_drop$pc_cumulative_variance <- 100
only_cumulative <- cfg$doublet; only_cumulative$pc_elbow_drop <- 100
neither <- only_drop; neither$pc_elbow_drop <- 100
stopifnot(qc_pc_selection(c(8,5,3,1),only_drop,4)$drop_pick==4L,
          qc_pc_selection(rep(1,50),only_cumulative,50)$pc_used==46L,
          qc_pc_selection(rep(1,50),neither,50)$pc_used==20L,
          qc_pc_selection(rep(1,3),neither,3)$pc_used==3L,
          qc_pc_selection(c(2,1),neither,2)$pc_used==2L,
          inherits(try(qc_choose_pcs(c(2,1),neither,1),silent=TRUE),"try-error"),
          inherits(try(qc_choose_pcs(c(2,NA),neither,2),silent=TRUE),"try-error"),
          inherits(try(qc_choose_pcs(c(0,0),neither,2),silent=TRUE),"try-error"))
# Strict cumulative >90, not >=90; exact binary-representable uniform shares.
strict <- cfg$doublet; strict$pc_cumulative_variance <- 87.5; strict$pc_individual_variance <- 20
strict$pc_elbow_drop <- 100
stopifnot(qc_choose_pcs(rep(1,8),strict,8)==8L)
report <- qc_report_new()
low <- qc_pc_selection(c(100,rep(1,19)),cfg$doublet,20)
qc_report_pcs(report,"A",low)
qc_report_pcs(report,"ALL",selection,merged=TRUE)
records <- do.call(rbind,report$rows)
stopifnot(any(records$status=="PCA_LOW_DIM_WARNING"),
  all(c("max_pc","pc_used","variance_pick","drop_pick","selection_method",
    "selected_pc_cumulative_variance","merged_max_pc","merged_pc_used",
    "merged_variance_pick","merged_drop_pick") %in% records$metric))
# Real Seurat PCA exercises the dynamic bound and the small-object guard.
set.seed(700)
counts <- Matrix::Matrix(matrix(rpois(80*12,2),80,12,
  dimnames=list(paste0("g",1:80),paste0("c",1:12))),sparse=TRUE)
object <- Seurat::CreateSeuratObject(counts=counts)
object$pMT <- seq_len(12)
object$sample <- rep(c("A","B"),each=6)
small <- qc_preprocess_for_doublet(object,cfg,"small",cluster=FALSE)
stopifnot(small$max_pc==11L,small$pc_used>=2L,small$pc_used<=11L,
          length(small$object@graphs)==0L,!"umap" %in% names(small$object@reductions))
stopifnot(inherits(try(qc_preprocess_for_doublet(object[,1:3],cfg,"tiny",cluster=FALSE),silent=TRUE),"try-error"))
gene_limited <- qc_preprocess_for_doublet(object[1:6,],cfg,"few_genes",cluster=FALSE)
stopifnot(gene_limited$max_pc==min(50L,ncol(gene_limited$object)-1L,
  length(SeuratObject::VariableFeatures(gene_limited$object))-1L,nrow(gene_limited$object)-1L))
embedding <- matrix(seq_len(24),12,2,dimnames=list(colnames(object),c("UMAP_1","UMAP_2")))
object[["umap"]] <- SeuratObject::CreateDimReducObject(embeddings=embedding,key="UMAP_",assay="RNA")
object$DF <- rep(c("Doublet","Singlet","NotEvaluated"),4)
object$DF_adj <- object$DF
retained <- colnames(object)[object$DF_adj!="Doublet"]
panels <- sc_singlet_panels(object,retained)
stopifnot(identical(rownames(panels[[1]]$data),retained),
  identical(unname(panels[[1]]$data$x),unname(embedding[retained,1])),
  identical(unname(panels[[1]]$data$y),unname(embedding[retained,2])),
  "NotEvaluated" %in% panels[[2]]$data$DF_adj,
  all(vapply(panels,function(p) identical(p$coordinates$limits$x,sc_umap_limits(embedding)$x),logical(1))),
  identical(unname(sc_doublet_colors[c("Doublet","Singlet")]),c("#ea5455","#3f72af")))
object$pRP <- seq_len(12); object$pHB <- 0
colors <- sc_sample_colors(object$sample)
before <- sc_metrics_plot(object,colors)
after <- sc_metrics_plot(object[,object$sample=="B"],colors)
stopifnot(length(before)==5L,length(after)==5L,
  identical(vapply(seq_len(5),function(i) before[[i]]$labels$title,character(1)),
    c("pMT","pRP","pHB","nFeature_RNA","nCount_RNA")))
for(i in seq_len(5)) {
  ba <- ggplot2::ggplot_build(before[[i]]); aa <- ggplot2::ggplot_build(after[[i]])
  stopifnot(identical(ba$plot$scales$get_scales("fill")$map("B"),
                    aa$plot$scales$get_scales("fill")$map("B")))
}
cat("PASS: last elbow, strict cumulative, fallback, variance diagnostics, dynamic PCA, small sample, same-embedding subset, violin order and stable colors.\n")
