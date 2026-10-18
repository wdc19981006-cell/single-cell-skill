# Focused scientific safeguards, plus a reproducible synthetic QC object for CLI tests.
args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if(length(args)) args[1] else ".",winslash="/",mustWork=TRUE)
source(file.path(root,"analysis/run_analysis.R"))
cfg <- analysis_config(jsonlite::fromJSON(file.path(root,"analysis/default_config.json"),simplifyVector=FALSE))
meta <- data.frame(sample=rep(c("S1","S2","S3","S4"),each=10),
  group=rep(c("Normal","Tumor"),each=20),batch=rep(c("B1","B2"),20))
default <- analysis_design_check(meta,cfg)
stopifnot(!default$eligible,default$reason=="NO_CONFIRMED_TECHNICAL_BATCH")
cfg$harmony$technical_fields <- "batch"
crossed <- analysis_design_check(meta,cfg)
stopifnot(crossed$eligible,identical(crossed$selected,"batch"),!crossed$confounded)
confounded <- meta
confounded$batch <- confounded$group
stopifnot(analysis_design_check(confounded,cfg)$reason=="BATCH_GROUP_CONFOUNDED")
cfg$harmony$technical_fields <- character()
cfg$harmony$allow_sample <- TRUE
stopifnot(analysis_design_check(meta,cfg)$reason=="BATCH_GROUP_CONFOUNDED")
cfg$harmony$confounding_acknowledgement <- list(confirmed_by="user",user_statement="Synthetic test: compare with explicit confounding warning")
stopifnot(analysis_design_check(meta,cfg)$eligible,analysis_design_check(meta,cfg)$confounded)
bad <- cfg; bad$harmony$batch_field <- "group"
stopifnot(inherits(try(analysis_config(bad),silent=TRUE),"try-error"))
cfg$harmony$batch_field <- NULL
cfg$harmony$technical_fields <- "batch"
cfg$harmony$allow_sample <- FALSE
cfg$harmony$confounding_acknowledgement <- NULL
stopifnot(analysis_ari(c(1,1,2,2),c("b","b","a","a"))==1,
          abs(analysis_ari(c(1,1,2,2),c(1,2,1,2))+.5)<1e-10,
          analysis_ari(rep(1,4),rep("a",4))==1)
if(requireNamespace("clustree",quietly=TRUE)) {
  tree_data <- data.frame(cluster_res.0.2=c(0,0,0,1,1,1),cluster_res.0.6=c(0,0,1,2,2,2),check.names=FALSE)
  tree <- analysis_clustree_plot(tree_data)
  stopifnot(inherits(ggplot2::ggplotGrob(tree),"gtable"))
}
set.seed(813)
n <- 160L; genes <- 400L
biology <- rep(c(0,1),each=n/2)
batch <- rep(c(0,1),n/2)
means <- matrix(2,genes,n)
means[1:45,biology==1] <- 8
means[46:90,biology==0] <- 8
means[91:115,batch==1] <- 7
counts <- Matrix::Matrix(matrix(stats::rpois(genes*n,as.vector(means)),genes,n,
  dimnames=list(paste0("g",seq_len(genes)),paste0("cell",seq_len(n)))),sparse=TRUE)
object <- Seurat::CreateSeuratObject(counts=counts)
object$sample <- rep(c("S1","S2","S3","S4"),each=n/4)
object$group <- ifelse(biology==1,"Tumor","Normal")
object$batch <- ifelse(batch==1,"B2","B1")
object$database <- "GSE999999"
object$pMT <- runif(n,1,10)
object$S.Score <- rnorm(n)
object$G2M.Score <- rnorm(n)
object$Phase <- "G1"
original <- analysis_counts(object)
# A diagnostic QC PCA must disappear and be freshly recomputed.
object[["pca"]] <- SeuratObject::CreateDimReducObject(
  embeddings=matrix(0,n,3,dimnames=list(colnames(object),paste0("PC_",1:3))),key="PC_",assay="RNA")
object$seurat_clusters <- "stale_qc_cluster"
cfg$nfeatures <- 200L; cfg$max_pcs <- 12L
fresh <- analysis_pca(object,cfg)
stopifnot(identical(analysis_counts(fresh),original),!"seurat_clusters" %in% names(fresh@meta.data),
          any(SeuratObject::Embeddings(fresh[["pca"]])!=0),
          identical(fresh@misc$workflow_analysis$regress,"pMT"),
          fresh@misc$workflow_analysis$pc_selection$pc_used>=2L,
          fresh@misc$workflow_analysis$pc_selection$pc_used<=12L)
analysis_validate_object(fresh,"pca")
fractional <- object
fractional[["RNA"]] <- SeuratObject::CreateAssay5Object(counts=counts*0.5)
stopifnot(inherits(try(analysis_counts(fractional),silent=TRUE),"try-error"))
# An explicitly selected corrected assay may be fractional, original RNA persists.
object[["decontX"]] <- SeuratObject::CreateAssay5Object(counts=counts*.9)
corrected <- cfg; corrected$assay <- "decontX"
fresh_corrected <- analysis_pca(object,corrected)
stopifnot(identical(analysis_counts(fresh_corrected),original),
          identical(fresh_corrected@misc$workflow_analysis$assay,"decontX"),
          identical(SeuratObject::DefaultAssay(fresh_corrected),"decontX"))
if(length(args)>1L) {
  out <- args[2]
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  saveRDS(object,file.path(out,"synthetic_qc.rds"))
  jsonlite::write_json(list(nfeatures=200L,max_pcs=12L,
    harmony=list(technical_fields="batch"),markers=c("g1","g50","g100","absent_marker"),
    resolutions=c(.2,.6,1)),file.path(out,"analysis_config.json"),auto_unbox=TRUE,pretty=TRUE)
}
cat("PASS: design confounding/role/group guards, ARI, fresh formal PCA, counts preservation and corrected-assay selection\n")
