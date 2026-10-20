# Real Seurat 5 layer/reduction/metadata compaction; no biological parameter changes.
args <- commandArgs(trailingOnly=TRUE)
source(file.path(args[1],"analysis/run_analysis.R"))
set.seed(91)
x <- Matrix::Matrix(matrix(rpois(60*40,3),60,40,
  dimnames=list(paste0("gene",1:60),paste0("cell",1:40))),sparse=TRUE)
object <- Seurat::CreateSeuratObject(x)
object$sample <- rep(c("A","B"),each=20); object$pMT <- seq_len(40)/10
object$group <- rep(c("Control","Tumor"),each=20)
object$S.Score <- 0; object$G2M.Score <- 0; object$Phase <- "G1"
object$decontX_contamination <- .125
object[["decontX"]] <- SeuratObject::CreateAssay5Object(counts=x*.875)
SeuratObject::DefaultAssay(object) <- "decontX"
object <- Seurat::NormalizeData(object,verbose=FALSE)
object <- Seurat::FindVariableFeatures(object,nfeatures=30,verbose=FALSE)
object <- Seurat::ScaleData(object,verbose=FALSE)
SeuratObject::LayerData(object[["decontX"]],layer="data.copy") <- SeuratObject::LayerData(object,assay="decontX",layer="data")
object <- Seurat::RunPCA(object,npcs=5,verbose=FALSE)
object[["harmony"]] <- object[["pca"]]
embedding <- matrix(rnorm(80),40,2,dimnames=list(colnames(object),c("UMAP_1","UMAP_2")))
object[["umap"]] <- SeuratObject::CreateDimReducObject(embedding,key="UMAP_",assay="decontX")
object[["umap_harmony"]] <- object[["umap"]]
object[["umap_unintegrated"]] <- SeuratObject::CreateDimReducObject(embedding+1,key="RAW_",assay="decontX")
object$cluster_res.0.3 <- rep(c("0","1"),each=20)
object$cluster_res.0.6 <- rep(c("0","1","2","3"),each=10)
object$seurat_clusters <- factor(object$cluster_res.0.3)
SeuratObject::Idents(object) <- object$seurat_clusters
object@misc$workflow_analysis <- list(assay="decontX",formal_pca=TRUE,regress="pMT",
  pc_selection=list(pc_used=3L),harmony_completed=TRUE,
  cluster_review=list(resolutions=c(.3,.6),columns=c("cluster_res.0.3","cluster_res.0.6"),umap="umap_harmony",selected=TRUE),
  resolution_decision=list(confirmed_by="user",user_statement="SYNTHETIC TEST ONLY",resolution=.3))
out <- args[2]; dir.create(out,recursive=TRUE,showWarnings=FALSE)
original <- sc_expression_fingerprints(object)
Sys.setenv(SINGLE_CELL_STORAGE_MODE="compact")
compact <- sc_compact_formal(object,file.path(out,"compact.rds"),"analysis")
stopifnot(identical(original,sc_expression_fingerprints(compact)),
  !"scale.data" %in% SeuratObject::Layers(compact[["decontX"]]),
  !"data.copy" %in% SeuratObject::Layers(compact[["decontX"]]),
  !"umap_harmony" %in% names(compact@reductions),
  all(c("umap","umap_unintegrated","pca","harmony") %in% names(compact@reductions)),
  !any(grepl("^cluster_res",names(compact@meta.data))),
  identical(compact@misc$workflow_analysis$cluster_review$umap,"umap"))
analysis_validate_object(compact,"finalize")
saveRDS(compact,file.path(out,"compact.rds")); sc_verify_formal(readRDS(file.path(out,"compact.rds")),operations=TRUE)
bad <- compact; bad[["decontX"]]@layers$counts@x[1] <- bad[["decontX"]]@layers$counts@x[1]+.001
stopifnot(inherits(try(sc_verify_formal(bad),silent=TRUE),"try-error"))
Sys.setenv(SINGLE_CELL_STORAGE_MODE="full")
full <- sc_compact_formal(object,file.path(out,"full.rds"),"analysis")
stopifnot("scale.data" %in% SeuratObject::Layers(full[["decontX"]]),
  all(c("cluster_res.0.3","cluster_res.0.6") %in% names(full@meta.data)),
  identical(original,sc_expression_fingerprints(full)))
analysis_validate_object(full,"finalize")
qc <- sc_compact_formal(object,file.path(out,"qc.rds"),"qc")
saveRDS(qc,file.path(out,"qc.rds")); sc_verify_formal(readRDS(file.path(out,"qc.rds")))
route_root <- file.path(out,"routing"); dir.create(route_root,showWarnings=FALSE)
dataset <- file.path(route_root,"results","StorageRouting")
work <- file.path(dataset,".work","checkpoints","qc","rebuilds","test")
reports <- file.path(dataset,"02_QC","rebuilds","test")
Sys.setenv(SINGLE_CELL_QC_CHECKPOINT_ROOT=work,SINGLE_CELL_QC_RUN_ROOT=reports)
routes <- qc_paths(route_root,"StorageRouting")
stopifnot(identical(routes$qc,reports),
  identical(qc_stage_paths(routes,"qc1")$final,file.path(work,"seurat_QC-1.rds")),
  identical(qc_stage_paths(routes,"qc3")$final,file.path(reports,"03_ambient_RNA","seurat_QC-3.rds")),
  identical(routes$decision,file.path(dataset,"02_QC","qc_decision.json")))
Sys.unsetenv(c("SINGLE_CELL_QC_CHECKPOINT_ROOT","SINGLE_CELL_QC_RUN_ROOT"))
cat("PASS: Seurat 5 lossless compact/full, fractional corrected counts, independent fingerprints, FeaturePlot, FindMarkers, annotation API\n")
