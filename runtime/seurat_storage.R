# Conservative, lossless compaction only at formal-save boundaries.
sc_object_inventory <- function(object) {
  rows <- list()
  add <- function(kind,name,x) rows[[length(rows)+1L]] <<- data.frame(
    kind=kind,name=name,bytes=as.numeric(object.size(x)),stringsAsFactors=FALSE)
  for(assay in names(object@assays)) {
    add("assay",assay,object[[assay]])
    for(layer in SeuratObject::Layers(object[[assay]]))
      add("layer",paste(assay,layer,sep="/"),SeuratObject::LayerData(object,assay=assay,layer=layer))
  }
  for(kind in c("reductions","graphs","neighbors"))
    for(name in names(slot(object,kind))) add(kind,name,slot(object,kind)[[name]])
  add("metadata","meta.data",object@meta.data)
  do.call(rbind,rows)
}

sc_expression_fingerprints <- function(object) {
  result <- list()
  for(assay in names(object@assays)) for(layer in intersect(c("counts","data"),SeuratObject::Layers(object[[assay]])))
    result[[paste(assay,layer,sep="/")]] <- digest::digest(
      SeuratObject::LayerData(object,assay=assay,layer=layer),algo="sha256")
  result
}

sc_compact_formal <- function(object,path,kind) {
  mode <- Sys.getenv("SINGLE_CELL_STORAGE_MODE")
  if(!mode %in% c("compact","full")) return(object)
  before <- sc_object_inventory(object)
  expression <- sc_expression_fingerprints(object)
  protected_meta <- object@meta.data
  removed <- character()
  if(kind=="analysis" && mode=="compact") {
    for(assay in names(object@assays)) {
      # Keep canonical counts/data exactly, including fractional decontX values.
      layers <- SeuratObject::Layers(object[[assay]])
      for(layer in layers) {
        if(layer %in% c("counts","data")) next
        duplicate <- any(vapply(intersect(c("counts","data"),layers),function(canonical)
          identical(SeuratObject::LayerData(object,assay=assay,layer=layer),
                    SeuratObject::LayerData(object,assay=assay,layer=canonical)),logical(1)))
        if(layer=="scale.data" || duplicate) {
          SeuratObject::LayerData(object[[assay]],layer=layer) <- NULL
          removed <- c(removed,paste(assay,layer,sep="/"))
        }
      }
    }
    # Preserve unique unintegrated UMAP, PCA and Harmony; remove numeric aliases.
    if("umap" %in% names(object@reductions)) for(name in setdiff(grep("umap",names(object@reductions),value=TRUE),"umap")) {
      a <- SeuratObject::Embeddings(object[[name]]); b <- SeuratObject::Embeddings(object[["umap"]])
      if(identical(rownames(a),rownames(b)) && identical(unname(a),unname(b))) {
        object[[name]] <- NULL; removed <- c(removed,name)
        if(identical(object@misc$workflow_analysis$cluster_review$umap,name))
          object@misc$workflow_analysis$cluster_review$umap <- "umap"
      }
    }
    candidates <- grep("^cluster_res\\.|^analysis_snn_res\\.",names(object@meta.data),value=TRUE)
    for(field in candidates) object[[field]] <- NULL
    protected_meta <- protected_meta[,setdiff(names(protected_meta),candidates),drop=FALSE]
    removed <- c(removed,candidates,names(object@graphs),names(object@neighbors))
    object@graphs <- list(); object@neighbors <- list()
  }
  if(!identical(expression,sc_expression_fingerprints(object)) || !identical(protected_meta,object@meta.data))
    stop("Compaction changed protected expression or metadata")
  object@misc$storage <- list(schema_version=1L,mode=mode,kind=kind,removed=removed,
    expression_sha256=expression,metadata_sha256=digest::digest(object@meta.data,algo="sha256"),
    clusters_sha256=if("seurat_clusters" %in% names(object@meta.data)) digest::digest(as.character(object$seurat_clusters),algo="sha256") else NULL)
  after <- sc_object_inventory(object)
  inventory <- merge(before,after,by=c("kind","name"),all=TRUE,suffixes=c("_before","_after"))
  utils::write.csv(inventory,file.path(dirname(path),"object_storage_inventory.csv"),row.names=FALSE,na="")
  object
}

sc_verify_formal <- function(object,operations=FALSE) {
  stored <- object@misc$storage
  if(is.null(stored)) stop("Formal storage signature missing")
  if(!identical(stored$expression_sha256,sc_expression_fingerprints(object)) ||
     !identical(stored$metadata_sha256,digest::digest(object@meta.data,algo="sha256")))
    stop("Formal expression or metadata readback differs")
  methods::validObject(object)
  if(operations) {
    assay <- object@misc$workflow_analysis$assay
    if(!"data" %in% SeuratObject::Layers(object[[assay]])) stop("Normalized data missing")
    SeuratObject::DefaultAssay(object) <- assay
    gene <- rownames(object[[assay]])[1]
    plot <- Seurat::FeaturePlot(object,features=gene,reduction="umap")
    ggplot2::ggplot_build(plot)
    # Bounded API smoke test, using actual formal labels and unchanged data.
    ids <- as.character(object$seurat_clusters)
    groups <- names(which(table(ids)>=3L))
    technical_contrast <- length(groups)<2L
    if(technical_contrast && ncol(object)<6L) stop("Too few cells for the bounded Marker API smoke test")
    cells <- if(technical_contrast) head(colnames(object),60L) else
      unlist(lapply(groups[1:2],function(g) head(colnames(object)[ids==g],30L)),use.names=FALSE)
    small <- subset(object,cells=cells)
    if(technical_contrast) {
      # Only an in-memory API test; no clusters/results are changed or reported.
      groups <- c("api_a","api_b")
      SeuratObject::Idents(small) <- rep(groups,length.out=ncol(small))
    } else SeuratObject::Idents(small) <- small$seurat_clusters
    suppressWarnings(Seurat::FindMarkers(small,ident.1=groups[1],ident.2=groups[2],
      features=head(rownames(small[[assay]]),20L),min.pct=0,logfc.threshold=0,test.use="wilcox",verbose=FALSE))
    small$storage_annotation_test <- as.character(small$seurat_clusters)
    stopifnot(identical(unname(as.character(small$storage_annotation_test)),unname(as.character(small$seurat_clusters))))
  }
  invisible(TRUE)
}
