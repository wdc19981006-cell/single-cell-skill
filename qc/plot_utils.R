sc_plot_dimensions <- function(type) {
  sizes <- list(umap_single=c(6,5), umap_2panel=c(11,5), umap_3panel=c(16,5),
                violin=c(10,6), qc_violin=c(20,8), large_sample_violin=c(20,8), feature_3panel=c(15,5))
  size <- sizes[[type]]
  if (is.null(size)) stop("Unknown single-cell plot type: ",type)
  setNames(size,c("width","height"))
}

get_umap_point_size <- function(n_cells) {
  if (length(n_cells)!=1L || !is.finite(n_cells) || n_cells<0) stop("Invalid UMAP cell count")
  if (n_cells<10000) 0.5 else if (n_cells<=50000) 0.3 else if (n_cells<=150000) 0.15 else 0.08
}

theme_sc <- function() {
  ggplot2::theme_bw(base_size=10) + ggplot2::theme(
    plot.title=ggplot2::element_text(size=13),
    axis.title=ggplot2::element_text(size=10),
    legend.title=ggplot2::element_text(size=10),
    legend.text=ggplot2::element_text(size=9))
}

theme_sc_umap <- function() {
  theme_sc() + ggplot2::theme(panel.grid=ggplot2::element_blank())
}

sc_umap_limits <- function(embedding) {
  if (!is.matrix(embedding) || ncol(embedding)<2L || !all(is.finite(embedding[,1:2])))
    stop("UMAP embedding must have two finite dimensions")
  padded <- function(x) {
    extent <- range(x)
    margin <- max(diff(extent)*0.03, .Machine$double.eps^0.5)
    extent+c(-margin,margin)
  }
  list(x=padded(embedding[,1]),y=padded(embedding[,2]))
}

sc_sample_legend <- function(n_samples) {
  if (n_samples<=10L) return(ggplot2::guides(color=ggplot2::guide_legend(ncol=1,override.aes=list(size=3))))
  list(ggplot2::theme(legend.position="bottom"),
       ggplot2::guides(color=ggplot2::guide_legend(ncol=min(10L,ceiling(n_samples/4)),
                                                 byrow=TRUE,override.aes=list(size=3))))
}

sc_sample_colors <- function(samples) {
  samples <- sort(unique(as.character(samples)))
  stats::setNames(grDevices::hcl.colors(length(samples),palette="Dynamic"),samples)
}
sc_doublet_colors <- c(Doublet="#ea5455",Singlet="#3f72af",NotEvaluated="#858585")

sc_umap_panel <- function(frame,field,limits,n_samples,sample_colors=sc_sample_colors(frame$sample)) {
  if (!field %in% c("sample","DF","DF_adj")) stop("Unsupported UMAP panel field")
  points <- frame
  if (field!="sample") points <- points[order(points[[field]]=="Doublet"),,drop=FALSE]
  panel <- ggplot2::ggplot(points,ggplot2::aes(x=x,y=y,color=.data[[field]])) +
    ggplot2::geom_point(size=0.2) +
    ggplot2::labs(title=field,color=field,x="UMAP_1",y="UMAP_2") +
    ggplot2::coord_fixed(xlim=limits$x,ylim=limits$y,expand=FALSE) + theme_sc_umap() +
    ggplot2::theme(panel.border=ggplot2::element_blank(),
      plot.title=ggplot2::element_text(hjust=0.5,face="bold")) +
    ggplot2::guides(color=ggplot2::guide_legend(override.aes=list(size=3)))
  if (field=="sample") panel + ggplot2::scale_color_manual(values=sample_colors) + sc_sample_legend(n_samples)
  else panel + ggplot2::scale_color_manual(
    values=sc_doublet_colors,drop=FALSE)
}

save_sc_plot <- function(plot,name,type,output_dir,n_cells=NULL,n_samples=NULL,report=NULL) {
  if (!grepl("^[A-Za-z0-9_]+$",name)) stop("Plot name must be a filename stem")
  if (!dir.exists(output_dir)) stop("Plot output directory is missing")
  size <- sc_plot_dimensions(type)
  if(type %in% c("qc_violin","large_sample_violin") && !is.null(n_samples))
    size[["width"]] <- max(20,0.6*n_samples)
  if(type %in% c("umap_2panel","umap_3panel") && !is.null(n_samples) && n_samples>10L)
    size[["height"]] <- size[["height"]]+1+0.22*ceiling(n_samples/min(10L,ceiling(n_samples/4)))
  if(type=="umap_2panel" && !is.null(n_samples) && n_samples>25L)
    size[["width"]] <- 16
  paths <- file.path(output_dir,paste0(name,c(".pdf",".png")))
  for (path in paths) {
    ggplot2::ggsave(path,plot=plot,width=size[["width"]],height=size[["height"]],
                    units="in",dpi=300,limitsize=FALSE,bg="white")
    if (!file.exists(path) || file.info(path)$size<=0) stop("Plot file missing or empty: ",path)
    if (!is.null(report)) qc_report_add(report,"PLOT",metric=basename(path),value=path,
      note=paste0("type=",type,"; cells=",n_cells %||% "NA","; samples=",n_samples %||% "NA"))
  }
  invisible(paths)
}

sc_metrics_plot <- function(object,sample_colors) {
  metrics <- c("pMT","pRP","pHB","nFeature_RNA","nCount_RNA")
  object$sample <- factor(as.character(object$sample),levels=names(sample_colors))
  Seurat::VlnPlot(object,group.by="sample",features=metrics,pt.size=0,ncol=3,
                  cols=sample_colors) & Seurat::NoLegend() &
    ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5,size=16,face="bold"),
      axis.text.x=ggplot2::element_text(angle=50,size=8,vjust=1,hjust=1))
}

qc_plot_metrics <- function(object,paths,name,title,report) {
  # The pre-filter sample universe fixes colors even if an entire sample disappears.
  if(is.null(report$sample_colors)) report$sample_colors <- sc_sample_colors(object$sample)
  plot <- sc_metrics_plot(object,report$sample_colors)
  n_samples <- length(report$sample_colors)
  save_sc_plot(plot,name,"qc_violin",
               paths$qc,ncol(object),n_samples,report)
}

qc_plot_doublets <- function(object,paths,report) {
  embedding <- SeuratObject::Embeddings(object[["umap"]])
  metadata <- object@meta.data[rownames(embedding),,drop=FALSE]
  frame <- data.frame(x=embedding[,1],y=embedding[,2],sample=as.character(metadata$sample),
                      DF=as.character(metadata$DF),DF_adj=as.character(metadata$DF_adj))
  limits <- sc_umap_limits(embedding)
  n_samples <- length(unique(frame$sample))
  sample_colors <- report$sample_colors %||% sc_sample_colors(frame$sample)
  panels <- lapply(c("sample","DF","DF_adj"),function(field)
    sc_umap_panel(frame,field,limits,n_samples,sample_colors))
  combined <- patchwork::wrap_plots(panels,ncol=3,guides=if(n_samples>10L) "collect" else "keep")
  if(n_samples>10L) combined <- combined & ggplot2::theme(legend.position="bottom")
  if(n_samples>25L) combined <- combined & ggplot2::theme(legend.box="vertical")
  save_sc_plot(combined,"doublet_umap","umap_3panel",
               paths$qc,nrow(frame),n_samples,report)
  violin <- rbind(data.frame(class=metadata$DF,call="DF",nFeature_RNA=metadata$nFeature_RNA),
                  data.frame(class=metadata$DF_adj,call="DF_adj",nFeature_RNA=metadata$nFeature_RNA))
  violin$call <- factor(violin$call,levels=c("DF","DF_adj"))
  violin$class <- factor(violin$class,levels=names(sc_doublet_colors))
  plot <- ggplot2::ggplot(violin,ggplot2::aes(x=class,y=nFeature_RNA,fill=class)) +
    ggplot2::geom_violin(scale="width") + ggplot2::geom_boxplot(width=0.15,outlier.shape=NA,fill="white") +
    ggplot2::scale_fill_manual(values=sc_doublet_colors) +
    ggplot2::facet_wrap(~call,ncol=2) + theme_sc() +
    ggplot2::labs(x="DoubletFinder call",y="nFeature_RNA")
  save_sc_plot(plot,"doublet_vlnplot","violin",paths$qc,nrow(frame),n_samples,report)
}

sc_singlet_panels <- function(object,retained,sample_colors=sc_sample_colors(object$sample)) {
  embedding <- SeuratObject::Embeddings(object[["umap"]])
  stopifnot(length(retained)>0L,!anyDuplicated(retained),all(retained %in% rownames(embedding)))
  metadata <- object@meta.data[retained,,drop=FALSE]
  if(any(!metadata$DF_adj %in% c("Singlet","NotEvaluated"))) stop("Singlet plot includes rejected doublets")
  frame <- data.frame(x=embedding[retained,1],y=embedding[retained,2],
    sample=as.character(metadata$sample),DF_adj=as.character(metadata$DF_adj),row.names=retained)
  # Use the full pre-filter limits and the exact subset, never refit UMAP.
  limits <- sc_umap_limits(embedding)
  lapply(c("sample","DF_adj"),function(field)
    sc_umap_panel(frame,field,limits,length(sample_colors),sample_colors))
}

qc_plot_singlets <- function(object,retained,paths,report) {
  colors <- report$sample_colors %||% sc_sample_colors(object$sample)
  panels <- sc_singlet_panels(object,retained,colors)
  combined <- patchwork::wrap_plots(panels,ncol=2,guides=if(length(colors)>10L) "collect" else "keep")
  if(length(colors)>10L) combined <- combined & ggplot2::theme(legend.position="bottom")
  if(length(colors)>25L) combined <- combined & ggplot2::theme(legend.box="vertical")
  save_sc_plot(combined,"singlet_dimplot","umap_2panel",paths$qc,length(retained),length(colors),report)
  qc_report_add(report,"RUN",metric="singlet_umap_embedding",value="EXACT_PRE_FILTER_SUBSET",status="PASS",
    note="Same coordinates and limits as doublet_umap; approved NotEvaluated cells keep their status")
}

qc_plot_cycle <- function(object,paths,report) {
  embedding <- SeuratObject::Embeddings(object[["umap"]])
  limits <- sc_umap_limits(embedding)
  n_samples <- length(unique(object$sample))
  phase <- sc_cycle_phase_panel(object,limits)
  save_sc_plot(phase,"cell_cycle_phase","umap_single",paths$qc,ncol(object),n_samples,report)
  metrics <- c("S.Score","G2M.Score","CC.Difference")
  panels <- lapply(metrics,function(metric) sc_cycle_feature_panel(object,metric,limits))
  save_sc_plot(patchwork::wrap_plots(panels,ncol=3),"cell_cycle_score","feature_3panel",
               paths$qc,ncol(object),n_samples,report)
}

sc_cycle_phase_panel <- function(object,limits) {
  object$Phase <- factor(as.character(object$Phase),levels=c("G1","S","G2M"))
  Seurat::DimPlot(object,reduction="umap",group.by="Phase",pt.size=get_umap_point_size(ncol(object))) +
    ggplot2::scale_color_manual(values=c(G1="#4C78A8",S="#F2A541",G2M="#59A14F"),drop=FALSE) +
    ggplot2::labs(title="Cell-cycle phase",x="UMAP_1",y="UMAP_2") +
    ggplot2::coord_fixed(xlim=limits$x,ylim=limits$y,expand=FALSE) + theme_sc_umap()
}

sc_cycle_feature_panel <- function(object,metric,limits) {
  Seurat::FeaturePlot(object,features=metric,reduction="umap",combine=FALSE,
                      pt.size=get_umap_point_size(ncol(object)))[[1]] +
    ggplot2::labs(title=metric,x="UMAP_1",y="UMAP_2") +
    ggplot2::coord_fixed(xlim=limits$x,ylim=limits$y,expand=FALSE) + theme_sc_umap()
}
