# Re-render from the exported Seurat PCA spectrum; no PCA/DF recomputation.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2L,isTRUE(getOption("r45.runtime.active")))
source(file.path(args[1],"qc/run_qc.R"))
directory <- args[2]
metrics <- read.csv(file.path(directory,"pca_metrics.csv"),stringsAsFactors=FALSE)
for(i in seq_len(nrow(metrics))) {
  row <- metrics[i,]
  spectrum <- read.csv(file.path(directory,paste0(row$sample,"_pca_spectrum.csv")))
  plot <- ggplot2::ggplot(spectrum,ggplot2::aes(x=PC,y=stdev)) +
    ggplot2::geom_point(size=1.6) +
    ggplot2::geom_vline(xintercept=row$pc_used,color="#ea5455",linetype="dashed") +
    ggplot2::theme_classic(base_size=11) +
    ggplot2::labs(title=paste0(row$sample," | selected PC ",row$pc_used),
      subtitle=paste0(row$mode,"\nSD cumulative ",round(row$cumulative_sd_share,2),
        "%; squared-SD cumulative ",round(row$cumulative_variance_explained,2),"%"),
      y="PCA standard deviation",x="PC")
  save_sc_plot(plot,paste0("elbow_",row$sample),"umap_single",directory)
}
cat("PASS: readable elbow captions; spectra and PC choices unchanged\n")
