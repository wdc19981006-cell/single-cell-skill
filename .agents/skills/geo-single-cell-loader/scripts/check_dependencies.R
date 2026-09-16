if (.Platform$OS.type == "windows" && !nzchar(Sys.getenv("PROCESSOR_ARCHITECTURE"))) {
  arch <- switch(R.version$arch, x86_64="AMD64", aarch64="ARM64", i386="x86", R.version$arch)
  Sys.setenv(PROCESSOR_ARCHITECTURE=arch)
}
cat(R.version.string,"\n")
routes <- list(
  base=c("Seurat","SeuratObject","Matrix","jsonlite"),
  text=c("data.table"),
  `10x_h5`=c("hdf5r"),
  h5ad_native=c("zellkonverter","SingleCellExperiment","SummarizedExperiment")
)
for (route in names(routes)) {
  cat("[",route,"]\n",sep="")
  for (p in routes[[route]]) cat(p,if(requireNamespace(p,quietly=TRUE)) as.character(packageVersion(p)) else "MISSING","\n")
}
cat("Optional packages are required only for their listed route; missing optional packages do not block other readers.\n")
cat("CRAN install: install.packages(c('Seurat','data.table','Matrix','jsonlite','hdf5r','BiocManager'))\n")
cat("Bioconductor install: BiocManager::install(c('zellkonverter','SingleCellExperiment'), ask=FALSE, update=FALSE)\n")
cat("Fallback: set GEO_SINGLE_CELL_PYTHON to a verified Python containing anndata, numpy, scipy\n")
