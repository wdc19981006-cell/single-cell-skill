if (.Platform$OS.type == "windows" && !nzchar(Sys.getenv("PROCESSOR_ARCHITECTURE"))) {
  arch <- switch(R.version$arch, x86_64="AMD64", aarch64="ARM64", i386="x86", R.version$arch)
  Sys.setenv(PROCESSOR_ARCHITECTURE=arch)
}
packages <- c("Seurat","SeuratObject","dplyr","data.table","Matrix","jsonlite","hdf5r","zellkonverter","SingleCellExperiment","SummarizedExperiment")
cat(R.version.string,"\n")
for (p in packages) cat(p,if(requireNamespace(p,quietly=TRUE)) as.character(packageVersion(p)) else "MISSING","\n")
cat("CRAN install: install.packages(c('Seurat','dplyr','data.table','Matrix','jsonlite','hdf5r','BiocManager'))\n")
cat("Bioconductor install: BiocManager::install(c('zellkonverter','SingleCellExperiment'), ask=FALSE, update=FALSE)\n")
cat("Fallback: set GEO_SINGLE_CELL_PYTHON to a verified Python containing anndata, numpy, scipy\n")
