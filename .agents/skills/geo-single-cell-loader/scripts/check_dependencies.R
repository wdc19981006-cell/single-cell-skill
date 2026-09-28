if (!isTRUE(getOption("r45.runtime.active"))) stop("Use runtime/r45/run_r45.py")
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
cat("Packages are never installed or changed by this workflow.\n")
cat("Fallback: set GEO_SINGLE_CELL_PYTHON to a verified Python containing anndata, numpy, scipy\n")
