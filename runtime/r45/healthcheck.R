args <- commandArgs(trailingOnly=TRUE)
if (length(args) > 1L) stop("Usage: healthcheck.R [base|qc|text|10x_h5|h5ad_native]")
route <- if (length(args)) args[1] else "base"
routes <- list(base=character(), qc="DoubletFinder", text="data.table",
               `10x_h5`="hdf5r",
               h5ad_native=c("zellkonverter", "SingleCellExperiment", "SummarizedExperiment"))
if (!route %in% names(routes)) stop("Unknown healthcheck route: ", route)
expected_home <- Sys.getenv("R45_EXPECTED_HOME")
expected_library <- Sys.getenv("R45_EXPECTED_LIBRARY")
if (!isTRUE(getOption("r45.runtime.active")) ||
    !startsWith(R.version.string, paste0("R version ", Sys.getenv("R45_EXPECTED_VERSION"), " ")) ||
    !identical(tolower(normalizePath(R.home(), winslash="/")), tolower(expected_home)) ||
    !tolower(expected_library) %in% tolower(normalizePath(.libPaths(), winslash="/"))) {
  stop("R_RUNTIME_MISMATCH")
}
required <- c("Seurat", "SeuratObject", "Matrix", "jsonlite", "yaml", "ggplot2", "patchwork")
packages <- unique(c(required, routes[[route]]))
for (package in packages) {
  version <- if (requireNamespace(package, quietly=TRUE))
    as.character(utils::packageVersion(package)) else "MISSING"
  cat("R45_PACKAGE\t", package, "\t", version, "\n", sep="")
}
cat("R45_VERSION\t", R.version.string, "\n", sep="")
cat("R45_SHORT_VERSION\t", as.character(getRversion()), "\n", sep="")
cat("R45_HOME\t", normalizePath(R.home(), winslash="/"), "\n", sep="")
cat("R45_LIBRARIES\t", paste(.libPaths(), collapse=";"), "\n", sep="")
cat("R45_ROUTE\t", route, "\n", sep="")
missing <- packages
missing <- missing[!vapply(missing, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) stop("R_PACKAGE_MISSING: ", paste(missing, collapse=", "))
cat("R45_HEALTH_PASS\n")
