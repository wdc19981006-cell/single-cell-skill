# Batch-only workaround for cli's native thread cleanup on this Windows R 4.5.0 host.
# The real architecture is retained throughout analysis; only .Last changes
# the cleanup branch. CLI_NO_THREAD prevents starting the background timer.
if (!identical(R.version.string, "R version 4.5.0 (2025-04-11 ucrt)") ||
    !identical(normalizePath(R.home(), winslash="/"), "D:/R/R-4.5.0") ||
    !identical(Sys.getenv("CLI_NO_THREAD"), "1")) {
  stop("R 4.5.0 QC exit profile requires the designated R and CLI_NO_THREAD=1")
}

.Last <- function() {
  if ("cli" %in% loadedNamespaces()) {
    Sys.setenv(PROCESSOR_ARCHITECTURE="ARM64")
  }
}

# Rscript does not call .Last after an unhandled error. Preserve the failure
# exit code while selecting the same cleanup branch before termination.
if (!is.null(getOption("error"))) {
  stop("R 4.5.0 QC exit profile will not replace an existing error handler")
}
options(error=function() {
  if ("cli" %in% loadedNamespaces()) {
    Sys.setenv(PROCESSOR_ARCHITECTURE="ARM64")
  }
  quit(save="no",status=1L)
})
