# Project-local R 4.5.0 profile. Never source this from another installation.
expected_home <- Sys.getenv("R45_EXPECTED_HOME")
expected_version <- Sys.getenv("R45_EXPECTED_VERSION")
actual_home <- normalizePath(R.home(), winslash="/", mustWork=TRUE)
fail_runtime <- function(message) {
  cat(message, "\n", file=stderr())
  quit(save="no", status=1L, runLast=FALSE)
}
if (!nzchar(expected_home) || !nzchar(expected_version) ||
    !startsWith(R.version.string, paste0("R version ", expected_version, " ")) ||
    !identical(tolower(actual_home), tolower(expected_home)) ||
    !identical(Sys.getenv("CLI_NO_THREAD"), "1")) {
  fail_runtime("R45_PROFILE_MISMATCH: designated R and CLI_NO_THREAD=1 required")
}
expected_library <- Sys.getenv("R45_EXPECTED_LIBRARY")
.libPaths(expected_library)
if (!identical(tolower(normalizePath(.libPaths(), winslash="/")),
               tolower(normalizePath(expected_library, winslash="/")))) {
  fail_runtime("R45_LIBRARY_MISMATCH")
}
if (!is.null(getOption("error"))) fail_runtime("R45_PROFILE_ERROR_HANDLER_CONFLICT")
if (!nzchar(Sys.getenv("PROCESSOR_ARCHITECTURE"))) {
  actual_arch <- switch(R.version$arch, x86_64="AMD64", aarch64="ARM64", i386="x86", "")
  if (!nzchar(actual_arch)) fail_runtime("R45_UNSUPPORTED_ARCHITECTURE")
  Sys.setenv(PROCESSOR_ARCHITECTURE=actual_arch)
}
options(r45.runtime.active=TRUE, r45.cli_exit_workaround=TRUE)

# The architecture seen by calculations is unchanged. Only cli's process-exit
# cleanup branch is selected after calculations have finished.
r45_exit_cleanup <- function() {
  if ("cli" %in% loadedNamespaces()) {
    Sys.setenv(PROCESSOR_ARCHITECTURE="ARM64")
  }
}
.Last <- r45_exit_cleanup
options(error=function() {
  r45_exit_cleanup()
  quit(save="no", status=1L)
})
