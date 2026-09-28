if (!isTRUE(getOption("r45.runtime.active"))) stop("Use runtime/r45/run_r45.py")
cat(R.version.string, "\n", R.home(), "\n", paste(.libPaths(), collapse=";"), "\n", sep="")
