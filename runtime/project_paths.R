sc_paths_files <- unlist(lapply(sys.frames(),function(frame) frame$ofile),use.names=FALSE)
sc_paths_dir <- dirname(normalizePath(tail(sc_paths_files,1L),winslash="/",mustWork=TRUE))
sc_dataset_name <- function(value) {
  if(length(value)!=1L || is.na(value) || !grepl("^[A-Za-z0-9][A-Za-z0-9_.-]{0,99}$",value) ||
     grepl("[. ]$",value) || grepl("^(CON|PRN|AUX|NUL|COM[1-9]|LPT[1-9])($|\\.)",toupper(value)))
    stop("Unsafe dataset/run name")
  value
}
sc_flat_storage <- function(root) {
  declared <- Sys.getenv("SINGLE_CELL_OUTPUT_ROOT")
  configured <- jsonlite::fromJSON(file.path(sc_paths_dir,"output_config.json"))$output_root
  normalized <- tolower(normalizePath(root,winslash="/",mustWork=FALSE))
  normalized %in% tolower(vapply(c(configured,if(nzchar(declared)) declared else character()),
    normalizePath,character(1),winslash="/",mustWork=FALSE))
}
sc_project_path <- function(root,name) {
  sc_dataset_name(name)
  if(sc_flat_storage(root)) file.path(root,name) else file.path(root,"results",name)
}
