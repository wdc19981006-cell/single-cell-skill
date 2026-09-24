source(".agents/skills/geo-single-cell-loader/scripts/seurat_common.R")

root <- normalizePath(tempfile("geo-utf8-manifest-"), winslash="/", mustWork=FALSE)
workflow <- file.path(root, "data", "GSE999999998", ".workflow")
dir.create(workflow, recursive=TRUE)
dir.create(file.path(root, "data", "GSE999999998", "raw", "pooled"), recursive=TRUE)
on.exit(unlink(root, recursive=TRUE), add=TRUE)

labels <- c(intToUtf8(c(32925L, 32454L, 32990L, 30284L)),
            intToUtf8(c(32925L, 20869L, 32966L, 31649L, 30284L)))
header <- paste(c("database", "sample", "tissue", "disease", "source_type",
                  "group", "local_path", "file_type", "count_source",
                  "count_evidence", "metadata_evidence"), collapse=",")
row <- function(sample, group) paste(c("GSE999999998", sample, "Liver",
                                       "Liver cancer", "Tissue", group,
                                       "data/GSE999999998/raw/pooled",
                                       "10x_mtx", "counts", "GEO", "GEO"),
                                     collapse=",")
body <- paste(c(header, row("GSM1", labels[1]), row("GSM2", labels[2])), collapse="\n")
path <- file.path(workflow, "sample_manifest.csv")
receipt <- file.path(workflow, "sample_manifest.confirmation.json")

for (with_bom in c(FALSE, TRUE)) {
  bytes <- charToRaw(enc2utf8(paste0(body, "\n")))
  if (with_bom) bytes <- c(as.raw(c(239L, 187L, 191L)), bytes)
  con <- file(path, "wb")
  writeBin(bytes, con)
  close(con)
  confirmation <- list(
    confirmed_by="user", user_statement="User confirms disease groups",
    groups=list(GSM1=labels[1], GSM2=labels[2]),
    manifest_md5=unname(tools::md5sum(path))
  )
  writeLines(jsonlite::toJSON(confirmation, auto_unbox=TRUE, pretty=TRUE),
             receipt, useBytes=TRUE)
  actual <- read_manifest(path, root)
  stopifnot(nrow(actual) == 2L, identical(actual$group, labels),
            identical(names(actual)[1], "database"))
}
cat("UTF-8 manifest labels pass with and without BOM\n")
