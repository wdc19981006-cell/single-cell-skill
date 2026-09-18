args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if (length(args)) args[1] else ".", winslash="/", mustWork=TRUE)
source(file.path(root, ".agents/skills/geo-single-cell-loader/scripts/seurat_common.R"))
need(c("Seurat", "Matrix"))

fixture_root <- tempfile("tenx-literal-na-")
dir.create(file.path(fixture_root, "data/GSE999999999/raw"), recursive=TRUE)
on.exit(unlink(fixture_root, recursive=TRUE), add=TRUE)
fixture_root <- normalizePath(fixture_root, winslash="/", mustWork=TRUE)

make_trio <- function(name, middle) {
  directory <- file.path(fixture_root, "data/GSE999999999/raw", name)
  dir.create(directory)
  counts <- Matrix::Matrix(matrix(c(1, 0, 2, 0, 3, 4), nrow=3), sparse=TRUE)
  Matrix::writeMM(counts, file.path(directory, "matrix.mtx"))
  writeLines(c("ID1\tGeneA", middle, "ID3\tGeneC"), file.path(directory, "features.tsv"))
  writeLines(c("cell1", "cell2"), file.path(directory, "barcodes.tsv"))
  data.frame(database="GSE999999999", sample="GSM1",
             local_path=paste0("data/GSE999999999/raw/", name),
             file_type="10x_mtx", count_source="counts", stringsAsFactors=FALSE)
}

literal_row <- make_trio("literal_na", "NA\tNA")
literal_na <- read_expression(literal_row, fixture_root)
stopifnot(identical(rownames(literal_na$counts), c("GeneA", "unannotated_feature_2", "GeneC")),
          identical(dim(literal_na$counts), c(3L, 2L)),
          identical(as.numeric(literal_na$counts[2, ]), c(0, 3)))

ordinary <- read_expression(make_trio("ordinary", "ID2\tGeneB"), fixture_root)
stopifnot(identical(rownames(ordinary$counts), c("GeneA", "GeneB", "GeneC")))
cat("PASS: placeholder 10x feature IDs receive stable names while ordinary trios remain unchanged.\n")
