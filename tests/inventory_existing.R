args <- commandArgs(trailingOnly=TRUE)
rows <- lapply(args, function(path) {
  x <- readRDS(path)
  data.frame(path=path, cells=ncol(x), genes=nrow(x),
    samples=paste(names(table(x$sample)),as.integer(table(x$sample)),sep=":",collapse=";"),
    groups=paste(names(table(x$group)),as.integer(table(x$group)),sep=":",collapse=";"),
    metadata=paste(names(x@meta.data),collapse=";"))
})
write.csv(do.call(rbind,rows),stdout(),row.names=FALSE)
