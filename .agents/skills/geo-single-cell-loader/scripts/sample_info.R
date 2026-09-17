# Preserve the user report and append results from the validated Seurat object.
update_sample_info <- function(out,status,m=NULL,object=NULL,inputs=character(),warnings=character(),failure=NULL) {
  path <- file.path(out,"sample_info.txt")
  lines <- if(file.exists(path)) readLines(path,encoding="UTF-8",warn=FALSE) else c(basename(out),"raw/ contains verified source expression files and exact byte reconstructions, not necessarily sequencing FASTQ raw reads.")
  # Remove prior status/failure/result footer on retries, leaving sample facts intact.
  headings <- which(lines %in% c("PROCESSING SUMMARY","Failure:"))
  if(length(headings)) {
    start <- min(headings)
    if(start > 1L && lines[start-1L] == strrep("=",60)) start <- start-1L
    lines <- if(start > 1L) lines[seq_len(start-1L)] else character()
  }
  lines <- lines[!grepl("^STATUS:",lines)]
  lines <- c(paste("STATUS:",status),lines)
  if(!is.null(failure)) lines <- c(lines,"","Failure:",failure)
  if(status == "COMPLETE") {
    # Manifest is the sole source for sample/group values in the final report.
    sample_heading <- which(lines == "SAMPLE INFORMATION")
    if(length(sample_heading)) lines <- lines[seq_len(max(1L,sample_heading[1]-2L))]
    lines <- c(lines,strrep("=",60),"SAMPLE INFORMATION",strrep("=",60),"")
    optional <- intersect(optional_fields,names(m))
    optional <- optional[vapply(m[optional],function(x) any(!blank(x)),logical(1))]
    for(i in seq_len(nrow(m))) {
      lines <- c(lines,paste("Sample",i),strrep("-",60))
      for(field in c("sample","author_sample","tissue","disease","source_type","sample_description",optional,"group")) lines <- c(lines,paste0(field,":"),value(m[i,,drop=FALSE],field,"Not reported"),"")
    }
    lines <- c(lines,strrep("=",60),"GROUP SUMMARY",strrep("=",60),"")
    for(group in unique(m$group)) lines <- c(lines,group,strrep("-",60),m$sample[m$group == group],"")
    lines <- c(lines,strrep("=",60),"PROCESSING SUMMARY",strrep("=",60),"",
      "Downloaded source data:","raw/","","Processed Seurat object:","seurat_raw.rds","",
      "Total samples:",as.character(nrow(m)),"","Total cells:",as.character(ncol(object)),"","Total genes:",as.character(nrow(object)),"",
      "Construction:","All samples merged at count-matrix level before CreateSeuratObject.","",
      "CreateSeuratObject:","min.cells = 3","min.features = 200","",
      "Additional QC:","None","","Normalization:","None","","Clustering:","None","",
      "Created:",format(Sys.time(),"%Y-%m-%d %H:%M %Z"),"","Seurat version:",as.character(packageVersion("Seurat")),"",
      "Cells per sample:",capture.output(table(object$sample)),"","Cells per group:",capture.output(table(object$group)),"",
      "Reading formats:",paste(unique(m$file_type),collapse=", "),"","Counts source:",inputs,"",
      "Warnings:",if(length(warnings)) unique(warnings) else "None")
  }
  writeLines(enc2utf8(lines),path,useBytes=TRUE)
}
