# Read-only report renderer. All scientific evidence comes from saved tables.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2L) stop("Usage: render_audit.R AUDIT_DIRECTORY GSE")
area <- normalizePath(args[1],winslash="/",mustWork=TRUE)
repo <- normalizePath(Sys.getenv("R_PROFILE_USER"),winslash="/",mustWork=TRUE)
repo <- dirname(dirname(dirname(repo)))
source(file.path(repo,"qc/plot_utils.R"))
audit <- read.csv(file.path(area,"workflow_audit.csv"),stringsAsFactors=FALSE,check.names=FALSE)
warnings <- read.csv(file.path(area,"workflow_warnings.csv"),stringsAsFactors=FALSE,check.names=FALSE)
output <- file.path(area,"workflow_audit.pdf")
if(file.exists(output)) stop("Refusing to overwrite existing audit PDF")
grDevices::pdf(output,width=11,height=8.5,onefile=TRUE)
page <- 0L
text_page <- function(title,lines) {
  chunks <- split(lines,ceiling(seq_along(lines)/31L))
  if(!length(chunks)) chunks <- list("No observed records.")
  for(chunk in chunks) {
    page <<- page+1L
    grid::grid.newpage()
    grid::grid.text(title,x=.06,y=.94,just=c("left","top"),gp=grid::gpar(fontsize=16,fontface="bold"))
    grid::grid.text(paste(chunk,collapse="\n"),x=.06,y=.86,just=c("left","top"),
                    gp=grid::gpar(fontsize=9,lineheight=1.4,fontfamily="sans"))
    grid::grid.text(paste(args[2],"| existing evidence only |",page),x=.94,y=.04,just="right",gp=grid::gpar(fontsize=8,col="grey40"))
  }
}
wrap <- function(x) unlist(lapply(x,function(s) {
  # Paths and serialized vectors may contain no whitespace. Break those
  # tokens before wrapping so no evidence runs beyond the page boundary.
  s <- gsub("([^[:space:]]{80})", "\\1 ", s, perl=TRUE)
  strwrap(s,width=112)
}))
stages <- c("download","qc1","qc2","qc3","pca","harmony","clustering")
overview <- c(paste("Generated UTC:",format(Sys.time(),tz="UTC",usetz=TRUE)),
  "This audit summarizes saved evidence. File presence alone is not scientific validation.",
  "Missing evidence and biological decisions remain explicit. No analysis was rerun.","",
  vapply(stages,function(s) paste(s,":",sum(audit$stage==s),"observed fields;",sum(warnings$stage==s & !warnings$status %in% c("RESOLVED","PASS")),"current review items"),character(1)))
text_page("Single-cell workflow audit",wrap(overview))
for(stage in stages) {
  rows <- audit[audit$stage==stage,,drop=FALSE]
  # Large tables stay in the CSV. The PDF favors sample summaries and
  # diagnostic/decision records, with the exact provenance retained in CSV.
  lines <- character()
  for(src in unique(rows$source_file)) {
    table <- rows[rows$source_file==src,,drop=FALSE]
    ids <- unique(table$row)
    # Show strongest measured PC associations first, retaining variable,
    # PC, rho, risk and status together rather than disconnected fields.
    rho <- table[table$metric=="rho",,drop=FALSE]
    if(nrow(rho)) ids <- unique(c(rho$row[order(abs(as.numeric(rho$value)),decreasing=TRUE)],ids))
    lines <- c(lines,paste("Source:",basename(src)))
    for(id in head(ids,24L)) {
      record <- table[table$row==id,,drop=FALSE]
      values <- paste(record$metric,record$value,sep="=")
      lines <- c(lines,paste(record$sample[1],paste(values,collapse="; "),sep=": "))
    }
    if(length(ids)>24L) lines <- c(lines,paste(length(ids)-24L,"additional records retained in workflow_audit.csv."))
    lines <- c(lines,"")
  }
  if(!length(lines)) lines <- "No saved summary metrics. Consult workflow_warnings.csv."
  text_page(paste("Observed evidence:",stage),wrap(lines))
}
if(nrow(warnings)) {
  current <- warnings[!warnings$status %in% c("RESOLVED","PASS"),,drop=FALSE]
  keys <- unique(paste(current$stage,current$category,current$status,sep=" | "))
  lines <- unlist(lapply(keys,function(key) {
    rows <- current[paste(current$stage,current$category,current$status,sep=" | ")==key,,drop=FALSE]
    c(paste0("[",key,"] ",nrow(rows)," record(s)"),head(rows$message,3L),
      if(nrow(rows)>3L) "Additional records preserved in workflow_warnings.csv",paste("Action:",rows$recommendation[1]),"")
  }))
  text_page("Current warnings and required decisions",wrap(if(length(lines)) lines else "No unresolved warnings."))
  historical <- warnings[warnings$status=="RESOLVED",,drop=FALSE]
  if(nrow(historical)) text_page("Resolved historical decisions",wrap(paste("[RESOLVED]",unique(historical$message))))
} else text_page("Warnings and required review","No warnings were present in the collected tables; this is not proof that all biological risks are absent.")
grDevices::dev.off()
cat("Audit PDF pages:",page,"\n")
