# Read complete RDS objects sequentially; compare expression in sparse cell blocks.
# Usage: run_r45.py compare_objects.R INPUT_CSV OUTPUT_DIRECTORY
args <- commandArgs(trailingOnly=TRUE)
stopifnot(isTRUE(getOption("r45.runtime.active")), length(args)==2L)
inputs <- read.csv(args[1],check.names=FALSE,stringsAsFactors=FALSE,encoding="UTF-8")
out <- args[2]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
rows <- list()
add <- function(gse,metric,a,b,status="IDENTICAL",note="") {
  rows[[length(rows)+1L]] <<- data.frame(dataset=gse,metric=metric,
    original=paste(a,collapse=";"),skill=paste(b,collapse=";"),classification=status,note=note)
}
write_report <- function() {
  result <- do.call(rbind,rows)
  write.csv(result,file.path(out,"object_comparison.csv"),row.names=FALSE,na="")
  lines <- c("# Seurat object comparison", "",
    "Complete RDS reads. Historical inputs are verified against data-QC Git LFS SHA-256; new QC outputs have separate raw-count and artifact validation. Expression is compared exactly in sparse blocks of 512 cells; no dense expression conversion.", "",
    "EXPECTED_DIFFERENCE denotes a known strategy/state difference, not numerical equivalence. UNEXPECTED_DIFFERENCE requires investigation; a plausible explanation alone is not proof.", "",
    "Reference code filters DF; Skill filters DF_adj. Reference code repeats a gene prevalence filter (>3 cells); Skill preserves raw genes. Reference code regresses nFeature_RNA and pMT; Skill regresses pMT only. The reference code's intended final step retains normalized data/PCA and cell-cycle scores, but the supplied original RDS artifacts actually contain only counts and no cell-cycle scores. This artifact/protocol mismatch is UNEXPECTED_DIFFERENCE, not a reason to change Skill behavior. Cell-cycle cross-version concordance is therefore unavailable.", "",
    "| Dataset | Metric | Original | Skill | Classification | Evidence / interpretation |",
    "|---|---|---|---|---|---|")
  for(i in seq_len(nrow(result))) {
    r <- result[i,]
    brief <- function(x) {x <- gsub("|","/",x,fixed=TRUE); if(nchar(x)>220) paste0(substr(x,1,220),"... (full CSV)") else x}
    lines <- c(lines,paste0("| ",paste(vapply(r,brief,character(1)),collapse=" | ")," |"))
  }
  writeLines(lines,file.path(out,"object_comparison.md"),useBytes=TRUE)
}
for(gse in unique(inputs$dataset)) {
  records <- inputs[inputs$dataset==gse,]
  cat("READING",gse,"\n"); flush.console()
  a <- readRDS(records$path[records$version=="original"])
  skill_record <- records[records$version!="original",]
  stopifnot(nrow(skill_record)==1L)
  b <- readRDS(skill_record$path)
  stopifnot(inherits(a,"Seurat"),inherits(b,"Seurat"))
  structures <- list(class=list(class(a),class(b)),cells=list(ncol(a),ncol(b)),genes=list(nrow(a),nrow(b)),
    assays=list(names(a@assays),names(b@assays)),RNA_assay_class=list(class(a[["RNA"]]),class(b[["RNA"]])),
    RNA_layers=list(SeuratObject::Layers(a[["RNA"]]),SeuratObject::Layers(b[["RNA"]])),
    VariableFeatures=list(SeuratObject::VariableFeatures(a),SeuratObject::VariableFeatures(b)),
    reductions=list(names(a@reductions),names(b@reductions)),graphs=list(names(a@graphs),names(b@graphs)),
    neighbors=list(names(a@neighbors),names(b@neighbors)),metadata_columns=list(names(a@meta.data),names(b@meta.data)))
  for(metric in names(structures)) {
    x <- structures[[metric]]
    add(gse,metric,x[[1]],x[[2]],if(identical(x[[1]],x[[2]])) "IDENTICAL" else "EXPECTED_DIFFERENCE",
      "Structural/filter strategy difference; see quantitative checks below")
  }
  if(!all(c("S.Score","G2M.Score","Phase","CC.Difference") %in% names(a@meta.data)) ||
     !"pca" %in% names(a@reductions) || !"data" %in% SeuratObject::Layers(a[["RNA"]]))
    add(gse,"original_artifact_vs_reference_final_step","Missing expected final analysis state","Counts-only final by design",
        "UNEXPECTED_DIFFERENCE","The actual original RDS does not contain the final cell-cycle/PCA state described by reference code; no values inferred or fabricated")
  cells <- intersect(colnames(a),colnames(b)); genes <- intersect(rownames(a),rownames(b))
  for(kind in c("cell","gene")) {
    x <- if(kind=="cell") colnames(a) else rownames(a)
    y <- if(kind=="cell") colnames(b) else rownames(b)
    add(gse,paste0(kind,"_intersection"),length(intersect(x,y)),length(intersect(x,y)))
    only_a <- setdiff(x,y); only_b <- setdiff(y,x)
    add(gse,paste0(kind,"_exclusive"),length(only_a),length(only_b),
      if(length(only_a)+length(only_b)) "EXPECTED_DIFFERENCE" else "IDENTICAL",
      if(kind=="cell") "Different DF/DF_adj filters and PCA/regression; exclusives do not prove cause" else "Original extra gene prevalence filter; verify using raw input where available")
    write.csv(data.frame(version=c(rep("original",length(only_a)),rep(skill_record$version,length(only_b))),id=c(only_a,only_b)),
      file.path(out,paste0(gse,"_",kind,"_differences.csv")),row.names=FALSE)
  }
  stopifnot(length(cells)>0,length(genes)>0)
  ac <- SeuratObject::LayerData(a,assay="RNA",layer="counts")
  bc <- SeuratObject::LayerData(b,assay="RNA",layer="counts")
  stopifnot(inherits(ac,"sparseMatrix"),inherits(bc,"sparseMatrix"))
  mismatches <- 0; max_delta <- 0
  calculated <- list(nCount_RNA=list(a=numeric(length(cells)),b=numeric(length(cells))),
                     nFeature_RNA=list(a=numeric(length(cells)),b=numeric(length(cells))))
  for(start in seq.int(1L,length(cells),by=512L)) {
    indexes <- start:min(start+511L,length(cells))
    block <- cells[indexes]
    delta <- ac[genes,block,drop=FALSE]-bc[genes,block,drop=FALSE]
    mismatches <- mismatches+sum(delta@x!=0)
    if(length(delta@x)) max_delta <- max(max_delta,abs(delta@x))
    for(version in c("a","b")) {
      full_block <- if(version=="a") ac[,block,drop=FALSE] else bc[,block,drop=FALSE]
      calculated$nCount_RNA[[version]][indexes] <- Matrix::colSums(full_block)
      calculated$nFeature_RNA[[version]][indexes] <- Matrix::colSums(full_block>0)
    }
  }
  add(gse,"counts_mismatched_entries",0,mismatches,if(mismatches==0) "IDENTICAL" else "UNEXPECTED_DIFFERENCE",
    paste0("Exact sparse comparison on ",length(genes)," genes x ",length(cells)," cells; max_abs_delta=",max_delta))
  for(field in c("sample","group")) {
    ta <- table(a@meta.data[[field]]); tb <- table(b@meta.data[[field]])
    for(level in union(names(ta),names(tb))) {
      va <- if(level %in% names(ta)) unname(ta[level]) else 0L
      vb <- if(level %in% names(tb)) unname(tb[level]) else 0L
      add(gse,paste(field,level,sep=":"),va,vb,if(va==vb) "IDENTICAL" else "EXPECTED_DIFFERENCE","Different final retained cell sets")
    }
  }
  for(field in c("sample","group","nFeature_RNA","nCount_RNA","pMT","pRP","pHB","DF","DF_adj","S.Score","G2M.Score","Phase","CC.Difference")) {
    if(!field %in% names(a@meta.data) || !field %in% names(b@meta.data)) {
      add(gse,field,field %in% names(a@meta.data),field %in% names(b@meta.data),"UNEXPECTED_DIFFERENCE","Required field missing")
      next
    }
    x <- a@meta.data[cells,field]; y <- b@meta.data[cells,field]
    numeric <- is.numeric(x) && is.numeric(y)
    same <- if(numeric) abs(x-y)<=1e-10 else as.character(x)==as.character(y)
    same[is.na(same)] <- is.na(x[is.na(same)]) & is.na(y[is.na(same)])
    count <- sum(!same)
    expected <- field %in% c("DF","DF_adj","S.Score","G2M.Score","Phase","CC.Difference")
    note <- if(numeric) paste0("common-cell mismatches (tol 1e-10); max_abs_delta=",max(abs(x-y),na.rm=TRUE)) else "common-cell category mismatches"
    if(field %in% names(calculated)) {
      proven <- all(x==calculated[[field]]$a) && all(y==calculated[[field]]$b) && mismatches==0
      expected <- proven
      note <- paste0(note,if(proven) "; both metadata vectors exactly equal sparse per-object counts statistics; differences explained by nonshared gene rows" else "; metadata/count reconciliation failed")
    }
    add(gse,field,0,count,if(count==0) "IDENTICAL" else if(expected) "EXPECTED_DIFFERENCE" else "UNEXPECTED_DIFFERENCE",note)
    if(field %in% c("DF","DF_adj","Phase")) {
      write.csv(as.data.frame(table(original=x,skill=y,useNA="ifany")),file.path(out,paste0(gse,"_",field,"_cross.csv")),row.names=FALSE)
    }
  }
  for(field in c("DF","DF_adj")) add(gse,paste0(field,"_all_cells"),
    paste(names(table(a@meta.data[[field]])),table(a@meta.data[[field]]),sep="="),
    paste(names(table(b@meta.data[[field]])),table(b@meta.data[[field]]),sep="="),"EXPECTED_DIFFERENCE","Final objects already filtered: cannot reconstruct all discarded calls")
  write_report()
  cat("COMPARED",gse,"cells",length(cells),"genes",length(genes),"counts mismatches",mismatches,"\n"); flush.console()
  rm(a,b,ac,bc,delta); gc()
}
