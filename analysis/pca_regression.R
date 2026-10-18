analysis_run_pca <- function(input,out,cfg,provenance) {
  analysis_require(c("Seurat","SeuratObject","Matrix","jsonlite","yaml","ggplot2","patchwork"))
  warnings <- analysis_warning_table()
  original <- readRDS(input)
  baseline <- analysis_pca(original,cfg,"pMT")
  baseline@misc$workflow_analysis$provenance <- provenance
  selection <- baseline@misc$workflow_analysis$pc_selection
  diagnostics <- analysis_regression_diagnostics(baseline,cfg)
  analysis_write_csv(diagnostics,file.path(out,"regression_diagnostics.csv"))
  analysis_write_csv(analysis_pc_selection_table(selection),file.path(out,"pca_selection.csv"))
  yaml::write_yaml(list(assay=cfg$assay,vars_to_regress="pMT",normalization="LogNormalize",
    nfeatures=cfg$nfeatures,seed=cfg$seed,pc_selection=selection),file.path(out,"pca_config.yaml"))
  elbow <- Seurat::ElbowPlot(baseline,ndims=selection$max_pc) + theme_sc() +
    ggplot2::geom_vline(xintercept=selection$pc_used,linetype=2,color="#B24332")
  analysis_plot(elbow,"elbow_plot","violin",out)
  heat <- ggplot2::ggplot(diagnostics,ggplot2::aes(x=pc,y=variable,fill=rho)) +
    ggplot2::geom_tile() + ggplot2::scale_fill_gradient2(limits=c(-1,1),na.value="grey85") +
    ggplot2::labs(title="Metadata versus selected PCs (Spearman)",x="PC",y=NULL) + theme_sc()
  panel <- analysis_embedding_panel(baseline,"sample","pca","Formal PCA: sample")
  library_panel <- analysis_embedding_panel(baseline,"nCount_RNA","pca","Formal PCA: RNA total")
  analysis_plot(patchwork::wrap_plots(panel,library_panel,heat,ncol=3),"pca_diagnostics","umap_3panel",out)
  if(selection$max_pc>=10L && selection$pc_used<=3L)
    warnings <- analysis_warning(warnings,"PARAMETER_UNCERTAINTY","Reference PC rule selected only 2-3 PCs despite >=10 available",
      "Review elbow and metadata relationships; do not force a fixed PC count")
  if(any(diagnostics$risk & diagnostics$variable!="pMT"))
    warnings <- analysis_warning(warnings,"PARAMETER_UNCERTAINTY","Metadata correlate with one or more selected PCs",
      "Correlation is diagnostic evidence, not proof of technical confounding; consider a targeted sensitivity comparison")
  if(cfg$assay!="RNA") warnings <- analysis_warning(warnings,"BIOLOGICAL_TRADEOFF",
    paste("Formal analysis explicitly uses",cfg$assay,"counts; original RNA counts remain available"),
    "Compare marker and sample structure against original RNA where biologically relevant",severity="INFO")
  requests <- unique(cfg$sensitivity)
  sensitivity_rows <- data.frame(target=character(),status=character(),variables=character(),max_abs_pc_correlation=double())
  for(target in requests) {
    extra <- if(target=="rna_total") "nCount_RNA" else c("S.Score","G2M.Score")
    rows <- diagnostics[diagnostics$variable %in% extra,,drop=FALSE]
    measurable <- any(is.finite(rows$rho))
    maxrho <- if(measurable) max(abs(rows$rho),na.rm=TRUE) else NA_real_
    risk <- measurable && maxrho>=cfg$risk_correlation
    status <- if(!all(extra %in% names(original@meta.data))) "SKIPPED_MISSING_METADATA" else
      if(!risk) "SKIPPED_NO_CLEAR_PC_RISK" else "COMPLETE"
    sensitivity_rows <- rbind(sensitivity_rows,data.frame(target=target,status=status,
      variables=paste(c("pMT",extra),collapse=";"),max_abs_pc_correlation=maxrho))
    if(status!="COMPLETE") next
    sensitivity_out <- file.path(out,"regression_sensitivity",target)
    dir.create(sensitivity_out,recursive=TRUE,showWarnings=FALSE)
    comparison <- analysis_pca(original,cfg,c("pMT",extra))
    comparison@misc$workflow_analysis$provenance <- provenance
    sensitivity_diagnostics <- analysis_regression_diagnostics(comparison,cfg)
    analysis_write_csv(sensitivity_diagnostics,file.path(sensitivity_out,"regression_diagnostics.csv"))
    analysis_write_csv(analysis_pc_selection_table(comparison@misc$workflow_analysis$pc_selection),file.path(sensitivity_out,"pca_selection.csv"))
    yaml::write_yaml(list(target=target,vars_to_regress=c("pMT",extra),trigger_abs_rho=maxrho,
      threshold=cfg$risk_correlation,baseline="pMT"),file.path(sensitivity_out,"sensitivity_config.yaml"))
    dims <- seq_len(min(selection$pc_used,comparison@misc$workflow_analysis$pc_selection$pc_used))
    # Compare loading subspaces without pretending PC signs/rotations are aligned.
    a <- SeuratObject::Embeddings(baseline[["pca"]])[,dims,drop=FALSE]
    b <- SeuratObject::Embeddings(comparison[["pca"]])[,dims,drop=FALSE]
    overlap <- sum(svd(crossprod(qr.Q(qr(a)),qr.Q(qr(b))),nu=0,nv=0)$d^2)/length(dims)
    analysis_write_csv(data.frame(metric="PC_subspace_overlap",value=overlap,
      note="Mean squared canonical cosine between PC subspaces; 1 means identical subspaces"),file.path(sensitivity_out,"sensitivity_comparison.csv"))
    plot <- patchwork::wrap_plots(analysis_embedding_panel(baseline,"sample","pca","Baseline: pMT"),
      analysis_embedding_panel(comparison,"sample","pca",paste("Sensitivity:",paste(extra,collapse=" + "))),ncol=2)
    analysis_plot(plot,"pca_comparison","umap_2panel",sensitivity_out)
    analysis_save(comparison,file.path(sensitivity_out,"seurat_sensitivity.rds"))
  }
  if(length(requests)) analysis_write_csv(sensitivity_rows,file.path(out,"sensitivity_summary.csv"))
  analysis_save(baseline,file.path(out,"seurat_pca.rds"))
  summary <- data.frame(metric=c("cells","genes","assay","baseline_regress","pc_used","selection_method"),
    value=c(ncol(baseline),nrow(baseline),cfg$assay,"pMT",selection$pc_used,selection$selection_method))
  analysis_finish("pca","COMPLETE",out,summary,warnings)
}
