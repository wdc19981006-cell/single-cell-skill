analysis_run_pca <- function(input,out,cfg,provenance) {
  analysis_require(c("Seurat","SeuratObject","Matrix","jsonlite","yaml","ggplot2","patchwork"))
  warnings <- analysis_warning_table()
  original <- readRDS(input)
  regress <- if(is.null(cfg$regression_decision)) "pMT" else c("pMT",switch(cfg$regression_decision$scheme,
    feature="nFeature_RNA",count="nCount_RNA",cell_cycle=c("S.Score","G2M.Score")))
  baseline <- analysis_pca(original,cfg,regress)
  baseline@misc$workflow_analysis$provenance <- provenance
  selection <- baseline@misc$workflow_analysis$pc_selection
  diagnostics <- analysis_regression_diagnostics(baseline,cfg)
  analysis_write_csv(diagnostics,file.path(out,"regression_diagnostics.csv"))
  analysis_write_csv(analysis_pc_selection_table(selection),file.path(out,"pca_selection.csv"))
  yaml::write_yaml(list(assay=cfg$assay,vars_to_regress=regress,normalization="LogNormalize",
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
  analysis_save(baseline,file.path(out,"seurat_pca.rds"))
  summary <- data.frame(metric=c("cells","genes","assay","baseline_regress","pc_used","selection_method"),
    value=c(ncol(baseline),nrow(baseline),cfg$assay,paste(regress,collapse=";"),selection$pc_used,selection$selection_method))
  analysis_finish("pca","COMPLETE",out,summary,warnings,cfg,provenance)
}
