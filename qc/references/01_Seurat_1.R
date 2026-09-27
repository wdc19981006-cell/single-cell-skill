#单细胞数据处理上Seurat分析流 基础概念介绍，质控，细胞周期影响及双细胞去除
#Guan
rm(list=ls())
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(cowplot)
library(viridis)
library(patchwork)
library(tibble)
#rownames(seurat@meta.data) <- Cells(seurat)
#一、数据读取----
##标准10x三联文件读取
path<-'GSE231993'#设置一个路径
dir <- list.dirs(path)[-1] #list.dirs会把根目录也列出来，所以要-1. （读取这个样本下面的所有文件名字）
names(dir) <- list.files(path,recursive = F)
seurat<-CreateSeuratObject(counts = Read10X(dir),min.cells = 3, min.features = 200)
#读取H5文件直接变成Read10X_h5,即可
seurat$database<-rep(path,ncol(seurat))
#对seurat取子集相当于对meta取子集（database就是说明样本，（这里指的是细胞）来自哪个的数据集，可能方便后面的多数据集
seurat$sample<-seurat$orig.ident
#类似于提示每个细胞来自哪个样本，注意不是数据集
saveRDS(seurat,'seurat_raw.rds')

table(Idents(seurat))#检查每个样本中的细胞
meta = seurat@meta.data
#下面是创建分组信息
# 1. 检查当前 meta.data 中的样本列（假设列名为 'sample'）
head(seurat@meta.data$sample)  # 确认样本编号是否正确
# 2. 添加分组信息（根据样本编号判断 Tumor/Control）
seurat$group <- ifelse(seurat$sample %in% c("1", "2", "3"), "Tumor", "Control")
# 3. 检查是否添加成功
table(seurat$group)  # 应该输出 Tumor 和 Control 的细胞数


#类似于读取csv或者tst文件的数据集
scelist = list()
for (i in 1:length(dir)) {
  count = read.table(paste0(dir[i],"matrix.txt"),header = T)
  meta = read.table(paste0(dir[i],"meta.txt"),header = T)
  sce = CreateSeuratObject(counts = count,min.cells = 3, min.features = 200,meta.data = meta)
  scelist[[i]] = sce
}
seurat = merge(scelist[[1]],scelist[[-1]],merge.data = T)
seurat = JoinLayers(seurat)
#JoinLayers合并layers
x = colnames(seurat)
y = rownames(seurat)
#取子集对meta操作，col row 对表达矩阵操作
str_split()

seurat = AddMetadata(seurat)
#创建分组信息
group = data.frame(dir = dir,
                   group = c(rep("HC",6),rep("CA",6)))
group1 = str_split(dir,"/",simplify = T)
#删掉dir
group$dir = group1[,2]
meta = seurat@meta.data
meta1 = merge(meta,group,by.x = "sample",by.y = "dir")
seurat$group = meta1$group


sce1$group = "HC"
sce2$group = "CA"
sce3

sce1$database = "GSE123435"
sce2$database = "GSE343524"
sce3

scelist = c(sce1,sce2,sce3)
scelist = list()

seurat = merge(scelist[[1]],scelist[[-1]],merge.data = T)
seurat = JoinLayers(seurat)

#更新seurat（用于subset以后去掉冗余metadata及其他冗余slot时）
if(F){
  seurat<-CreateSeuratObject(counts=LayerData(seurat,layers='counts'),
                             meta.data =FetchData(seurat,vars = c('celltype','group','sample')) )
}

#二、质控-----
getwd()
#获得当下的工作目录
#人类
if(T){
  seurat=PercentageFeatureSet(seurat, "^MT-", col.name = "pMT")
  seurat=PercentageFeatureSet(seurat, "^RP[SL]", col.name = "pRP")
  seurat=PercentageFeatureSet(seurat, "^HB[^(P)]", col.name = "pHB")
}
#PercentageFeatureSet基因的表达比例
#小鼠
if(F){
  seurat=PercentageFeatureSet(seurat, "^mt-", col.name = "pMT")
  seurat=PercentageFeatureSet(seurat, "^Rp[sl]", col.name = "pRP")
  seurat=PercentageFeatureSet(seurat, "^Hb[^(p)]", col.name = "pHB")
}

#质控小提琴图：
feats <- c("pMT", "pRP", "pHB","nFeature_RNA", "nCount_RNA")
VlnPlot(seurat, group.by = "sample", features = feats, pt.size = 0, ncol = 3) &
  NoLegend() &
  theme(legend.position="none",#可以用ggplot的逻辑改变坐标轴摆放方式：
        plot.title=element_text(hjust=.5,size=16,face="bold"),
        axis.text.x=element_text(angle=50,size=8,vjust=1,hjust=1))
ggsave(filename="before_QC.pdf",height = 8, width = 20)


#细胞表达的基因在200-4000，mRNA数在500-30000之间
#UMI指的是一条mRNA,也就是ncount
nFeature_lower <- 200
nFeature_upper <- 4000
nCount_lower <- 500
nCount_upper <- 30000
pMT_upper <- 25
pRP_upper <- 100
pHB_upper<-1
selected_c <- WhichCells(seurat, expression =
                           nFeature_RNA > nFeature_lower & nFeature_RNA < nFeature_upper
                         & nCount_RNA > nCount_lower & nCount_RNA < nCount_upper )
selected_mito <- WhichCells(seurat, expression = pMT < pMT_upper )
selected_ribo <- WhichCells(seurat, expression = pRP < pRP_upper)
selected_hb <- WhichCells(seurat, expression = pHB < pHB_upper)

seurat.filt <- subset(seurat, cells = c(selected_c)) %>%
  subset(cells = c(selected_hb)) %>%
  subset(cells = c(selected_mito))%>%
  subset(cells = c(selected_ribo))
#过滤3%-5%数据很好，10%-15%可以接受
#上面是筛选细胞，接下来筛选基因
f<-JoinLayers(seurat)[["RNA"]]$counts
selected_f <- rownames(seurat)[Matrix::rowSums(f > 0 ) > 3]
#基因至少在3个细胞中表达
seurat.filt <- subset(seurat.filt, features = selected_f)

feats <- c("pMT", "pRP", "pHB","nFeature_RNA", "nCount_RNA")
VlnPlot(seurat.filt, group.by = "sample", features = feats, pt.size = 0, ncol = 3) &
  NoLegend() &
  theme(legend.position="none",#可以用ggplot的逻辑改变坐标轴摆放方式：
        plot.title=element_text(hjust=.5,size=16,face="bold"),
        axis.text.x=element_text(angle=50,size=8,vjust=1,hjust=1))
ggsave(filename="after_QC.pdf",height = 8, width = 20)

FeatureScatter(seurat, feature1 = 'nFeature_RNA', feature2 = 'nCount_RNA')
FeatureScatter(seurat, feature1 = 'nFeature_RNA', feature2 = 'nCount_RNA')

seurat<-seurat.filt
rm(seurat.filt)
saveRDS(seurat, "seurat.qc.rds")
#三、双细胞-----
library(tidyverse)
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(cowplot)
library(viridis)
library(patchwork)
library(tibble)
library(DoubletFinder)
library(tidydr)
library(randomcoloR)
library(SCP)
###DoubletFinder算法----
#安装
#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
sample.ident='sample'
sample.cols=distinctColorPalette(unique(seurat@meta.data[,sample.ident])%>%length)
###核心关键点：确保数据经过了标准质控后
###样本分割进行双细胞去除
seurat = JoinLayers(seurat)

# 双细胞只能单一样本跑，不可以整体跑
# 检测Doublets----先找高变、归一化数据，再RunPCA，跑完UMAP；才能计算双细胞比例
# DoubletFinder的计算需要利用PCA降维信息
## 去除双细胞前的准备
seurat.list <- SplitObject(seurat, split.by = sample.ident)
for (i in 1:length(seurat.list)) {
  print(i)
  seurat.list[[i]] <- NormalizeData(seurat.list[[i]])
  seurat.list[[i]] <- FindVariableFeatures(seurat.list[[i]])
  seurat.list[[i]] <- ScaleData(seurat.list[[i]], vars.to.regress = c("nFeature_RNA", "pMT"))
  seurat.list[[i]] = RunPCA(seurat.list[[i]], npcs = 50)
  pct<-seurat.list[[i]][['pca']]@stdev / sum(seurat.list[[i]][['pca']]@stdev)*100
  cumu<-cumsum(pct)
  pc.use<-min(which(cumu>90&pct<5)[1],sort(which((pct[1:length(pct)-1]-pct[2:length(pct)])>0.1),decreasing=T)[1]+1)
  rm(pct,cumu)
  #seurat.list[[i]] = RunTSNE(seurat.list[[i]], npcs = 20)
  seurat.list[[i]] = RunUMAP(seurat.list[[i]], dims = 1:pc.use)
  seurat.list[[i]] <- FindNeighbors(seurat.list[[i]], reduction = "pca",dims = 1:pc.use)
  seurat.list[[i]] <- FindClusters(seurat.list[[i]],resolution = 1)
}
#根据地一个循环的结果，去除双细胞
for (i in 1:length(seurat.list)) {
  print(i)
  pct<-seurat.list[[i]][['pca']]@stdev / sum(seurat.list[[i]][['pca']]@stdev)*100
  cumu<-cumsum(pct)
  pc.use<-min(which(cumu>90&pct<5)[1],sort(which((pct[1:length(pct)-1]-pct[2:length(pct)])>0.1),
                                           decreasing=T)[1]+1)
  rm(pct,cumu)
  sweep.res.list <- paramSweep(seurat.list[[i]], PCs = 1:pc.use, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  #ggplot(bcmvn, aes(pK, BCmetric, group = 1)) + geom_point() +geom_line()
  bcmvn <- find.pK(sweep.stats)
  pK <- bcmvn %>%
    filter(BCmetric == max(BCmetric)) %>%
    select(pK)
  pK <- as.numeric(as.character(pK[[1]]))

  homotypic.prop <- modelHomotypic(seurat.list[[i]]$seurat_clusters) ##同源双细胞
  nExp_poi <- round(ncol(seurat.list[[i]]) * 0.04)  # expect 4% doublets;估计有4%的双细胞 依据经验设置 0.075
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) # 减去同源

  seurat.list[[i]] <- doubletFinder(seurat.list[[i]], pN = 0.25, pK = pK, nExp = nExp_poi, PCs = 1:pc.use)
  col = ncol(seurat.list[[i]]@meta.data)
  colnames(seurat.list[[i]]@meta.data)[col] = "DF"
  colnames(seurat.list[[i]]@meta.data)[col-1] = "pANN"

  seurat.list[[i]] <- doubletFinder(seurat.list[[i]], pN = 0.25, pK = pK, nExp = nExp_poi.adj, PCs = 1:pc.use,reuse.pANN ="pANN")
  col = ncol(seurat.list[[i]]@meta.data)
  colnames(seurat.list[[i]]@meta.data)[col] = "DF_adj"
}
seurat = merge(x = seurat.list[[1]], y = seurat.list[-1],
               merge.data = TRUE)
seurat = JoinLayers(seurat)
seurat <- CreateSeuratObject(counts = GetAssayData(seurat,layer = "counts"),meta.data = seurat@meta.data)
seurat <- NormalizeData(seurat, verbose = FALSE)
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
seurat <- ScaleData(seurat, vars.to.regress = c("nFeature_RNA", "pMT"))
seurat <- RunPCA(seurat, npcs = 50)
pct<-seurat[['pca']]@stdev / sum(seurat[['pca']]@stdev)*100
cumu<-cumsum(pct)
pc.use<-min(which(cumu>90&pct<5)[1],sort(which((pct[1:length(pct)-1]-pct[2:length(pct)])>0.1),decreasing=T)[1]+1)
rm(pct,cumu)
#seurat = RunTSNE(seurat, npcs = 20)
seurat = RunUMAP(seurat, dims = 1:pc.use)

#DF 像初步筛查：把“所有长得像双细胞的人”都抓出来（可能误抓）。
#DF_adj 像二次审查：放回那些“长得像但其实是好人”的细胞，只留真正的双细胞。
table(seurat$DF)
table(seurat$DF_adj)

cols<-c('Doublet'='#ea5455','Singlet'='#3f72af')
reduction='umap'
plot_list<-list()
for(g in c(sample.ident,'DF','DF_adj')){
  if(g==sample.ident){cols_touse=sample.cols}else{
    cols_touse=cols}
  plot_list[[g]]<-data.frame(
    x=seurat@reductions[[reduction]]@cell.embeddings[,1],
    y=seurat@reductions[[reduction]]@cell.embeddings[,2],
    category=seurat@meta.data[,g])%>%
    ggplot(aes(x =x,y = y))+
    geom_point(aes(color=category),size=0.2)+
    scale_color_manual(values = cols_touse)+
    labs(x=paste0(reduction,'_1'),y=paste0(reduction,'_2'),color='')+
    guides(color=guide_legend(override.aes = list(size = 3)))+ ggtitle(g)+
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_blank(),
          plot.background = element_rect(fill = "transparent", colour = NA),
          panel.background = element_rect(fill = "transparent", colour = NA))+
    theme_noaxis(axis.line.x.bottom = element_line2(id = 1, xlength = 0.2,
                                                    arrow =grid::arrow(length = unit(0.15,"inches"),
                                                                       type = "closed")),
                 axis.line.y.left = element_line2(id = 2,ylength =0.2,
                                                  arrow = grid::arrow(length = unit(0.15,"inches"),
                                                                      type = "closed")),
                 axis.title = element_text(hjust = 0.01))
}
plot_list%>%patchwork::wrap_plots(ncol=3)
ggsave(filename="doublet_dimplot.pdf",height = 4,width =15)

## 绘制单细胞与双细胞表达谱，双细胞的基因数明显高于单细胞
FeatureStatPlot( srt = scCustomize::Convert_Assay(seurat_object=seurat,convert_to="V3"),
                 group.by="DF", stat.by="nFeature_RNA" ,palcolor=cols,bg_palcolor=cols, #sig_label ='p.format',
                 add_box = TRUE, bg.by="DF")
ggsave(filename="doublet_vlnplot.pdf",height = 4,width = 4)
#过滤doublet
if(T){
  seurat=seurat[, seurat@meta.data[, "DF"] == "Singlet"]##只保留Singlet
  cowplot::plot_grid(ncol = 2, DimPlot(seurat, group.by = "orig.ident") + NoAxes(),
                                     DimPlot(seurat, group.by = "DF") + NoAxes())
  ggsave(filename="Singlet_dimplot.pdf",height = 6,width = 8)
}
seurat <- CreateSeuratObject(counts = GetAssayData(seurat,layer = "counts"),meta.data = seurat@meta.data)

#scDoubletFinder------
#BiocManager::install("scDblFinder")
#library(Seurat)
#library(scDblFinder)
#library(BiocParallel)

#sce<- scDblFinder(as.SingleCellExperiment(seurat),
                  samples="sample", BPPARAM=MulticoreParam(3))
#sce$scDblFinder.score
#sce$scDblFinder.class%>%table()




#四、添加细胞周期分数------
#单独赋值给另一个变量用于计算细胞周期分数，并加入到原对象中
seurat<-seurat %>%
  JoinLayers() %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()
seurat[["RNA"]] <- as(object = seurat[["RNA"]], Class = "Assay")
s.genes <- cc.genes$s.genes#Seurat内置了基因集，也可以自己设定
g2m.genes <- cc.genes$g2m.genes

seurat<- CellCycleScoring(seurat,
                          s.features = s.genes, #指定s期基因集
                          g2m.features = g2m.genes, #指定g2m基因集
                          set.ident = TRUE)#直接帮你设置好active.ident为细胞周期分数
identical(colnames(seurat),colnames(seurat))
seurat$CC.Difference <- seurat$S.Score - seurat$G2M.Score#仅去除G2M和S之间的差别时用，造血系统适用


DimPlot(seurat,group.by = 'Phase',cols=c('G1'='#89c3eb','S'='#dccb18','G2M'='#e45e32'))+
  guides(color=guide_legend(override.aes = list(size = 3)))+
  theme_noaxis(axis.line.x.bottom = element_line2(id = 1, xlength = 0.2,
                                                  arrow =grid::arrow(length = unit(0.15,"inches"),
                                                                     type = "closed")),
               axis.line.y.left = element_line2(id = 2,ylength =0.2,
                                                arrow = grid::arrow(length = unit(0.15,"inches"),
                                                                    type = "closed")),
               axis.title = element_text(hjust = 0.01))
ggsave('cellcycling_phase.pdf',width=5,height=4)

FeaturePlot(seurat,features=c('S.Score','G2M.Score','CC.Difference'),cols=c('#fbfaf5','#2a83a2'),ncol=3)&
  theme_noaxis(axis.line.x.bottom = element_line2(id = 1, xlength = 0.2,
                                                  arrow =grid::arrow(length = unit(0.15,"inches"),
                                                                     type = "closed")),
               axis.line.y.left = element_line2(id = 2,ylength =0.2,
                                                arrow = grid::arrow(length = unit(0.15,"inches"),
                                                                    type = "closed")),
               axis.title = element_text(hjust = 0.01))

ggsave('cellcycling_score.pdf',width=15,height=4)
saveRDS(seurat, "seurat.all_qc.rds")
