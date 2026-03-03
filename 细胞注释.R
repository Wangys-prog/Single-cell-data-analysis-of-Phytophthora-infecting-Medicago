
####加载包----
library(Seurat)
library(tidyverse)
library(harmony)
library(dplyr)
library(patchwork)
library(SeuratWrappers)

scRNA_harmony<- readRDS("E:/苜蓿数据分析/Medicago_data/rds数据/Merge_Medicago_注释.rds")

scRNA_harmony = UpdateSeuratObject(scRNA_harmony)
Idents(scRNA_harmony)<-"cells_type"

levels(Idents(scRNA_harmony))
View(scRNA_harmony)

####读取矩阵数据----
assays <- dir("D:/wangys_uestc_data2/苜蓿数据分析/单细胞矩阵数据/")
dir<- paste0("D:/wangys_uestc_data2/苜蓿数据分析/单细胞矩阵数据/",assays)
samples_name<- c("I_18_1","I_18_2","I_9_1","I_9_2","M_18_1","M_18_2","M_9_1","M_9_2")
scRNAlist<- list()
for(i in 1:length(dir)){
  counts<- Read10X(data.dir = dir[i])
  scRNAlist[[i]]<- CreateSeuratObject(counts,project = samples_name[i],
                                      min.cells=3,min.features=200)
  ##给细胞barcode加入前缀，防止合并后barcode重名
  scRNAlist[[i]]<-RenameCells(scRNAlist[[i]],add.cell.id=sample_name[i])
  ##计算线粒体基因比例
  if(T){
    scRNAlist[[i]][["percent.mt"]]<-PercentageFeatureSet(scRNAlist[[i]],pattern = "^MT-")
  }
  ##计算核糖体基因比例
  if(T){
    scRNAlist[[i]][["percent.rb"]]<-PercentageFeatureSet(scRNAlist[[i]],pattern = "^RP[SL]")
  }
  ##计算红细胞基因比例
  if(T){
    HB.genes<-c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
    HB.genes<- CaseMatch(HB.genes,rownames(scRNAlist[[i]]))
    scRNAlist[[i]][["percent.HB"]]<-PercentageFeatureSet(scRNAlist[[i]],features = HB.genes)
  }
}
##给列表命名并保存
# dir.create("Integrate_Medicago")
# setwd("./Integrate_Medicago")

names(scRNAlist)<- samples_name
##只是一个一个单独的seurat
system.time(saveRDS(scRNAlist,file="D:/wangys_uestc_data2/苜蓿数据分析/Medicago_data/rds数据/Integrate_Medicago.rds"))

####merge将scRNAlist合成一个Seurat的对象----
scRNA<-merge(scRNAlist[[1]],scRNAlist[2:length(scRNAlist)])
scRNA
table(scRNA$orig.ident)
saveRDS(scRNA,file="D:/wangys_uestc_data2/苜蓿数据分析/Medicago_data/rds数据/Merge_Medicago.rds")

####加载rds数据----
Medicago_scRNA_object<- readRDS("D:/wangys_uestc_data2/苜蓿数据分析/Medicago_data/rds数据/Merge_Medicago.rds")
Medicago_scRNA_object
View(Medicago_scRNA_object)

####convert a v5 assay to a v3 assay
Medicago_scRNA_object[["RNA4"]] <- as(object = Medicago_scRNA_object[["RNA"]], Class = "Assay")

####convert a v3 assay to a v5 assay

Medicago_scRNA_object[["RNA5"]] <- as(object = Medicago_scRNA_object[["RNA3"]], Class = "Assay5")


####数据质控----
setwd("D:/wangys_uestc_data2/苜蓿数据分析/")
dir.create("1.SingleCell_QC")
devtools::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)
##去除双细胞
obj = SplitObject(Medicago_scRNA_object, split.by = "orig.ident")
obj_rm=list() #创建空list用于存放去除双细胞以后的seurat对象
doublets_plot = list() #创建空list用于存放双细胞的分布情况（图片信息）
pc.num = 1:20 #设置的经验PC维数

#RemoveDoublets函数
#设置对应的double.rate，pN值，pc的维数
RemoveDoublets <-function(
    object,
    doublet.rate,
    pN=0.25,
    pc.num=1:30
){
  ## 寻找最优pK值
  sweep.res.list <- paramSweep(object, PCs = pc.num, sct = F)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
  bcmvn <- find.pK(sweep.stats) #求出最大bcmvn值
  pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric() #最大bcmvn值对应的pk值为最优
  ## 排除不能检出的同源doublets，优化期望的doublets数量
  homotypic.prop <- modelHomotypic(object$seurat_clusters) 
  nExp_poi <- round(doublet.rate*ncol(object)) 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  seu.scored <- doubletFinder(object, PCs = pc.num, pN = 0.25, pK = pK_bcmvn, 
                              nExp = nExp_poi.adj, reuse.pANN = F, sct = F)
  # 选出储存双细胞预测结果的列
  cname <-colnames(seu.scored[[]])
  DF<-cname[grep('^DF',cname)]
  seu.scored[["doublet"]] <- as.numeric(seu.scored[[DF]]=="Doublet")# 去除双细胞
  seu.removed <- subset(seu.scored, subset = doublet != 1)
  p1 <- DimPlot(seu.scored, group.by = DF)
  res.list <- list("plot"=p1, "obj"=seu.removed)
  return(res.list) #返回对应的信息(plot是对应的图片，obj是对应的seurat的list)
}
#循环进行双细胞的去除
for( i in names(obj)){
  print(i)
  #归一化标准化
  obj[[i]] <- NormalizeData(obj[[i]])
  obj[[i]] <- FindVariableFeatures(obj[[i]], selection.method = "vst", nfeatures = 2000)
  obj[[i]] <- ScaleData(obj[[i]])
  
  #降维聚类
  obj[[i]] <- RunPCA(obj[[i]])
  obj[[i]] <- RunUMAP(obj[[i]], dims = 1:20)
  obj[[i]] <- FindNeighbors(obj[[i]], dims = pc.num) %>% FindClusters(resolution = 0.3)#此处选择了一个较小的resolution，避免over-clustering的情况
  #去除双细胞***
  tmp <- RemoveDoublets(obj[[i]], doublet.rate=0.008,pc.num=pc.num)# 1000细胞对应的doublets rate是0.8%， 如果细胞数大于1000，可以根据每1000个细胞增加0.8%的方式进行估算和替换
  obj_rm[[i]] <- tmp$obj #seurat对象信息 
  doublets_plot[[i]] <- tmp$plot #图片信息
}


####质量控制----
######质量控制之前----
scRNA<-merge(obj_rm[[1]],obj_rm[2:length(obj_rm)])
#合并去除完双细胞的seurat对象
Idents(scRNA) <- 'orig.ident' #更改Idents信息
table(Idents(scRNA))
#percent.mt计算***
scRNA[["percent.mt"]] = PercentageFeatureSet(scRNA,pattern = "^MT-") #因为人的线粒体基因是MT-开头的，所以可以通过这种方式匹配，如果是小鼠的话，是以“mt-”开头的。其他特殊物种可以直接提供线粒体基因
#scRNA[["percent.mt"]]等于scRNA@meta.data$percent.mt

#percent.HB计算***
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ") #直接提供血红细胞基因
HB_m <- match(HB.genes, rownames(scRNA@assays$RNA)) #血红细胞基因和seurat对象的基因匹配，返回对应的基因定位信息

HB.genes <- rownames(scRNA@assays$RNA)[HB_m] #将位置信息定位回去
HB.genes <- HB.genes[!is.na(HB.genes)] #去除NA信息
scRNA[["percent.HB"]]<-PercentageFeatureSet(scRNA, features=HB.genes) #根据提供的红细胞基因列表计算红细胞基因表达的比例

#数据过滤前的可视化***
beforeQC_vlnplot = VlnPlot(scRNA, 
                           features = c("nFeature_RNA", 
                                        "nCount_RNA", 
                                        "percent.mt",
                                        "percent.HB"), 
                           ncol = 4, 
                           pt.size = 0) #设置对应的图片排版和点的大小
beforeQC_vlnplot
#图片保存
ggsave("D:/wangys_uestc_data2/苜蓿数据分析/1.SingleCell_QC/BeforeQC_nFeature_nCount_percent.mt_percent.HB_vlnplot.pdf", plot = beforeQC_vlnplot)
ggsave("D:/wangys_uestc_data2/苜蓿数据分析/1.SingleCell_QC/BeforeQC_nFeature_nCount_percent.mt_percent.HB_vlnplot.png", plot = beforeQC_vlnplot)

#统计分析
summary(scRNA@meta.data[,c("nFeature_RNA","nCount_RNA", "percent.mt", "percent.HB")]) #最大最小值，平均数中位数等

theme.set2=theme(axis.title.x=element_blank())
plot.features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb","percent.HB")
group="orig.ident"
i=1
plots=list()
for(i in 1:length(plot.features)){
  feature <- plot.features[i]
  plots[[i]]=VlnPlot(scRNA,group.by=group,pt.size=0,
                     features=feature)+theme.set2+NoLegend()
}

violin <- wrap_plots(plots=plots,nrow=2)
ggsave("D:/wangys_uestc_data2/苜蓿数据分析/1.SingleCell_QC/BeforeQC_nFeature_nCount_percent.mt_percent.HB_vlnplot.jpg",
       plot = violin,width=9,height = 8)


######质量控制----
minGene=200
#maxGene=4000
#maxUMI=15000
pctMT=10
pctHB=1

scRNA = subset(scRNA,subset = nFeature_RNA > minGene & percent.mt < pctMT &percent.HB< pctHB)

plot.features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb","percent.HB")
group="orig.ident"
plots=list()
for(i in seq_along(plot.features)){
  plots[[i]]=VlnPlot(scRNA,group.by=group,pt.size=0,
                     features=plot.features[i])+theme.set2+NoLegend()
}
violin <- wrap_plots(plots=plots,nrow=2)
ggsave("D:/wangys_uestc_data2/苜蓿数据分析/1.SingleCell_QC/AftherQC_nFeature_nCount_percent.mt_percent.HB_vlnplot.jpg",
       plot = violin,width=9,height = 8)

#统计分析
summary(scRNA@meta.data[,c("nFeature_RNA","nCount_RNA", "percent.mt", "percent.HB")])

saveRDS(scRNA,"D:/wangys_uestc_data2/苜蓿数据分析/Medicago_data/rds数据/Merge_Medicago_After_QC.rds")
scRNA<-readRDS("D:/wangys_uestc_data2/苜蓿数据分析/Medicago_data/rds数据/Merge_Medicago_After_QC.rds")

####归一化和标准化----
#归一化与标准化***
scRNA <- NormalizeData(scRNA) #归一化
scRNA <- FindVariableFeatures(scRNA,  selection.method = "vst",nfeatures = 2000) #选择高变基因，特征选择找出包含信息最多的基因，通常选择2000个高变基因
scRNA <- ScaleData(scRNA, features = VariableFeatures(scRNA)) #标准化

####计算每个样本平均UMI和gene----

scRNA_harmony<- JoinLayers(scRNA_harmony)
saveRDS(scRNA_harmony,"D:/wangys_uestc_data2/苜蓿数据分析/Medicago_data/rds数据/Merge_Medicago_注释.rds")
GetAssayData(scRNA_harmony, slot = "counts")
scRNA_harmony[["UMI_counts"]] <- colSums(GetAssayData(scRNA_harmony, assay = "RNA", slot = "counts"))

# 提取样本信息（假设样本信息存储在meta.data中的"sample"列中）
sample_info <- scRNA_harmony@meta.data$orig.ident

# 将UMI数按样本分组
umi_counts_per_sample <- tapply( scRNA_harmony$UMI_counts, sample_info, mean)

# 打印每个样本的平均UMI数
print(umi_counts_per_sample)

####降维聚类和批次矫正----

#创建文件夹
dir.create("D:/wangys_uestc_data2/苜蓿数据分析/2.Clustering")

#PCA分析***
scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA)) #使用高变基因进行PCA降维（原始表达矩阵太大）
scRNA <- RunUMAP(scRNA, reduction = "pca", dims = 1:20)
#散点图
plot1 <- DimPlot(scRNA, reduction = "pca", group.by="orig.ident") #每个样本在pca降维结果中的分布情况
head(scRNA@reductions$pca@cell.embeddings)

#碎石图
plot2 <- ElbowPlot(scRNA, ndims=50, reduction="pca") #选取前50维绘制elbowplot，纵坐标表示标准差，代表这个维度的贡献度

#保存
plot3 <- plot1+plot2
plot3
ggsave("D:/wangys_uestc_data2/苜蓿数据分析/2.Clustering/pca.pdf", plot = plot3, width = 8, height = 4)
ggsave("D:/wangys_uestc_data2/苜蓿数据分析/2.Clustering/pca.png", plot = plot3, width = 8, height = 4)

####聚类分析----
#KNN + SNN***
scRNA <- FindNeighbors(scRNA, dims = 1:20) #取前20个维度

#Louvain***
for (res in c(0.01, 0.05, 0.1, 0.2 ,0.3, 0.5, 0.6,0.8, 1)) {
  scRNA = FindClusters(scRNA,
                       resolution = res, algorithm = 1)
  
}

#设置resolution，resolution越大，分的亚群越多，一般免疫细胞T和B细胞的resolution设定为0.6或者0.8，其他细胞在0.4左右

#保存聚类信息
metadata <- scRNA@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters) #将聚类结果另存为data.frame
write.csv(cell_cluster,'D:/wangys_uestc_data2/苜蓿数据分析/2.Clustering/cell_cluster.csv',row.names = F, quote = F)

######非线性降维----

#UMAP***
scRNA <- RunUMAP(scRNA, dims = 1:20) #UMAP二次降维，取前20个维度

#散点图
plot3 = DimPlot(scRNA, reduction = "umap", label=T) #查看clusters在UMAP降维图中的分布
plot3
ggsave("D:/wangys_uestc_data2/苜蓿数据分析/2.Clustering/UMAP.pdf",plot3,width = 8,height = 6)
ggsave("D:/wangys_uestc_data2/苜蓿数据分析/2.Clustering/UMAP.png",plot3,width = 8,height = 6)

plot4 = DimPlot(scRNA, reduction = "umap", group.by='orig.ident') #查看每个样本在UMAP降维图中的分布
plot5 = DimPlot(scRNA, reduction = "umap", split.by='orig.ident') #查看每个样本在UMAP降维图中的分面图
plotc <- plot4+plot5 
plotc
ggsave("2.Clustering/UMAP_cluster_sample.pdf", plot = plotc, width = 10, height = 5)
ggsave("2.Clustering/UMAP_cluster_sample.png", plot = plotc, width = 10, height = 5)

#保存rds
saveRDS(scRNA,"D:/wangys_uestc_data2/苜蓿数据分析/Medicago_data/rds数据/Merge_Medicago_降维聚类.rds")
scRNA <-readRDS("D:/wangys_uestc_data2/苜蓿数据分析/Medicago_data/rds数据/Merge_Medicago_降维聚类.rds")

####批次矫正----
#mnn不需要PCA的结果，相当于替代了PCA，但是要进行后续的聚类和非线性降维（UMAP）

#按样本分析
scRNAlist <- SplitObject(scRNA, split.by = "orig.ident") #将合并的seurat对象重新拆分
scRNAlist <- lapply(scRNAlist, FUN = function(x) NormalizeData(x)) #对每一个样本进行归一化
scRNAlist <- lapply(scRNAlist, FUN = function(x) FindVariableFeatures(x)) #对每一个样本寻找高变基因

#fastMNN去批次***
scRNA_mnn <- RunFastMNN(object.list = scRNAlist) 

#聚类和UMAP降维
scRNA_mnn <- FindVariableFeatures(scRNA_mnn) #对去完批次以及合并后的对象重新进行高变基因的选取
scRNA_mnn <- RunUMAP(scRNA_mnn, reduction = "mnn", dims = 1:20,reduction.name = "umap.mnn") #UMAP降维
scRNA_mnn<- RunTSNE(object = scRNA_mnn,reduction = "mnn", reduction.key = "tsneharmony_",dims.use = 1:20,reduction.name = "tsne.harmony")
scRNA_mnn <- FindNeighbors(scRNA_mnn, reduction = "mnn", dims = 1:20) #SNN + KNN
scRNA_mnn <- FindClusters(scRNA_mnn) #Louvain
View(scRNA_mnn)
#散点图
p1 <- DimPlot(scRNA_mnn, reduction = "mnn", group.by = "orig.ident", pt.size=0.1) + 
  ggtitle("Integrated by MNN")#去批次后的可视化
p2 <- DimPlot(scRNA, group.by="orig.ident", reduction = "pca",pt.size=0.1) + 
  ggtitle("No integrated")#去批次前的可视化
p = p1 + p2 + plot_layout(guides='collect')
#合并两者结果，将legend统一放到右边
#比较去批次前后聚类结果的变化
ggsave('D:/wangys_uestc_data2/苜蓿数据分析/2.Clustering/MNN_PCA_sample.pdf', p, width=8, height=4)
ggsave('D:/wangys_uestc_data2/苜蓿数据分析/2.Clustering/MNN_PCA_sample.png', p, width=8, height=4)


p3 <- DimPlot(scRNA_mnn, reduction = "umap.mnn", pt.size=0.1,label = T) #去批次后重降维聚类结果的可视化
p4 <- DimPlot(scRNA_mnn, reduction = "umap.mnn", split.by="orig.ident", pt.size=0.1) #去批次后重降维聚类结果的分样本可视化展示
p = p3 + p4 + plot_layout(guides='collect')  #比较各个样本内cluster的组成情况
ggsave('D:/wangys_uestc_data2/苜蓿数据分析/2.Clustering/2.Clustering/MNN.pdf', p, width=8, height=4)
ggsave('D:/wangys_uestc_data2/苜蓿数据分析/2.Clustering/2.Clustering/MNN.png', p, width=8, height=4)

#保存rds
View(scRNA_mnn)

saveRDS(scRNA_mnn,"D:/wangys_uestc_data2/苜蓿数据分析/Medicago_data/rds数据/Merge_Medicago_mnn去批次.rds")
scRNA<- readRDS("D:/wangys_uestc_data2/苜蓿数据分析/Medicago_data/rds数据/Merge_Medicago_mnn去批次.rds")

####harmony去批次----
#harmony分析***
scRNA_harmony <- RunHarmony(scRNA, group.by.vars = "orig.ident")#harmony需要pca降维的结果，而scRNA对象此前已经做过PCA降维了，所以此处不再重复

#进行UMAP降维和聚类
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:20,reduction.name = "umap.harmony")#UMAP降维
scRNA_harmony <- RunTSNE(scRNA_harmony,reduction = "harmony", reduction.key = "tsneharmony_",dims.use = 1:20,reduction.name = "tsne.harmony")
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:20)
scRNA_harmony <- FindClusters(scRNA_harmony)#标准聚类
p3 = DimPlot(scRNA_harmony, reduction = "tsne.harmony", group.by ='orig.ident')+ggtitle("tsne.harmony")
p4 = DimPlot(scRNA_harmony, reduction = "umap.harmony", group.by ='orig.ident')+ggtitle("umap.harmony")
p5 = DimPlot(scRNA, reduction = "umap", group.by ='orig.ident')+ggtitle("umap")
p6 = DimPlot(scRNA_mnn, reduction = "umap.mnn", group.by ='orig.ident')+ggtitle("mnn")
p7 = DimPlot(scRNA_mnn, reduction = "tsne.harmony", group.by ='orig.ident')+ggtitle("tsne.mnn")
plotc <- p5+p3+p4+p6+p7+plot_layout(guides='collect')
ggsave("D:/wangys_uestc_data2/苜蓿数据分析/2.Clustering/harmony_mnn.pdf", plot = plotc, width = 10, height = 5)
ggsave("D:/wangys_uestc_data2/苜蓿数据分析/2.Clustering/harmony_mnn.png", plot = plotc, width = 10, height = 5)

#保存RDS
saveRDS(scRNA_harmony,"D:/wangys_uestc_data2/苜蓿数据分析/Medicago_data/rds数据/Merge_Medicago_harmony去批次.rds")

####去批次之后细胞注释----
scRNA_harmony = readRDS("D:/wangys_uestc_data2/苜蓿数据分析/Medicago_data/rds数据/Merge_Medicago_harmony去批次.rds")

for (res in c(0.01, 0.05, 0.1, 0.2 ,0.4,0.3, 0.5, 0.6,0.8, 1)) {
  scRNA_harmony = FindClusters(scRNA_harmony,
                               resolution = res, algorithm = 1)
  
}
scRNA_harmony = FindClusters(scRNA_harmony,
                             resolution = 0.4, algorithm = 1)

scRNA_harmony <- FindVariableFeatures(scRNA_harmony,selection.method = "vst",nfeatures = 3000)
View(scRNA_harmony)

saveRDS(scRNA_harmony,"D:/wangys_uestc_data2/苜蓿数据分析/Medicago_data/rds数据/Merge_Medicago_harmony去批次.rds")

gene_data <- VariableFeatures(object = scRNA_harmony)


Idents(scRNA_harmony)<-"RNA_snn_res.0.4"
View(scRNA_harmony)
Idents(scRNA_harmony)<- "RNA"

levels(Idents(scRNA_harmony))

DimPlot(scRNA_harmony, reduction = "umap.harmony", group.by ='RNA_snn_res.0.4')+ggtitle("umap.harmony")

####chooseR 选择res----

library(clustree)
library(patchwork)
# seq<-seq(0.1,1,by=0.1)
# for (res in seq){
#   scRNA_harmony<- FindClusters(scRNA_harmony,resolution = res)
# }

p1<-clustree(scRNA_harmony,prefix="RNA_snn_res.")+coord_flip()
p2<-DimPlot(scRNA_harmony,group.by = "RNA_snn_res.0.3",label = T)
p1+p2+plot_layout(widths = c(3,1))

Idents(scRNA_harmony)<-"RNA_snn_res.0.4"
View(scRNA_harmony)
Idents(scRNA_harmony)<- "RNA"

levels(Idents(scRNA_harmony))

DimPlot(scRNA_harmony, reduction = "umap.harmony", group.by ='RNA_snn_res.0.4',label = T)+ggtitle("umap.harmony")

DimPlot(scRNA_harmony, reduction = "tsne.harmony", group.by ='RNA_snn_res.0.4')+ggtitle("tsne.harmony")

####小提琴图marker----

cell_marker<- read.csv("D:/wangys_uestc_data2/苜蓿数据分析/3.人工注释/marker_gene.csv")
######删除找不到的基因
nofinder<- c("MTR-6g029180", "MTR-5g083030", "MTR-8g06995",
             "MTR-4g059310", "MTR-1291s0010", "MTR-4g059840", 
             "MTR-8g087470", "MTR-8g027290", "MTR-1g093030", "MTR-5g080050",
             "MTR-7g027190", "MTR-1g044650", "MTR-1g008630",
             "MTR-1g113750", "MTR-3g101640", "MTR-8g094550", "MTR-2g101740",
             "MTR-1g087935", "MTR-3g111780", "MTR-2g093810", "MTR-3g097190", 
             "MTR-4g011850", "MTR-3g097370", "MTR-3g097390", "MTR-3g060400", 
             "MTR-3g060410", "MTR-3g060770", "MTR-3g060820","MTR-5g011900",
             "MTR-3g460990", "MTR-8g008690" ,"MTR-6g477820", "MTR-0198s0030", 
             "MTR-7g427710", "MTR-2g030600",
             "MTR-2g030660", "MTR-2g083350", "MTR-2g083330", "MTR-2g030690", 
             "MTR-7g027265", "MTR-7g027605","MTR-5g046040", "MTR-0653s0010" ,
             "MTR-1g103450", "MTR-1g103470", "MTR-6g084630"," MTR-1g097720","MTR-1g097720"
)

cell_marker_filtered <- cell_marker[!cell_marker$gene_name %in% nofinder, ]

markers <-cell_marker_filtered$gene_name
markers <- as.data.frame(markers)
Idents(scRNA_harmony)<-"RNA_snn_res.0.4"
unique(cell_marker_filtered$gene_name)
library(randomcoloR)

my36colors = distinctColorPalette(23)


my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175'
)


####ggplot2画小提琴图----
View(scRNA_harmony)
vln.dat=FetchData(scRNA_harmony,c(unique(cell_marker_filtered$gene_name),"RNA_snn_res.0.4"))
colnames(vln.dat)
#宽转长
vln.dat.melt <- reshape2::melt(vln.dat, id.vars = c("RNA_snn_res.0.4"), 
                               measure.vars = unique(cell_marker_filtered$gene_name),
                               variable.name = "gene", 
                               value.name = "Expr") 

vln.dat.melt$gene<- factor(vln.dat.melt$gene,levels = unique(vln.dat.melt$gene))

View(vln.dat.melt)

library(cowplot)
vln.dat.melt$RNA_snn_res.0.4<- factor(vln.dat.melt$RNA_snn_res.0.4,levels = c("0","1","2","3","4","5","6",
                                                                              "7","8","9","10","11","12","13","14","15"))
vln.dat.melt$cells_type <- factor(vln.dat.melt$cells_type,levels=c("Lateral root cap", "Cortex","Endodermis","Pericycle","Xylem",
                                                                   "Atrichoblast","Root hairs","Unknow" ))


markers_name<- read.csv("D:/wangys_uestc_data2/苜蓿数据分析/3.人工注释/marker_gene.csv",header = T)
markers_name$Our_Feature_ID
names(markers_name)
vln.dat.melt$feature_id = markers_name[match(vln.dat.melt$gene,markers_name$Our_Feature_ID),3]
vln.dat.melt$feature_id<- factor(vln.dat.melt$feature_id,levels = unique(vln.dat.melt$feature_id))


p1 <- ggplot(vln.dat.melt, aes(feature_id, Expr, fill = feature_id)) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE) +
  scale_y_continuous(expand = c(0, 0), position="right", labels = function(x)
    c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
  facet_grid(rows = vars(RNA_snn_res.0.4), scales = "free", switch = "y") +
  scale_fill_manual(values = my36colors) + 
  theme_cowplot(font_size = 20) +
  theme(legend.position = "none", panel.spacing = unit(0, "lines"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        plot.margin = margin(7, 7, 0, 7, "pt"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.y.left = element_text(angle = 0),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black")
  ) +ylab("Expression Level")
#ggtitle("Feature on x-axis with annotation") + ylab("Expression Level")
p1

ggsave("D:/wangys_uestc_data2/苜蓿数据分析/3.人工注释/细胞marker小提琴图RNA_snn_res.0.4.jpg",p1,width = 7,height = 7)

ggsave("D:/wangys_uestc_data2/苜蓿数据分析/3.人工注释/细胞marker小提琴图RNA_snn_res.0.4.pdf",p1,width = 7,height = 7)
####命名cluster----

old_to_new_groups <- c("0" = "Pericycle",
                       "1" = "Lateral root cap",
                       "10" = "Cortex",
                       "11" = "Endodermis",
                       "12" = "Pericycle",
                       "13" ="Pericycle",
                       "14" ="Atrichoblast",
                       "15"="Xylem",
                       "2"="Unknow",
                       "3"= "Cortex",
                       "4"="Root hairs",
                       "5"="Unknow",
                       "6"="Unknow",
                       "7"="Cortex",
                       "8" ="Cortex",
                       "9"="Pericycle")


scRNA_harmony@meta.data$cells_type <- plyr::mapvalues(scRNA_harmony@meta.data$RNA_snn_res.0.4,
                                                      from = names(old_to_new_groups), 
                                                      to = old_to_new_groups)

table(scRNA_harmony@meta.data$cells_type)

DimPlot(scRNA_harmony, group.by = "cells_type",reduction = "umap.harmony", pt.size=0.1)

saveRDS(scRNA_harmony,"D:/wangys_uestc_data2/苜蓿数据分析/Medicago_data/rds数据/Merge_Medicago_注释.rds")


####查找每个cluster的差异基因----

scRNA_harmony <- JoinLayers(scRNA_harmony, assay = "RNA",layers="RNA")
scRNA_harmony<-JoinLayers(scRNA_harmony)
all.markers = FindAllMarkers(scRNA_harmony,assay = "RNA",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)#此处的n可以替换为任意数字，如果只想看top1，就把10替换成1 #n=10是降序排列，如果需要升序排列，则是n=-10

#保存 Marker 基因
write.table(all.markers, 
            "3.Marker/all_Markers_of_each_clusters.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t")
write.table(top10, 
            "3.Marker/top10_Markers_of_each_clusters.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t")
saveRDS(scRNA_harmony,"D:/wangys_uestc_data2/苜蓿数据分析/Medicago_data/rds数据/Merge_Medicago_注释.rds")

####添加样本信息到metadata----
scRNA_harmony <-readRDS("D:/wangys_uestc_data2/苜蓿数据分析/Medicago_data/rds数据/Merge_Medicago_注释.rds")
library(tidyverse)
library(dplyr)
metadata_res<-data.table::fread("D:/wangys_uestc_data2/苜蓿数据分析/4.metadata信息/group.csv",header = TRUE)  
metadata<- FetchData(scRNA_harmony,"orig.ident")
metadata$cell_id <- rownames(metadata)
metadata<- left_join(x=metadata,y=metadata_res,by="orig.ident")
rownames(metadata)<-metadata$cell_id
scRNA_harmony <- AddMetaData(scRNA_harmony,metadata = metadata)
table(scRNA_harmony@meta.data$group1)
table(scRNA_harmony@meta.data$group2)
table(scRNA_harmony@meta.data$group3)

saveRDS(scRNA_harmony,"D:/wangys_uestc_data2/苜蓿数据分析/Medicago_data/rds数据/Merge_Medicago_注释.rds")

###细胞注释UMap 图 提取降维信息用ggplot 画图####
scRNA_harmony <-readRDS("G:/wangys_uestc_data2/苜蓿数据分析/Medicago_data/rds数据/Merge_Medicago_注释.rds")
library(tidyverse)
umap.harmonyredu<- scRNA_harmony@reductions$umap.harmony@cell.embeddings%>%
  as.data.frame() %>% 
  cbind(cell_type = scRNA_harmony@meta.data$cells_type)%>% 
  cbind(group1= scRNA_harmony@meta.data$group1) %>% 
  cbind(group2 = scRNA_harmony@meta.data$group2)%>%
  cbind(group3 = scRNA_harmony@meta.data$group3)%>%
  cbind(cluster = scRNA_harmony@meta.data$RNA_snn_res.0.4)


names(umap.harmonyredu)

###颜色设置###画小坐标轴####https://cloud.tencent.com/developer/article/1924260

cell_type_color <- c("#FF6347","#ffc556","#1E90FF","#9370DB","#0000FF","#00FFFF","#FF00ff","#20B2AA")

umap.harmonyredu$cell_type<- factor(umap.harmonyredu$cell_type, levels = c("Lateral root cap", "Cortex","Endodermis","Pericycle","Xylem",
                                                                           "Atrichoblast","Root hairs","Unknow" ))

write.csv(umap.harmonyredu,"G:/苜蓿数据分析/Medt_scRNA_manuscript/原始数据上传/画图数据/umap.harmonyredu.csv")

cell_type_plot <- ggplot(umap.harmonyredu,aes(x= umapharmony_1, y = umapharmony_2,color = cell_type)) +  
  geom_point(size = 0.5, alpha =1)  +  
  scale_color_manual(values = cell_type_color)+
  xlab("UMAP1")+
  ylab("UMAP2")+
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white"))+
  theme(legend.title = element_blank(), #去掉legend.title 
        legend.key=element_rect(fill='white'), #
        legend.text = element_text(size=15), #设置legend标签的大小
        legend.key.size=unit(1,'cm') ) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes=list(shape = 15,size=5)))+#设置legend中 点的大小
  geom_segment(aes(x = min(umap.harmonyredu$umapharmony_1) , y = min(umap.harmonyredu$umapharmony_2) ,
                   xend = min(umap.harmonyredu$umapharmony_1) +3, yend = min(umap.harmonyredu$umapharmony_2) ),
               colour = "black", size=1,arrow = arrow(length = unit(0.5,"cm")))+ 
  geom_segment(aes(x = min(umap.harmonyredu$umapharmony_1)  , y = min(umap.harmonyredu$umapharmony_2)  ,
                   xend = min(umap.harmonyredu$umapharmony_1) , yend = min(umap.harmonyredu$umapharmony_2) + 3),
               colour = "black", size=1,arrow = arrow(length = unit(0.5,"cm"))) +
  annotate("text", x = min(umap.harmonyredu$umapharmony_1) +1.5, y = min(umap.harmonyredu$umapharmony_2) -1, label = "UMAP_1",
           color="black",size = 5, fontface="bold" ) + 
  annotate("text", x = min(umap.harmonyredu$umapharmony_1) -1, y = min(umap.harmonyredu$umapharmony_2) + 1.5, label = "UMAP_2",
           color="black",size = 5, fontface="bold" ,angle=90)

cell_type_plot

cell_type_med <- umap.harmonyredu %>%
  group_by(cell_type) %>%
  summarise(
    umapharmony_1 = median(umapharmony_1),
    umapharmony_2 = median(umapharmony_2)
  )

library(ggrepel)
cell_type_med

cell_type2<-
  cell_type_plot +geom_label_repel(aes(label=cell_type),size=4,color="black",fontface="bold",data = cell_type_med,
                                   point.padding=unit(0.5, "lines"),fill = alpha(c("white"),0.5),
                                   segment.size=0.5, nudge_x=0, nudge_y = 0,direction="y",max.overlaps=50)+
  theme(legend.position = "right")


ggsave('D:/wangys_uestc_data2/苜蓿数据分析/3.人工注释/cell_type细胞注释.jpg',cell_type2,width =10,height = 6)
ggsave('D:/wangys_uestc_data2/苜蓿数据分析/3.人工注释/cell_type细胞注释.pdf',cell_type2,width =10,height = 6)

####umap分组注释----
##调整legend 顺序
cluster_type_color<- c("#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#FFFF00", "#808000","#FF00FF","#B9D3ee","#7B68EE",
                       "#9400D3","#FF1493","#A0522D","#b03060","seagreen4")
unique(umap.harmonyredu$cluster)

umap.harmonyredu$cluster <- factor(umap.harmonyredu$cluster, levels=c("0",'1', '2', '3', '4', 
                                                                      '5',"6","7","8","9","10","11","12","13","14","15"))

cluster_plot <- ggplot(umap.harmonyredu,aes(x= umapharmony_1, y = umapharmony_2,color = cluster)) +  
  geom_point(size = 0.2, alpha =1)  +  
  scale_color_manual(values = cluster_type_color)+
  xlab("UMAP1")+
  ylab("UMAP2")+
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white"))+
  theme(legend.title = element_blank(), #去掉legend.title 
        legend.key=element_rect(fill='white'), #
        legend.text = element_text(size=10), #设置legend标签的大小
        legend.key.size=unit(1,'cm')) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes = list(size=7)))+ #设置legend中 点的大小
  geom_segment(aes(x = min(umap.harmonyredu$umapharmony_1) , y = min(umap.harmonyredu$umapharmony_2) ,
                   xend = min(umap.harmonyredu$umapharmony_1) +3, yend = min(umap.harmonyredu$umapharmony_2) ),
               colour = "black", size=1,arrow = arrow(length = unit(0.5,"cm")))+ 
  geom_segment(aes(x = min(umap.harmonyredu$umapharmony_1), y = min(umap.harmonyredu$umapharmony_2)  ,
                   xend = min(umap.harmonyredu$umapharmony_1), yend = min(umap.harmonyredu$umapharmony_2) + 3),
               colour = "black", size=1,arrow = arrow(length = unit(0.5,"cm"))) +
  annotate("text", x = min(umap.harmonyredu$umapharmony_1) +1.5, y = min(umap.harmonyredu$umapharmony_2) -1, label = "UMAP_1",
           color="black",size = 5, fontface="bold" ) + 
  annotate("text", x = min(umap.harmonyredu$umapharmony_1) -1, y = min(umap.harmonyredu$umapharmony_2) + 1.5, label = "UMAP_2",
           color="black",size = 5, fontface="bold" ,angle=90)


cluster_plot

cell_type_med <- umap.harmonyredu %>%
  group_by(cluster) %>%
  summarise(
    umapharmony_1 = median(umapharmony_1),
    umapharmony_2 = median(umapharmony_2)
  )

library(ggrepel)
cell_type_med

cluster_plot2<-
  cluster_plot +geom_label_repel(aes(label=cluster),size=4,color="black",fontface="bold",data = cell_type_med,
                                 point.padding=unit(0.5, "lines"),fill = alpha(c("white"),0.5),
                                 segment.size=0.5, nudge_x=0, nudge_y = 0,direction="y",max.overlaps=50)+
  theme(legend.position = "none")

ggsave('D:/wangys_uestc_data2/苜蓿数据分析/3.人工注释/cell_type细胞注释cluster.jpg',cluster_plot2,width =8,height = 6)

ggsave('D:/wangys_uestc_data2/苜蓿数据分析/3.人工注释/cell_type细胞注释cluster.pdf',cluster_plot2,width =8,height = 6)
ggsave('D:/wangys_uestc_data2/苜蓿数据分析/3.人工注释/cell_type细胞注释cluster小.pdf',cluster_plot2,width =5,height = 4)



####dotplot----
devtools::install_github("junjunlab/jjAnno",force = TRUE)
library(jjAnno)

markers <- c("MTR-2g080260","MTR-4g059670","MTR-4g035590","MTR-4g059630",#Lateral root cap
             "MTR-4g088195","MTR-4g088160","MTR-1g074930","MTR-4g101330",#Cortex
             "MTR-1g058560","MTR-7g074650","MTR-4g073140",#Endodermis
             "MTR-8g083150","MTR-4g063090","MTR-8g070770",#Pericycle
             "MTR-1g101830","MTR-2g093990","MTR-7g015960","MTR-4g058860","MTR-7g101080",#Xylem
             "MTR-4g118810","MTR-4g047800", #Atrichoblast
             "MTR-1g075380","MTR-5g090250" #Root hairs
)   

markers <- c("Medtr2g080260","Medtr4g059670","MtBRN","MtRC7",#Lateral root cap
             "MtIFS1","MtIFS3","MtPT5","Medtr4g101330",#Cortex
             "MtCASPL-2","MtSCR","MtPTI11",#Endodermis
             "MtPHO1.1","MtPHO1.2","MtPHO1.3",#Pericycle
             "MtPrx13","Medtr2g093990","Medtr7g015960","Medtr4g058860", "Medtr7g101080",#Xylem
             "Medtr4g118810","Medtr4g047800",#Atrichoblast
             "Medtr1g075380", "Medtr5g090250"#Root hairs
) 


cluster<-c("Lateral root cap", "Cortex","Endodermis","Pericycle","Xylem",
           "Atrichoblast","Root hairs","Unknow")
Idents(scRNA_harmony) <-"cells_type"
levels(Idents(scRNA_harmony))
p1 <- DotPlot(scRNA_harmony, features = markers,
              assay='RNA')

p1
markers_count<- p1$data
write.csv(markers_count,"D:/wangys_uestc_data2/苜蓿数据分析/3.人工注释/dotlpot画图数据.csv")
markers_name<- read.csv("D:/wangys_uestc_data2/苜蓿数据分析/3.人工注释/marker_gene.csv",header = T)
markers_count$features.plot
markers_name$Our_Feature_ID
names(markers_name)
markers_count$Feature_id = markers_name[match(markers_count$features.plot,markers_name$Our_Feature_ID),3]


names(markers_count)

colnames(markers_count)<-c("avg.exp","pct.exp","features.plot","cell_type","avg.exp.scaled",
                           "Feature_id")

markers_count$cell_type<- factor(markers_count$cell_type, levels = rev(c("Lateral root cap", "Cortex","Endodermis","Pericycle","Xylem",
                                                                         "Atrichoblast","Root hairs","Unknow")))

markers_count$Feature_id <- factor(markers_count$Feature_id,levels=unique(markers_count$Feature_id))

cell_type_color <- c("#FF6347","#ffc556","#1E90FF","#9370DB","#0000FF","#00FFFF","#FF00ff","#20B2AA")
pdot <-
  ggplot(markers_count,aes(x = Feature_id,y = cell_type)) +
  geom_point(aes(fill =avg.exp.scaled,size = pct.exp+5),
             color = 'black',
             shape = 21) +
  theme_bw(base_size = 14) +
  xlab('') + ylab('') +
  scale_fill_gradient2(low = 'white',mid = '#EB1D36',high = '#990000',
                       midpoint = 5,
                       limits = c(0,10),
                       #breaks = seq(0,2,10),
                       #labels = seq(0,2,10),
                       name = 'Average expression') +
  scale_size(range = c(1,5),
             limits = c(0, 100),
             breaks = seq(20,100,20),
             labels = seq(20,100,20)
  ) +
  theme(panel.grid = element_blank(),
        legend.title = element_text(size=12,face="bold"),
        legend.text=element_text(size=9,face="bold"),
        axis.text = element_text(color = 'black'),
        aspect.ratio = 0.5,
        plot.margin = margin(t = 0.5,r = 0.5,b = 0.5,l = 0.5,unit = 'cm'),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5, face = 'bold.italic',size=12),
        axis.text.y = element_text(color = "black",size=12,face="bold")
  ) +
  coord_cartesian(clip = 'off') +
  guides(
    fill = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      barwidth = unit(5, "cm")
    ),
    size = guide_legend(
      title = "Percent expressed",
      direction = "horizontal",
      title.position = "top",
      label.position = "bottom",
      override.aes = list(
        color = "black",
        fill = "grey"
      )
    )
  )


pdot
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#处理后有73种差异还比较明显的颜色，基本够用
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) 


p2 <- annoRect(object = pdot,
               annoPos = 'left', # 添加方框的位置
               annoManual = T, # 可以手动设置
               # 设置方框的位置，由于刻度是整数，从下往上依次是1，2，3，4
               # 所以设置方框的位置是0.5，1.5，等，上下包含要框的内容
               # list 表示添加控件的个数，比如下面表示添加 [0.5,4.5]和[4.5,8.5]
               yPosition = list(c(0.5,1.5),c(1.5,2.5)),
               # 个性化多次调整，多尝试，找一个合适的数据值
               xPosition = c(-6,0.3),
               pCol = rep('white',10), # 边框颜色
               pFill = c("#20B2AA","#FF00ff"), # 填充颜色
               alpha = 0.3 # 填充颜色透明度
)

p2

p2_1 <- annoRect(object = p2,
                 annoPos = 'left', # 添加方框的位置
                 annoManual = T, # 可以手动设置
                 # 设置方框的位置，由于刻度是整数，从下往上依次是1，2，3，4
                 # 所以设置方框的位置是0.5，1.5，等，上下包含要框的内容
                 # list 表示添加控件的个数，比如下面表示添加 [0.5,4.5]和[4.5,8.5]
                 yPosition = list(c(2.5,3.5),c(3.5,4.5)),
                 # 个性化多次调整，多尝试，找一个合适的数据值
                 xPosition = c(-6,0.3),
                 pCol = rep('white',10), # 边框颜色
                 pFill = c("#00FFFF","#0000FF"), # 填充颜色
                 alpha = 0.3 # 填充颜色透明度
)

p2_1

p2_2 <- annoRect(object = p2_1,
                 annoPos = 'left', # 添加方框的位置
                 annoManual = T, # 可以手动设置
                 # 设置方框的位置，由于刻度是整数，从下往上依次是1，2，3，4
                 # 所以设置方框的位置是0.5，1.5，等，上下包含要框的内容
                 # list 表示添加控件的个数，比如下面表示添加 [0.5,4.5]和[4.5,8.5]
                 yPosition = list(c(4.5,5.5),c(5.5,6.5)),
                 # 个性化多次调整，多尝试，找一个合适的数据值
                 xPosition = c(-6,0.3),
                 pCol = rep('white',10), # 边框颜色
                 pFill = c("#9370DB","#1E90FF"), # 填充颜色
                 alpha = 0.3 # 填充颜色透明度
)


p2_2

p2_3 <- annoRect(object = p2_2,
                 annoPos = 'left', # 添加方框的位置
                 annoManual = T, # 可以手动设置
                 # 设置方框的位置，由于刻度是整数，从下往上依次是1，2，3，4
                 # 所以设置方框的位置是0.5，1.5，等，上下包含要框的内容
                 # list 表示添加控件的个数，比如下面表示添加 [0.5,4.5]和[4.5,8.5]
                 yPosition = list(c(6.5,7.5),c(7.5,8.5)),
                 # 个性化多次调整，多尝试，找一个合适的数据值
                 xPosition = c(-6,0.3),
                 pCol = rep('white',10), # 边框颜色
                 pFill = c("#ffc556","#FF6347"), # 填充颜色
                 alpha = 0.3 # 填充颜色透明度
)

p2_3


p3_1<- annoSegment(object = p2_3,
                   annoPos = 'top',
                   xPosition = c(2.5,6.5),
                   yPosition =c(8.8,8.8),
                   segWidth = 3,
                   addText = T,
                   textLabel =c("Lateral root cap", "Cortex"),
                   textRot = 45,
                   hjust = 0,
                   vjust = 0,
                   textCol = c(rep("black", 2)),
                   fontface = "bold",
                   pCol = c("#FF6347","#ffc556"),
                   textSize = 12)
p3_1


p3_2<- annoSegment(object = p3_1,
                   annoPos = 'top',
                   xPosition = c(10,13),
                   yPosition =c(8.8,8.8),
                   segWidth = 2,
                   addText = T,
                   textLabel =c("Endodermis","Pericycle"),
                   textRot = 45,
                   hjust = 0,
                   vjust = 0,
                   textCol = c(rep("black", 2)),
                   fontface = "bold",
                   pCol = c("#1E90FF","#9370DB"),
                   textSize = 12)
p3_2

p3_3<- annoSegment(object = p3_2,
                   annoPos = 'top',
                   xPosition = c(17),
                   yPosition =c(8.8,8.8),
                   segWidth = 4,
                   addText = T,
                   textLabel =c("Xylem"),
                   textRot = 45,
                   hjust = 0,
                   vjust = 0,
                   textCol = c(rep("black", 1)),
                   fontface = "bold",
                   pCol = c("#0000FF"),
                   textSize = 12)
p3_3

p3_4<- annoSegment(object = p3_3,
                   annoPos = 'top',
                   xPosition = c(20.5,22.5),
                   yPosition =c(8.8,8.8),
                   segWidth = 1,
                   addText = T,
                   textLabel =c("Atrichoblast","Root hairs","Unknow"),
                   textRot = 45,
                   hjust = 0,
                   vjust = 0,
                   textCol = c(rep("black", 3)),
                   fontface = "bold",
                   pCol = c("#00FFFF","#FF00ff","#20B2AA"),
                   textSize = 12)
p3_4

p4_1 <- annoRect(object = p3_4,
                 annoPos = 'left',
                 annoManual = T,
                 yPosition = c(0.5,8.5),
                 xPosition = c(0.6,4.4),
                 pCol = 'black',
                 pFill = 'transparent',
                 lty = 'dashed',
                 lwd = 1.5)

p4_1

p4_2 <- annoRect(object = p4_1,
                 annoPos = 'left',
                 annoManual = T,
                 yPosition = c(0.5,8.5),
                 xPosition = c(4.6,8.4),
                 pCol = 'black',
                 pFill = 'transparent',
                 lty = 'dashed',
                 lwd = 1.5)

p4_2

p4_3 <- annoRect(object = p4_2,
                 annoPos = 'left',
                 annoManual = T,
                 yPosition = c(0.5,8.5),
                 xPosition = c(8.6,11.4),
                 pCol = 'black',
                 pFill = 'transparent',
                 lty = 'dashed',
                 lwd = 1.5)
p4_3

p4_4 <- annoRect(object = p4_3,
                 annoPos = 'left',
                 annoManual = T,
                 yPosition = c(0.5,8.5),
                 xPosition = c(11.6,14.4),
                 pCol = 'black',
                 pFill = 'transparent',
                 lty = 'dashed',
                 lwd = 1.5)
p4_4 
p4_5 <- annoRect(object = p4_4,
                 annoPos = 'left',
                 annoManual = T,
                 yPosition = c(0.5,8.5),
                 xPosition = c(14.6,19.4),
                 pCol = 'black',
                 pFill = 'transparent',
                 lty = 'dashed',
                 lwd = 1.5)
p4_5
p4_6 <- annoRect(object = p4_5,
                 annoPos = 'left',
                 annoManual = T,
                 yPosition = c(0.5,8.5),
                 xPosition = c(19.6,21.4),
                 pCol = 'black',
                 pFill = 'transparent',
                 lty = 'dashed',
                 lwd = 1.5)
p4_6 

p4_7 <- annoRect(object = p4_6,
                 annoPos = 'left',
                 annoManual = T,
                 yPosition = c(0.5,8.5),
                 xPosition = c(21.6,23.4),
                 pCol = 'black',
                 pFill = 'transparent',
                 lty = 'dashed',
                 lwd = 1.5)
p4_7


ggsave("D:/wangys_uestc_data2/苜蓿数据分析/3.人工注释/marker_dotplot.jpg",p4_7,width =10,height =9)

ggsave("D:/wangys_uestc_data2/苜蓿数据分析/3.人工注释/marker_dotplot.pdf",p4_7,width =10,height =9)

####重新featureplot----


library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidydr)
library(stringr)
library(viridis)

markers <- c("MTR-2g080260","MTR-4g059670","MTR-4g035590","MTR-4g059630",#Lateral root cap
             "MTR-4g088195","MTR-4g088160","MTR-1g074930","MTR-4g101330",#Cortex
             "MTR-1g058560","MTR-7g074650","MTR-4g073140",#Endodermis
             "MTR-8g083150","MTR-4g063090","MTR-8g070770",#Pericycle
             "MTR-1g101830","MTR-2g093990","MTR-7g015960","MTR-4g058860","MTR-7g101080",#Xylem
             "MTR-4g118810","MTR-4g047800", #Atrichoblast
             "MTR-1g075380","MTR-5g090250" #Root hairs
)   

markers <- c("Medtr2g080260","Medtr4g059670","MtBRN","MtRC7",#Lateral root cap
             "MtIFS1","MtIFS3","MtPT5","Medtr4g101330",#Cortex
             "MtCASPL-2","MtSCR","MtPTI11",#Endodermis
             "MtPHO1.1","MtPHO1.2","MtPHO1.3",#Pericycle
             "MtPrx13","Medtr2g093990","Medtr7g015960","Medtr4g058860", "Medtr7g101080",#Xylem
             "Medtr4g118810","Medtr4g047800",#Atrichoblast
             "Medtr1g075380", "Medtr5g090250"#Root hairs
) 

markers <- c("MTR-2g080260",#Lateral root cap
             "MTR-4g088195",#Cortex
             "MTR-1g058560",#Endodermis
             "MTR-8g083150",#Pericycle
             "MTR-1g101830",#Xylem
             "MTR-4g118810", #Atrichoblast
             "MTR-1g075380" #Root hairs
) 


pal <- viridis(n = 10, option = "C")
# "#0D0887FF" "#47039FFF" "#7301A8FF" "#9C179EFF" "#BD3786FF" "#D8576BFF" "#ED7953FF" "#FA9E3BFF" "#FDC926FF" "#F0F921FF"
pal <- viridis(n = 15, option = "D", direction = -1)


p<-FeaturePlot(object =scRNA_harmony,features = markers,reduction = "umap.harmony",raster=FALSE, cols = pal, order = T)
p

i=1
plots=list()
for (i in 1:length(markers)){
  plots[[i]]=FeaturePlot(object=scRNA_harmony,features = markers[i],reduction = "umap.harmony",raster=FALSE, cols = pal, order = T)+
    xlab("UMAP_1")+
    ylab("UMAP_1")+
    theme(legend.position = "none")
}
library(patchwork)
p<-wrap_plots(plots, ncol = 4)+plot_annotation(tag_levels = "A");p

p_legend<-FeaturePlot(object=scRNA_harmony,features = "MTR-2g080260",reduction = "umap.harmony",raster=FALSE, cols = pal, order = T)+
  xlab("UMAP_1")+
  ylab("UMAP_1")+
  theme(legend.position = "right")

legend <- get_legend(p_legend+
                       #guides(color = guide_legend(nrow = 1)) +
                       theme(legend.position = "right"))
plt_featureplot2<- plot_grid(p,legend, ncol=2,rel_widths = c(1,0.1))

ggsave(file="D:/wangys_uestc_data2/苜蓿数据分析/3.人工注释/marker_featureplot/markerfeatureplot_legend.jpg",plt_featureplot2,width=15,height=9)
ggsave(file="D:/wangys_uestc_data2/苜蓿数据分析/3.人工注释/marker_featureplot/markerfeatureplot_legend.pdf",plt_featureplot2,width=15,height=9)

#ggsave("T_cells_CD3E.jpg",device = "pdf",width = 10,height = 10.5,units = "cm")


####细胞比例饼图
#画图使用数据为长数据，第一列group，第二列物种或者细胞类型，第三列为value
Cellratio <- prop.table(table(scRNA_harmony@meta.data$cells_type,scRNA_harmony@meta.data$group1), margin = 2)
Cellratio <- as.data.frame(Cellratio)
View(Cellratio)
colnames(Cellratio) <- c("Cluster","group","Freq")

write.csv(Cellratio,"D:/wangys_uestc_data2/苜蓿数据分析/3.人工注释/细胞比例计算/Cellratio_细胞比例饼图数据.csv")
###加载
library(tidyverse)
library(ggtext)
library(glue)
devtools::install_github("AllanCameron/VoronoiPlus")
library(VoronoiPlus) 
library(ggplot)


au_vor1 <- voronoi_treemap(Freq ~ Cluster,
                           data = Cellratio %>% filter(group=="Mock (9 hpi)")) %>% 
  mutate(type="Mock (9 hpi)")

View(au_vor1)

au_vor2 <- voronoi_treemap(Freq ~ Cluster,
                           data = Cellratio %>% filter(group=="Mock (18 hpi)")) %>% 
  mutate(type="Mock (18 hpi)")

au_vor3 <- voronoi_treemap(Freq ~ Cluster,
                           data = Cellratio %>% filter(group=="Infected (9 hpi)")) %>% 
  mutate(type="Infected (9 hpi)")

au_vor4 <- voronoi_treemap(Freq ~ Cluster,
                           data = Cellratio %>% filter(group=="Infected (18 hpi)")) %>% 
  mutate(type="Infected (18 hpi)")


groups <- au_vor1 %>% bind_rows(au_vor2,au_vor3,au_vor4)

write.csv(groups,"D:/wangys_uestc_data2/苜蓿数据分析/3.人工注释/细胞比例计算/细胞比例饼图数据2.csv")

####################
###修改了比例重新读入
groups<-read.csv("D:/wangys_uestc_data2/苜蓿数据分析/3.人工注释/细胞比例计算/细胞比例饼图数据2.csv",header = T,row.names = 1)

groups$type<- factor(groups$type, levels = c('Mock (9 hpi)', 'Mock (18 hpi)','Infected (9 hpi)',"Infected (18 hpi)"))
groups$group<- factor(groups$group, levels = c("Lateral root cap", "Cortex","Endodermis","Pericycle","Xylem",
                                               "Atrichoblast","Root hairs","Unknow"))

cell_type_color <- c("#FF6347","#ffc556","#1E90FF","#9370DB","#0000FF","#00FFFF","#FF00ff","#20B2AA")


label_group<-groups%>% group_by(group,type,value) %>% summarise(x_median=median(x),y_median=median(y))

nrow(label_group)
View(label_group)

label_group$value <- paste(round(label_group$value *100,2),sep="","%")
pie_plot<- ggplot() +
  geom_polygon(data = groups,mapping = aes(x = x, y = y, group = group, fill = group),
               colour = "white",linewidth =0.8) +
  facet_wrap(.~type,scale="free")+
  scale_fill_manual(values =cell_type_color)+
  geom_text(aes(x_median,y_median,label=value),size=5,color="black",fontface="bold",data = label_group)+
  labs(x=NULL,y=NULL)+
  theme_test()+
  theme(axis.text=element_blank(),
        axis.ticks = element_blank(),
        plot.background = element_rect(fill = "white", colour = "white"),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.spacing.x=unit(-0.1,"cm"),
        panel.spacing.y=unit(-0.1,"cm"),
        plot.margin = margin(2,2,0,0),
        legend.key.height = unit(0.2,"in"), 
        legend.key.width = unit(0.2,"in"),
        legend.title = element_blank(),
        legend.text = element_text(margin = margin(l=0,unit="cm"),size=12),
        strip.background = element_rect(fill="grey90"),
        strip.text = element_text(face="bold",size=15))+theme(legend.position = 'bottom')+
  guides(color = guide_legend(nrow = 4))

pie_plot


ggsave("D:/wangys_uestc_data2/苜蓿数据分析/3.人工注释/细胞比例计算/pie_细胞比例.jpg",pie_plot,width=8,height=9)
ggsave("D:/wangys_uestc_data2/苜蓿数据分析/3.人工注释/细胞比例计算/pie_细胞比例.pdf",pie_plot,width=8,height=9)



####unknow UMAP图----

umap.harmonyredu<- scRNA_harmony@reductions$umap.harmony@cell.embeddings%>%
  as.data.frame() %>% 
  cbind(cell_id = scRNA_harmony@meta.data$cell_id)%>% 
  cbind(cell_type = scRNA_harmony@meta.data$cells_type)%>% 
  cbind(group1= scRNA_harmony@meta.data$group1) %>% 
  cbind(group2 = scRNA_harmony@meta.data$group2)%>%
  cbind(group3 = scRNA_harmony@meta.data$group3)%>%
  cbind(cluster = scRNA_harmony@meta.data$RNA_snn_res.0.4)

names(umap.harmonyredu)

cell_type_color <- c("#FF6347","#ffc556","#1E90FF","#9370DB","#0000FF","#00FFFF","#FF00ff","#20B2AA")

umap.harmonyredu$cell_type<- factor(umap.harmonyredu$cell_type, levels = c("Lateral root cap", "Cortex","Endodermis","Pericycle","Xylem",
                                                                           "Atrichoblast","Root hairs","Unknow" ))


umap.harmonyredu_mock9<-umap.harmonyredu[umap.harmonyredu$group1 %in% c("Mock (9 hpi)"),]

umap.harmonyredu_mock9$cell_type<- factor(umap.harmonyredu_mock9$cell_type,levels = c("Lateral root cap", "Cortex","Endodermis","Pericycle","Xylem",
                                                                                      "Atrichoblast","Root hairs","Unknow"))
mock9_plot <- ggplot(umap.harmonyredu_mock9,aes(x= umapharmony_1, y = umapharmony_2,color = cell_type)) +  
  geom_point(size = 0.5, alpha =1)  +  
  scale_color_manual(values = cell_type_color)+
  xlab("UMAP1")+
  ylab("UMAP2")+
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white"))+
  theme(legend.title = element_blank(), #去掉legend.title 
        legend.key=element_rect(fill='white'), #
        legend.text = element_text(size=15), #设置legend标签的大小
        legend.key.size=unit(1,'cm') ) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes=list(shape = 15,size=5)))+#设置legend中 点的大小
  geom_segment(aes(x = min(umap.harmonyredu$umapharmony_1) , y = min(umap.harmonyredu$umapharmony_2) ,
                   xend = min(umap.harmonyredu$umapharmony_1) +3, yend = min(umap.harmonyredu$umapharmony_2) ),
               colour = "black", size=1,arrow = arrow(length = unit(0.5,"cm")))+ 
  geom_segment(aes(x = min(umap.harmonyredu$umapharmony_1)  , y = min(umap.harmonyredu$umapharmony_2)  ,
                   xend = min(umap.harmonyredu$umapharmony_1) , yend = min(umap.harmonyredu$umapharmony_2) + 3),
               colour = "black", size=1,arrow = arrow(length = unit(0.5,"cm"))) +
  annotate("text", x = min(umap.harmonyredu$umapharmony_1) +1.5, y = min(umap.harmonyredu$umapharmony_2) -1, label = "UMAP_1",
           color="black",size = 5, fontface="bold" ) + 
  annotate("text", x = min(umap.harmonyredu$umapharmony_1) -1, y = min(umap.harmonyredu$umapharmony_2) + 1.5, label = "UMAP_2",
           color="black",size = 5, fontface="bold" ,angle=90)+
  theme(legend.position = "none")

mock9_plot

ggsave('D:/wangys_uestc_data2/苜蓿数据分析/3.人工注释/unknow细胞注释/cell_type细胞注释mock9.jpg',mock9_plot,width =8,height = 8)


umap.harmonyredu_I9<-umap.harmonyredu[umap.harmonyredu$group1 %in% c("Infected (9 hpi)"),]

umap.harmonyredu_I9$cell_type<- factor(umap.harmonyredu_I9$cell_type,levels = c("Lateral root cap", "Cortex","Endodermis","Pericycle","Xylem",
                                                                                "Atrichoblast","Root hairs","Unknow"))
I9_plot <- ggplot(umap.harmonyredu_I9,aes(x= umapharmony_1, y = umapharmony_2,color = cell_type)) +  
  geom_point(size = 0.5, alpha =1)  +  
  scale_color_manual(values = cell_type_color)+
  xlab("UMAP1")+
  ylab("UMAP2")+
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white"))+
  theme(legend.title = element_blank(), #去掉legend.title 
        legend.key=element_rect(fill='white'), #
        legend.text = element_text(size=15), #设置legend标签的大小
        legend.key.size=unit(1,'cm') ) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes=list(shape = 15,size=5)))+#设置legend中 点的大小
  geom_segment(aes(x = min(umap.harmonyredu$umapharmony_1) , y = min(umap.harmonyredu$umapharmony_2) ,
                   xend = min(umap.harmonyredu$umapharmony_1) +3, yend = min(umap.harmonyredu$umapharmony_2) ),
               colour = "black", size=1,arrow = arrow(length = unit(0.5,"cm")))+ 
  geom_segment(aes(x = min(umap.harmonyredu$umapharmony_1)  , y = min(umap.harmonyredu$umapharmony_2)  ,
                   xend = min(umap.harmonyredu$umapharmony_1) , yend = min(umap.harmonyredu$umapharmony_2) + 3),
               colour = "black", size=1,arrow = arrow(length = unit(0.5,"cm"))) +
  annotate("text", x = min(umap.harmonyredu$umapharmony_1) +1.5, y = min(umap.harmonyredu$umapharmony_2) -1, label = "UMAP_1",
           color="black",size = 5, fontface="bold" ) + 
  annotate("text", x = min(umap.harmonyredu$umapharmony_1) -1, y = min(umap.harmonyredu$umapharmony_2) + 1.5, label = "UMAP_2",
           color="black",size = 5, fontface="bold" ,angle=90)+
  theme(legend.position = "none")

I9_plot

ggsave('D:/wangys_uestc_data2/苜蓿数据分析/3.人工注释/unknow细胞注释/cell_type细胞注释I9.jpg',I9_plot,width =8,height = 8)


umap.harmonyredu_mock18<-umap.harmonyredu[umap.harmonyredu$group1 %in% c("Mock (18 hpi)"),]

umap.harmonyredu_mock18$cell_type<- factor(umap.harmonyredu_mock18$cell_type,levels = c("Lateral root cap", "Cortex","Endodermis","Pericycle","Xylem",
                                                                                        "Atrichoblast","Root hairs","Unknow"))
mock18_plot <- ggplot(umap.harmonyredu_mock18,aes(x= umapharmony_1, y = umapharmony_2,color = cell_type)) +  
  geom_point(size = 0.5, alpha =1)  +  
  scale_color_manual(values = cell_type_color)+
  xlab("UMAP1")+
  ylab("UMAP2")+
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white"))+
  theme(legend.title = element_blank(), #去掉legend.title 
        legend.key=element_rect(fill='white'), #
        legend.text = element_text(size=15), #设置legend标签的大小
        legend.key.size=unit(1,'cm') ) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes=list(shape = 15,size=5)))+#设置legend中 点的大小
  geom_segment(aes(x = min(umap.harmonyredu$umapharmony_1) , y = min(umap.harmonyredu$umapharmony_2) ,
                   xend = min(umap.harmonyredu$umapharmony_1) +3, yend = min(umap.harmonyredu$umapharmony_2) ),
               colour = "black", size=1,arrow = arrow(length = unit(0.5,"cm")))+ 
  geom_segment(aes(x = min(umap.harmonyredu$umapharmony_1)  , y = min(umap.harmonyredu$umapharmony_2)  ,
                   xend = min(umap.harmonyredu$umapharmony_1) , yend = min(umap.harmonyredu$umapharmony_2) + 3),
               colour = "black", size=1,arrow = arrow(length = unit(0.5,"cm"))) +
  annotate("text", x = min(umap.harmonyredu$umapharmony_1) +1.5, y = min(umap.harmonyredu$umapharmony_2) -1, label = "UMAP_1",
           color="black",size = 5, fontface="bold" ) + 
  annotate("text", x = min(umap.harmonyredu$umapharmony_1) -1, y = min(umap.harmonyredu$umapharmony_2) + 1.5, label = "UMAP_2",
           color="black",size = 5, fontface="bold" ,angle=90)+
  theme(legend.position = "none")

mock18_plot

ggsave('D:/wangys_uestc_data2/苜蓿数据分析/3.人工注释/unknow细胞注释/cell_type细胞注释mock18.jpg',mock18_plot,width =8,height = 8)


umap.harmonyredu_I18<-umap.harmonyredu[umap.harmonyredu$group1 %in% c("Infected (18 hpi)"),]

umap.harmonyredu_I18$cell_type<- factor(umap.harmonyredu_I18$cell_type,levels = c("Lateral root cap", "Cortex","Endodermis","Pericycle","Xylem",
                                                                                  "Atrichoblast","Root hairs","Unknow"))
I18_plot <- ggplot(umap.harmonyredu_I18,aes(x= umapharmony_1, y = umapharmony_2,color = cell_type)) +  
  geom_point(size = 0.5, alpha =1)  +  
  scale_color_manual(values = cell_type_color)+
  xlab("UMAP1")+
  ylab("UMAP2")+
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white"))+
  theme(legend.title = element_blank(), #去掉legend.title 
        legend.key=element_rect(fill='white'), #
        legend.text = element_text(size=15), #设置legend标签的大小
        legend.key.size=unit(1,'cm') ) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes=list(shape = 15,size=5)))+#设置legend中 点的大小
  geom_segment(aes(x = min(umap.harmonyredu$umapharmony_1) , y = min(umap.harmonyredu$umapharmony_2) ,
                   xend = min(umap.harmonyredu$umapharmony_1) +3, yend = min(umap.harmonyredu$umapharmony_2) ),
               colour = "black", size=1,arrow = arrow(length = unit(0.5,"cm")))+ 
  geom_segment(aes(x = min(umap.harmonyredu$umapharmony_1)  , y = min(umap.harmonyredu$umapharmony_2)  ,
                   xend = min(umap.harmonyredu$umapharmony_1) , yend = min(umap.harmonyredu$umapharmony_2) + 3),
               colour = "black", size=1,arrow = arrow(length = unit(0.5,"cm"))) +
  annotate("text", x = min(umap.harmonyredu$umapharmony_1) +1.5, y = min(umap.harmonyredu$umapharmony_2) -1, label = "UMAP_1",
           color="black",size = 5, fontface="bold" ) + 
  annotate("text", x = min(umap.harmonyredu$umapharmony_1) -1, y = min(umap.harmonyredu$umapharmony_2) + 1.5, label = "UMAP_2",
           color="black",size = 5, fontface="bold" ,angle=90)+
  theme(legend.position = "none")

I18_plot

ggsave('D:/wangys_uestc_data2/苜蓿数据分析/3.人工注释/unknow细胞注释/cell_type细胞注释I18.jpg',I18_plot,width =8,height = 8)

library(cowplot)

cell_type_merged1 <-plot_grid(mock9_plot,mock18_plot,
                              I9_plot,I18_plot,
                              ncol = 2,
                              labels = c('Mock 9hpi','Mock 18hpi',
                                         "Infected 9hpi","Infected 18hpi"),
                              label_size = 20,
                              rel_widths = c(1,1,1,1),
                              rel_heights = c(1,1,1,1))


ggsave('D:/wangys_uestc_data2/苜蓿数据分析/3.人工注释/unknow细胞注释/侵染未侵染四张图cell_type_merged.jpg',cell_type_merged1,width =15,height = 15)
ggsave('D:/wangys_uestc_data2/苜蓿数据分析/3.人工注释/unknow细胞注释/侵染未侵染四张图cell_type_merged.pdf',cell_type_merged1,width =15,height = 15)


