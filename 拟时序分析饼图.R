ks_run_Monocle2 <- function(object, #seurat obj or expression matrix (建议数据格式转为matrix,如果数据量大转化为稀疏矩阵as(as.matrix(data), "sparseMatrix"))
                            layer, #used when object is a seurat obj
                            assay, #used when object is a seurat obj
                            lowerDetectionLimit = 0.1, 
                            VARgenesM=c("dispersionTable","seurat","differentialGeneTest"),
                            cellAnno=NULL, 
                            define_root=F,
                            root_state,
                            reverse=NULL
){
  
  
  if(class(object)[1] == 'Seurat') {
    
    data <- GetAssayData(object=object, layer=layer, assay=assay)#get expression matrix data
    
    pd <- new("AnnotatedDataFrame", data = object@meta.data)
    fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
    fd <- new("AnnotatedDataFrame", data = fData)
    
    if(all(data == floor(data))) {
      expressionFamily <- negbinomial.size()
    } else if(any(data < 0)){
      expressionFamily <- uninormal()
    } else {
      expressionFamily <- tobit()
    }
    
    #Creates a new CellDateSet object.
    monocle_cds <- newCellDataSet(data,
                                  phenoData = pd, 
                                  featureData = fd,
                                  lowerDetectionLimit=0.1,
                                  expressionFamily=expressionFamily)
    
  }else{
    
    print("This fucntions only apply for a seurat obj")
  }
  
  
  
  #Estimate size factors and dispersions
  #数据处理
  monocle_cds <- estimateSizeFactors(monocle_cds)#size facotr标准化细胞之间的mRNA差异
  monocle_cds <- estimateDispersions(monocle_cds)
  
  #质量控制-filter cells
  monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)
  # print(head(fData(monocle_cds)))
  # print(head(pData(monocle_cds)))
  # expressed_genes <- row.names(subset(fData(mouse_monocle), num_cells_expressed >= 10))
  monocle_cds <- monocle_cds[fData(monocle_cds)$num_cells_expressed >= 10, ]
  
  
  #select methods for VariableFeatures
  if(VARgenesM=="dispersionTable"){
    
    disp_table <- dispersionTable(monocle_cds)
    ordering_genes <- subset(disp_table,
                             mean_expression >= 0.1 &
                               dispersion_empirical >= 1.5* dispersion_fit)$gene_id
    
  }
  
  
  if(VARgenesM=="seurat"){
    
    ordering_genes <- VariableFeatures(FindVariableFeatures(object, assay = "RNA"), assay = "RNA")
    
  }
  
  
  if(VARgenesM=="differentialGeneTest"){
    
    diff_test_res <- differentialGeneTest(monocle_cds,fullModelFormulaStr = paste0("~",cellAnno))##~后面是表示对谁做差异分析的变量
    diff_test_res_sig <- diff_test_res[order(diff_test_res$qval,decreasing=F),]
    
    ordering_sce <- diff_test_res_sig[diff_test_res_sig$qval< 0.01,]
    
    if(nrow(ordering_sce)>3000){
      
      ordering_genes <- ordering_sce$gene_short_name[1:3000]
      
    }else{
      
      ordering_genes <- rdering_sce$gene_short_name
    }
    
  }
  
  #Marks genes for clustering
  monocle_cds <- setOrderingFilter(monocle_cds, ordering_genes)
  plot_ordering_genes(monocle_cds)
  
  
  #cluster
  monocle_cds <- reduceDimension(monocle_cds, max_components = 2,reduction_method = 'DDRTree')
  
  #order cells
  monocle_cds <- orderCells(monocle_cds, reverse=reverse)
  
  if(define_root){
    monocle_cds <- monocle_cds <- orderCells(monocle_cds,root_state = root_state)
  }
  
  
  return(monocle_cds)
  
}



#monocle2 拟时序分析----

scRNA<- readRDS("G:/苜蓿数据分析/Medicago_data/rds数据/Merge_Medicago_注释.rds")
Cortex_Unknow_scrna <- subset(scRNA, subset = cells_type %in% c("Cortex", "Unknow"))

expressed_genes <- rowSums(GetAssayData(Cortex_Unknow_scrna, slot = "counts") > 0) > 30
scRNA_filtered <- subset(Cortex_Unknow_scrna, features = names(expressed_genes[expressed_genes]))
scRNA_filtered <- subset(scRNA_filtered, subset = nFeature_RNA > 500)
filtered_counts2 <- GetAssayData(scRNA_filtered, slot = "counts")
ncol(filtered_counts2)

set.seed(123)  # 保证可复现
sample_cells <- sample(colnames(filtered_counts2), 20000)
final_counts <- filtered_counts2[, sample_cells]

meta_data <- scRNA_filtered@meta.data[sample_cells, , drop = FALSE]
new_seurat <- CreateSeuratObject(counts = final_counts, meta.data = meta_data)
ncol(new_seurat)
unique(new_seurat$cells_type)
sce_CDS <- ks_run_Monocle2(object =new_seurat,
                           layer = 'counts',
                           assay = "RNA",
                           VARgenesM="seurat",
                           cellAnno = "cells_type",
                           root_state ="Cortex")
saveRDS(sce_CDS,"G:/苜蓿数据分析/6.拟时序分析/monocle2/Cortex_Unknow/reroot_Cortex_Unknow拟时序20000.rds")


#monocle饼图绘制-----
library(monocle)
library(Seurat)
library(tidyverse)
library(harmony)
library(dplyr)
library(patchwork)
library(SeuratWrappers)

sce_CDS<-readRDS("G:/苜蓿数据分析/6.拟时序分析/monocle2/Cortex_Unknow/reroot_Cortex_Unknow拟时序20000.rds")


#提取数据=======================================================================
data_df <- t(reducedDimS(sce_CDS)) %>% as.data.frame() %>% #提取坐标  
  select_(Component_1 = 1, Component_2 = 2) %>% #重命名  
  rownames_to_column("cells") %>% #rownames命名  
  mutate(pData(sce_CDS)$State) %>% #添加State  
  mutate(pData(sce_CDS)$Pseudotime,          
         pData(sce_CDS)$orig.ident,          
         pData(sce_CDS)$cells_type)#将这些需要作图的有用信息都添加上

library(dplyr)  # Ensure tidyverse loaded

data_df <- t(reducedDimS(sce_CDS)) %>%
  as.data.frame() %>%
  select(Component_1 = 1, Component_2 = 2) %>%  # Fixed: select() not select_()
  rownames_to_column("cells") %>%
  left_join(  # Better: join pData once (faster/safer than mutate with external vectors)
    as.data.frame(pData(sce_CDS)) %>% 
      select(State, Pseudotime, orig.ident, cells_type) %>%
      rownames_to_column("cells"),
    by = "cells"
  )


colnames(data_df) <- c("cells","Component_1","Component_2","State",                       
                       "Pseudotime","orig.ident","celltype")

View(data_df)
#==============================================================================
#轨迹数据提取---完全摘录于monocle包原函数
dp_mst <- minSpanningTree(sce_CDS)
reduced_dim_coords <- reducedDimK(sce_CDS)
ica_space_df <- Matrix::t(reduced_dim_coords) %>% as.data.frame() %>%   
  select_(prin_graph_dim_1 = 1, prin_graph_dim_2 = 2) %>%   
  mutate(sample_name = rownames(.), sample_state = rownames(.))
#构建一个做轨迹线图的数据
edge_df <- dp_mst %>% igraph::as_data_frame() %>%   
  select_(source = "from", target = "to") %>%   
  left_join(ica_space_df %>% select_(source = "sample_name",                                      
                                     source_prin_graph_dim_1 = "prin_graph_dim_1",                                     
                                     source_prin_graph_dim_2 = "prin_graph_dim_2"), by = "source") %>%   
  left_join(ica_space_df %>% select_(target = "sample_name",                                      
                                     target_prin_graph_dim_1 = "prin_graph_dim_1",                                      
                                     target_prin_graph_dim_2 = "prin_graph_dim_2"), by = "target")

##这里以state为例，计算样本占比：
#计算细胞比例
Cellratio <- prop.table(table(data_df$State, data_df$celltype), margin = 2)
#计算各组样本不同细胞群比例
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c('State',"celltype","Freq")
write.csv(Cellratio,"G:/苜蓿数据分析/6.拟时序分析/monocle2/Cortex_Unknow/Cellratio_Cortex_Unknow拟时序retoot20000.csv")

##ggplot2作图：
#ggplot作图
library(ggplot2)
library(tidydr)
library(ggforce)
library(ggrastr)
library(viridis)
g <- ggplot() +   
  geom_point_rast(data = data_df, aes(x = Component_1,                                  
                                      y = Component_2,                                 
                                      color =Pseudotime)) + #散点图  
  scale_color_viridis()+#密度色  
  geom_segment(aes_string(x = "source_prin_graph_dim_1",                          
                          y = "source_prin_graph_dim_2",                           
                          xend = "target_prin_graph_dim_1",                           
                          yend = "target_prin_graph_dim_2"),                
               linewidth = 1,                
               linetype = "solid", na.rm = TRUE, data = edge_df)+#添加轨迹线  
  theme_dr(arrow = grid::arrow(length = unit(0, "inches")))+#坐标轴主题修改  
  theme(panel.grid.major = element_blank(),        
        panel.grid.minor = element_blank())+  
  # geom_arc(arrow = arrow(length = unit(0.15, "inches"), #曲线箭头                         
  #                        type = "closed",angle=30),           
  #          aes(x0=0,y0=-3,r=5, start=-0.4*pi, end=0.4*pi),lwd=1)+  
  geom_arc_bar(data=subset(Cellratio,State=='1'),stat = "pie",#添加饼图
               aes(x0=-10,y0=0,r0=0,r=2,amount=Freq,fill=celltype))+
  # geom_arc_bar(data=subset(Cellratio,State=='2'),stat = "pie",
  #              aes(x0=0,y0=9,r0=0,r=2.5,amount=Freq,fill=celltype))+
  geom_arc_bar(data=subset(Cellratio,State=='3'),stat = "pie",
               aes(x0=0,y0=5,r0=0,r=2,amount=Freq,fill=celltype))+
  geom_arc_bar(data=subset(Cellratio,State=='4'),stat = "pie",
               aes(x0=5,y0=-5,r0=0,r=2,amount=Freq,fill=celltype))+
  # geom_arc_bar(data=subset(Cellratio,State=='5'),stat = "pie",
  #              aes(x0=10,y0=-8,r0=0,r=2.5,amount=Freq,fill=celltype))+
  geom_arc_bar(data=subset(Cellratio,State=='6'),stat = "pie",
               aes(x0=10,y0=5,r0=0,r=2,amount=Freq,fill=celltype))+
  # geom_arc_bar(data=subset(Cellratio,State=='7'),stat = "pie",
  #              aes(x0=10,y0=-8,r0=0,r=2.5,amount=Freq,fill=celltype))+
  scale_fill_manual(values = c("#FFA500","#20B2AA"))

g

colour=c("#9370DB","#20B2AA")

colour=c("#FFA500","#98FB98","#8B008B","#40E0D0","#EE82EE")

ggsave("G:/苜蓿数据分析/6.拟时序分析/monocle2/Cortex_Unknow/Cortex_Unknow拟时序20000reroot.pdf",g,width = 8,height=7)



####Pericycle和Unknow拟时序分析----
library(monocle)
library(Seurat)
library(tidyverse)
library(harmony)
library(dplyr)
library(patchwork)
library(SeuratWrappers)

scRNA<- readRDS("E:/苜蓿数据分析/Medicago_data/rds数据/Merge_Medicago_注释.rds")

Idents(scRNA_harmony)<-"cells_type"   ###将cells_type 付给Idents

Pericycle_Unknow_scrna <- subset(scRNA, subset = cells_type %in% c("Pericycle", "Unknow"))

ncol(Pericycle_Unknow_scrna)

Pericycle_Unknow_scrna$cells_type

expressed_genes <- rowSums(GetAssayData(Pericycle_Unknow_scrna, slot = "counts") > 0) > 30
scRNA_filtered <- subset(Pericycle_Unknow_scrna, features = names(expressed_genes[expressed_genes]))
scRNA_filtered <- subset(scRNA_filtered, subset = nFeature_RNA > 500)

filtered_counts2 <- GetAssayData(scRNA_filtered, slot = "counts")

ncol(scRNA_filtered)#36436

# 随机抽样细胞（假设抽样 10000 个）
set.seed(123)  # 保证可复现
sample_cells <- sample(colnames(filtered_counts2), 36434)
final_counts <- filtered_counts2[, sample_cells]

meta_data <- scRNA_filtered@meta.data[sample_cells, , drop = FALSE]
new_seurat <- CreateSeuratObject(counts = final_counts, meta.data = meta_data)
ncol(new_seurat)
new_seurat$cells_type
sce_CDS <- ks_run_Monocle2(object =new_seurat,
                           layer = 'counts',
                           assay = "RNA",
                           VARgenesM="seurat",
                           cellAnno = "cells_type",
                           root_state ="Pericycle")

unique(sce_CDS$cells_type)
saveRDS(sce_CDS,"/home/zqwangyansu/muxu/sc_data/Pericycle_Unknow/Pericycle_Unknow拟时序36434.rds")


saveRDS(sce_CDS,"E:/苜蓿数据分析/6.拟时序分析/monocle2/Pericycle_Unknow/Pericycle_Unknow拟时序20000.rds")

sce_CDS<- readRDS("E:/苜蓿数据分析/6.拟时序分析/monocle2/Pericycle_Unknow/Pericycle_Unknow拟时序36434.rds")

sce_CDS<- readRDS("E:/苜蓿数据分析/6.拟时序分析/monocle2/Pericycle_Unknow/Pericycle_Unknow拟时序20000.rds")


sce_CDS<- readRDS("/home/zqwangyansu/muxu/sc_data/Pericycle_Unknow/Pericycle_Unknow拟时序36434.rds")

# sce_CDS<- readRDS("/home/zqwangyansu/muxu/sc_data/Pericycle_Unknow/Pericycle_Unknow拟时序20000.rds")

plot_cell_trajectory(sce_CDS,color_by = "Pseudotime")
plot_cell_trajectory(sce_CDS,color_by = "State")
plot_cell_trajectory(sce_CDS,color_by = "cells_type")

sce_CDS

####ggplot2个性可视化monocle2结果----

#提取数据=======================================================================
data_df <- t(reducedDimS(sce_CDS)) %>% as.data.frame() %>% #提取坐标  
  select_(Component_1 = 1, Component_2 = 2) %>% #重命名  
  rownames_to_column("cells") %>% #rownames命名  
  mutate(pData(sce_CDS)$State) %>% #添加State  
  mutate(pData(sce_CDS)$Pseudotime,          
         pData(sce_CDS)$orig.ident,          
         pData(sce_CDS)$cells_type)#将这些需要作图的有用信息都添加上
colnames(data_df) <- c("cells","Component_1","Component_2","State",                       
                       "Pseudotime","orig.ident","celltype")
#==============================================================================
#轨迹数据提取---完全摘录于monocle包原函数
dp_mst <- minSpanningTree(sce_CDS)
reduced_dim_coords <- reducedDimK(sce_CDS)
ica_space_df <- Matrix::t(reduced_dim_coords) %>% as.data.frame() %>%   
  select_(prin_graph_dim_1 = 1, prin_graph_dim_2 = 2) %>%   
  mutate(sample_name = rownames(.), sample_state = rownames(.))
#构建一个做轨迹线图的数据
edge_df <- dp_mst %>% igraph::as_data_frame() %>%   
  select_(source = "from", target = "to") %>%   
  left_join(ica_space_df %>% select_(source = "sample_name",                                      
                                     source_prin_graph_dim_1 = "prin_graph_dim_1",                                     
                                     source_prin_graph_dim_2 = "prin_graph_dim_2"), by = "source") %>%   
  left_join(ica_space_df %>% select_(target = "sample_name",                                      
                                     target_prin_graph_dim_1 = "prin_graph_dim_1",                                      
                                     target_prin_graph_dim_2 = "prin_graph_dim_2"), by = "target")

##这里以state为例，计算样本占比：
#计算细胞比例
Cellratio <- prop.table(table(data_df$State, data_df$celltype), margin = 2)
#计算各组样本不同细胞群比例
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c('State',"celltype","Freq")

write.csv(Cellratio,"G:/苜蓿数据分析/6.拟时序分析/monocle2/Pericycle_Unknow/Cellratio_Pericycle_Unknow拟时序36424_pp2.csv")
##ggplot2作图：
#ggplot作图
library(ggplot2)
library(tidydr)
library(ggforce)
library(ggrastr)
library(viridis)
g <- ggplot() +   
  geom_point_rast(data = data_df, aes(x = Component_1,                                  
                                      y = Component_2,                                 
                                      color =Pseudotime)) + #散点图  
  scale_color_viridis()+#密度色  
  geom_segment(aes_string(x = "source_prin_graph_dim_1",                          
                          y = "source_prin_graph_dim_2",                           
                          xend = "target_prin_graph_dim_1",                           
                          yend = "target_prin_graph_dim_2"),                
               linewidth = 1,                
               linetype = "solid", na.rm = TRUE, data = edge_df)+#添加轨迹线  
  theme_dr(arrow = grid::arrow(length = unit(0, "inches")))+#坐标轴主题修改  
  theme(panel.grid.major = element_blank(),        
        panel.grid.minor = element_blank())+  
  geom_arc(arrow = arrow(length = unit(0.15, "inches"), #曲线箭头                         
                         type = "closed",angle=30),           
           aes(x0=0,y0=-3,r=5, start=-0.4*pi, end=0.4*pi),lwd=1)+  
  geom_arc_bar(data=subset(Cellratio,State=='1'),stat = "pie",#添加饼图            
               aes(x0=-15,y0=0,r0=0,r=2.5,amount=Freq,fill=celltype))+  
  geom_arc_bar(data=subset(Cellratio,State=='2'),stat = "pie",               
               aes(x0=2,y0=9,r0=0,r=2.5,amount=Freq,fill=celltype))+  
  geom_arc_bar(data=subset(Cellratio,State=='3'),stat = "pie",               
               aes(x0=10,y0=-8,r0=0,r=2.5,amount=Freq,fill=celltype))+  
  scale_fill_manual(values = c("#9370DB","#20B2AA"))

g

colour=c("#9370DB","#20B2AA")

colour=c("#FFA500","#98FB98","#8B008B","#40E0D0","#EE82EE")

ggsave("G:/苜蓿数据分析/6.拟时序分析/monocle2/Pericycle_Unknow/Pericycle_Unknow拟时序36424_pp2.pdf",g,width = 8,height=7)


