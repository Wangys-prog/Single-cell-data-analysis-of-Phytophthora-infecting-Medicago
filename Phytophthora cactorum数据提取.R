####提取Phytophthora cactorum数据----
######苜蓿数据
umap.harmonyredu<- scRNA_harmony@reductions$umap.harmony@cell.embeddings%>%
  as.data.frame() %>% 
  cbind(cell_id = scRNA_harmony@meta.data$cell_id)%>% 
  cbind(cell_type = scRNA_harmony@meta.data$cells_type)%>% 
  cbind(group1= scRNA_harmony@meta.data$group1) %>% 
  cbind(group2 = scRNA_harmony@meta.data$group2)%>%
  cbind(group3 = scRNA_harmony@meta.data$group3)%>%
  cbind(cluster = scRNA_harmony@meta.data$RNA_snn_res.0.4)

phy_scRNA_object<- readRDS("G:/苜蓿数据分析/卵菌矩阵数据/raw_all.rds")
View(phy_scRNA_object)
metadata<- FetchData(phy_scRNA_object,"orig.ident")
metadata$cell_id <- rownames(metadata)
phy_scRNA_object <- AddMetaData(phy_scRNA_object,metadata = metadata)
phy_scRNA_object@meta.data$cell_id

Phytophthora_data<-phy_scRNA_object@meta.data %>%
  as.data.frame()

colnames(Phytophthora_data)
head(Phytophthora_data)

umap.harmonyredu$Phy_nCount_RNA<- Phytophthora_data[match(umap.harmonyredu$cell_id,Phytophthora_data$cell_id),2]
umap.harmonyredu$Phy_nFeature_RNA<- Phytophthora_data[match(umap.harmonyredu$cell_id,Phytophthora_data$cell_id),3]

setdiff(umap.harmonyredu$cell_id,Phytophthora_data$cell_id)

colnames(umap.harmonyredu)
head(umap.harmonyredu)

umap.harmonyredu[is.na(umap.harmonyredu)]=0

# umap.harmonyredu_filter<-umap.harmonyredu[!is.na(umap.harmonyredu$Phy_nCount_RNA),]
umap.harmonyredu_filter<-umap.harmonyredu[!umap.harmonyredu$Phy_nFeature_RNA=="0",]

library(viridis)

pal <- viridis(n = 15, option = "D", direction = -1)

####Infected (9 hpi)
umap.harmonyredu$group2
umap.harmonyredu_I9_phy<-umap.harmonyredu[umap.harmonyredu$group2 %in% c("Infected (9 hpi)_replicate2"),]

umap.harmonyredu_I9_phy_filter<-umap.harmonyredu_filter[umap.harmonyredu_filter$group2 %in% c("Infected (9 hpi)_replicate2"),]
View(umap.harmonyredu_I9_phy_filter)
write.csv(umap.harmonyredu_I9_phy_filter,"D:/wangys_uestc_data2/苜蓿数据分析/3.人工注释/unknow细胞注释/umap.harmonyredu_I9_phy_filter_rep2.csv")
umap.harmonyredu_I9_phy$cell_type<- factor(umap.harmonyredu_I9_phy$cell_type,levels = c("Lateral root cap", "Cortex","Endodermis","Pericycle","Xylem",
                                                                                        "Atrichoblast","Root hairs","Unknow"))
I9_plot_phy <- ggplot()+ 
  geom_point(data=umap.harmonyredu_I9_phy,aes(x=umapharmony_1, y = umapharmony_2),
             size = 1, alpha =1,fill="#c1cdc1",color="#c1cdc1")+
  
  theme_bw()+
  geom_point(data=umap.harmonyredu_I9_phy_filter,aes(x= umapharmony_1, y = umapharmony_2,color =Phy_nCount_RNA),
             size = 5, alpha =1)+
  scale_colour_gradient2(low = "yellow",mid="blue", high = "red")+
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
        legend.key.size=unit(0.5,'cm') ) +  # 设置legend标签之间的大小
  #guides(color = guide_legend(override.aes=list(size=5)))+#设置legend中 点的大小
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
  theme(legend.position = "right")

I9_plot_phy


######Infected (18 hpi)_rep1
colnames(umap.harmonyredu_I18_phy)
umap.harmonyredu_I18_phy_rep1<-umap.harmonyredu[umap.harmonyredu$group2 %in% c("Infected (18 hpi)_replicate1"),]

umap.harmonyredu_I18_phy_filter_rep1<-umap.harmonyredu_filter[umap.harmonyredu_filter$group2 %in% c("Infected (18 hpi)_replicate1"),]
write.csv(umap.harmonyredu_I18_phy_filter_rep1,"D:/wangys_uestc_data2/苜蓿数据分析/3.人工注释/unknow细胞注释/umap.harmonyredu_I18_phy_filter_rep1.csv")

I18_plot_phy_rep1 <- ggplot()+ 
  geom_point(data=umap.harmonyredu_I18_phy_rep1,aes(x=umapharmony_1, y = umapharmony_2),
             size = 1, alpha =1,fill="#c1cdc1",color="#c1cdc1")+
  
  theme_bw()+
  geom_point(data=umap.harmonyredu_I18_phy_filter_rep1,aes(x= umapharmony_1, y = umapharmony_2,color = Phy_nCount_RNA),
             size = 5, alpha =1)+
  scale_colour_gradient2(low = "yellow",mid="blue", high = "red")+
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
  theme(#legend.title = element_blank(), #去掉legend.title 
    legend.key=element_rect(fill='white'), #
    legend.text = element_text(size=10), #设置legend标签的大小
    legend.key.size=unit(0.5,'cm') ) +  # 设置legend标签之间的大小
  #guides(color = guide_legend(override.aes=list(size=5)))+#设置legend中 点的大小
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
  theme(legend.position = "right")

I18_plot_phy_rep1


######Infected (18 hpi)_rep1
colnames(umap.harmonyredu_I18_phy)
umap.harmonyredu_I18_phy_rep2<-umap.harmonyredu[umap.harmonyredu$group2 %in% c("Infected (18 hpi)_replicate2"),]

umap.harmonyredu_I18_phy_filter_rep2<-umap.harmonyredu_filter[umap.harmonyredu_filter$group2 %in% c("Infected (18 hpi)_replicate2"),]
write.csv(umap.harmonyredu_I18_phy_filter_rep2,"D:/wangys_uestc_data2/苜蓿数据分析/3.人工注释/unknow细胞注释/umap.harmonyredu_I18_phy_filter_rep2.csv")

I18_plot_phy_rep2 <- ggplot()+ 
  geom_point(data=umap.harmonyredu_I18_phy_rep2,aes(x=umapharmony_1, y = umapharmony_2),
             size = 1, alpha =1,fill="#c1cdc1",color="#c1cdc1")+
  
  theme_bw()+
  geom_point(data=umap.harmonyredu_I18_phy_filter_rep2,aes(x= umapharmony_1, y = umapharmony_2,color = Phy_nCount_RNA),
             size = 5, alpha =1)+
  scale_colour_gradient2(low = "yellow",mid="blue", high = "red")+
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
  theme(#legend.title = element_blank(), #去掉legend.title 
    legend.key=element_rect(fill='white'), #
    legend.text = element_text(size=10), #设置legend标签的大小
    legend.key.size=unit(0.5,'cm') ) +  # 设置legend标签之间的大小
  #guides(color = guide_legend(override.aes=list(size=5)))+#设置legend中 点的大小
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
  theme(legend.position = "right")

I18_plot_phy_rep2


library(cowplot)

cell_type_merged2 <-plot_grid(I9_plot_phy,I18_plot_phy_rep1,
                              I18_plot_phy_rep2,
                              ncol = 2,
                              labels = c('Infected 9hpi_replicate2',"Infected 18hpi_replicate1",
                                         "Infected 18hpi_replicate2"),
                              label_size = 20,
                              rel_widths = c(1.2,1.2,1.2),
                              rel_heights = c(1.2,1.2,1.2))


ggsave('D:/wangys_uestc_data2/苜蓿数据分析/3.人工注释/unknow细胞注释/Phytophthora_cactorum侵染未侵染四张图_merged.jpg',cell_type_merged2,width =15,height = 12)
ggsave('D:/wangys_uestc_data2/苜蓿数据分析/3.人工注释/unknow细胞注释/Phytophthora_cactorum侵染未侵染四张图_merged.pdf',cell_type_merged2,width =15,height = 12)

####Phytophthora cactorum所在位置堆砌柱状图----
library(ggplot2)
library(ggprism)
library(ggalluvial)
cell_data<- read.table("clipboard",header = T,sep="\t")
colnames(cell_data)
cell_data_long<-cell_data %>% 
  pivot_longer(-X)

colnames(cell_data_long)<-c("samples","cell_type","value")

cell_data_long$cell_type<- factor(cell_data_long$cell_type,levels= c("Lateral.root.cap", "Cortex","Endodermis","Pericycle","Xylem",
                                                                     "Unknow"))

cell_type_color <- c("#FF6347","#ffc556","#1E90FF","#9370DB","#0000FF","#20B2AA")


cell_data_long$samples <- factor(cell_data_long$samples,levels=c("I_9_2",
                                                                 "I_18_1",
                                                                 "I_18_2"))

cell_data_plot <- ggplot(data=cell_data_long,aes(samples,value,fill=cell_type,stratum =cell_type, alluvium = cell_type)) +
  geom_stratum(color="black",width=0.6,size=1)+
  geom_flow(alpha = 0.5) +  #绘制同类别之间的连接线
  scale_fill_manual(values=cell_type_color) +
  scale_y_continuous(name = "Relative abundance (%)")+
  theme(
    axis.title=element_text(size=15,face="Arial",color="black"),
    axis.text = element_text(size=15,face="Arial",color="black"),
    legend.title=element_text(size=15,face="Arial",color="black"),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black", size = 0.5))+theme_bw()+
  theme(text=element_text(family="B",size=20))+
  theme(axis.ticks.length=unit(-0.25, "cm"),
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"),angle=45),
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
  )+
  guides(fill=guide_legend(keywidth = 1, keyheight = 1)) +
  theme_prism(base_fontface = "plain",
              base_size = 18,
              base_line_size = 1,
              axis_text_angle = 45)+
  theme(legend.position = 'right',
        axis.title.x = element_blank())
cell_data_plot

#######################

cluster_data<- read.table("clipboard",header = T,sep="\t")
colnames(cluster_data)
cluster_data_long<-cluster_data %>% 
  pivot_longer(-X)

colnames(cluster_data_long)<-c("samples","cluster","value")

cluster_data_long$cluster<- factor(cluster_data_long$cluster,levels= c("sc_RNA_C0","sc_RNA_C1","sc_RNA_C2",	
                                                                       "sc_RNA_C3","sc_RNA_C5","sc_RNA_C6",	
                                                                       "sc_RNA_C7","sc_RNA_C8","sc_RNA_C10",	
                                                                       "sc_RNA_C11","sc_RNA_C15"))

cluster_type_color<- c("#20B2AA","#FFA500","#9370DB","#98FB98","#1E90FF","#FFFF00", "#808000","#FF00FF","#7B68EE",
                       "#9400D3","seagreen4")


cluster_data_long$samples <- factor(cluster_data_long$samples,levels=c("I_9_2",
                                                                       "I_18_1",
                                                                       "I_18_2"))

cluster_data_plot <- ggplot(data=cluster_data_long,aes(samples,value,fill=cluster,stratum =cluster, alluvium = cluster)) +
  geom_stratum(color="black",width=0.6,size=1)+
  geom_flow(alpha = 0.5) +  #绘制同类别之间的连接线
  scale_fill_manual(values=cluster_type_color) +
  scale_y_continuous(name = "Relative abundance (%)")+
  theme(
    axis.title=element_text(size=15,face="Arial",color="black"),
    axis.text = element_text(size=15,face="Arial",color="black"),
    legend.title=element_text(size=15,face="Arial",color="black"),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black", size = 0.5))+theme_bw()+
  theme(text=element_text(family="B",size=20))+
  theme(axis.ticks.length=unit(-0.25, "cm"),
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"),angle=45),
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
  )+
  guides(fill=guide_legend(keywidth = 1, keyheight = 1)) +
  theme_prism(base_fontface = "plain",
              base_size = 18,
              base_line_size = 1,
              axis_text_angle = 45)+
  theme(legend.position = 'right',
        axis.title.x = element_blank())
cluster_data_plot

library(cowplot)
cell_type_merged2 <-plot_grid(cell_data_plot,cluster_data_plot,
                              ncol = 2,
                              # labels = c('A',"Infected 18hpi_replicate1",
                              #            "Infected 18hpi_replicate2"),
                              label_size = 20,
                              rel_widths = c(1.5,1.5),
                              rel_heights = c(1,5,1.5))


ggsave('D:/wangys_uestc_data2/苜蓿数据分析/3.人工注释/unknow细胞注释/Phytophthora_cactorum位置堆砌柱状图_merged.jpg',cell_type_merged2,width =10,height = 5)
ggsave('D:/wangys_uestc_data2/苜蓿数据分析/3.人工注释/unknow细胞注释/Phytophthora_cactorum位置堆砌柱状图_merged.pdf',cell_type_merged2,width =10,height = 5)
