# 1.Loading packages

library(ggsci)
library(viridis)
library(ggplot2)
library(RCurl)
library(cowplot)
library(dplyr)
library(Seurat)

# 2.loadng ST-seq datasets

ST <- integrated_data
# ST <- readRDS("integrated_data.rds")
ST # 54582 features across 16993 samples within 3 assays 
ST$sample <-  factor(ST$sample, levels = c("Control_1", "Control_2","Control_3","Control_4","Control_5","Control_6"))
ST$type = factor(ST$type, levels = c("Spatial", "FFPE"))
ST$source = factor(ST$source, levels = c("This_dataset", "10x Genomics"))
ST$slide = factor(ST$slide, levels = c("Sagittal","Coronal"))
 
p1 <- DimPlot(ST,label = F,pt.size = 0.5,label.size = 5,group.by = "sample") + scale_color_npg()
p1

p2 <- DimPlot(ST,label = F,pt.size = 0.5,label.size = 5,group.by = "type",
              cols = c("#00AFBB","#FC4E07"))
p2

p3 <- DimPlot(ST,label = F,pt.size = 0.5,label.size = 5,group.by = "source",
              cols = c("#EE0000","#3B4992"))
p3

p4 <- DimPlot(ST,label = F,pt.size = 0.5,label.size = 5,group.by = "slide",
              cols = c("#53A85F","#E95C59"))
p4

p1 | p2 | p3 | p4 # 4*20

DimPlot(ST,label = F,pt.size = 0.5,label.size = 5,split.by = "sample",group.by = "sample")

p5 <- DimPlot(subset(ST,subset = sample %in% c("Control_1", "Control_2","Control_3","Control_4")),
              label = F,pt.size = 0.5,label.size = 5,split.by = "sample",group.by = "sample")
p5

# 3.Find markers of each cluster  

# all_DEGS <- FindAllMarkers(ST,only.pos = T,logfc.threshold = 0.5,min.pct = 0.5)
# write.csv(all_DEGS,"all_DEGS.csv")

# 4.Marker expression level of each cluster

# DoHeatmap(ST, slot = "data",features = c(
#   "Nphs1","Nphs2","Podxl","Kdr","Eng", # Glom
#   "Slc17a3","Slc22a8","Fxyd2", # proximal tubule(PT)
#   "Slc5a2","Slc5a12", # Proximal convoluted tubule cell(PCT)
#   "Atp11a", "Slc13a3", # proximal straight tubules(PST)
#   "Slc27a2","Lrp2","Slc22a8","Slc5a2","Slc5a12","Fxyd2","Slc17a3", # proximal tubule(PT)
#   "Atp11a", "Slc13a3","Slc34a1","Gpx3",# Proximal convoluted tubule cell(PTC)
#   "Aqp1","Bst1", # Descending loop of Henle (DLH)
#   "Slc12a1","Umod","Cldn8","Krt18","Krt8",  # Ascending loop of Henle(ALH)
#   "Slc12a3"  # Distal convoluted tubule(DCT)
#   "Atp6v0d2","Atp6v1g3","Slc4a1","Aqp6","Slc26a4","Hmx2", # Collecting duct intercalated cell(CD-IC)
#   "Aqp2","Hsd11b2" #  Collecting duct principal / epithelial cell (CD-PC)
#   "Rhbg","Insrr","Stmn1" # Collecting duct transitional cell (CD-TC)
#   "Cdca3","Mki67" # Novel cell
#   "Kdr","Ehd3","Plat","Vim","S100a4","Aqp1","Bst1","Pecam1","Eng","Cd34"  # Endo # Endothelial(Endo)
#   "Nphs1", "Nphs2", # Podocyte (Podo)
#   "Vim","S100a4", # Pericytes and vascular smooth muscle (Peri)
#   "Plac8", # Fibroblast (Fibro)
#   "C1qa","C1qb" # Macrophage (Macro)
#   "Cd79a", "Cd79b" # B lymphocyte (B lymph)
#   "Cxcr6","Ltb","Il7r","Cd3d","Cd3e","Ifng"  # T lymphocyte (T lymph)
#   "Gzma","Nkg7","Gnly", # Natural killer cell (NK)
#   "Lyz2","Cd14" # Monocytes (Mono)
#   "S100a8","S100a9" # Neutrophil (Neutro)
#   "Col1a1","Col1a2","Tagln","Acta2","C3","Vim","Myl9","S100a9","Thbs1" # Inter
#   "Slc14a1","Slc14a2","Ly6d","Muc20","Psca","Upk1b","Upk3b" # Uro
# )) + NoLegend()

# SpatialFeaturePlot

p1 <- 
  (SpatialFeaturePlot(ST,features = c("Nphs1"),alpha = 1,images = c("Control_1")) + color2) + # Glom 
  (SpatialFeaturePlot(ST,features = c("Slc22a6"),alpha = 1,images = c("Control_1")) + color2) + # PT-S2
  (SpatialFeaturePlot(ST,features = c("Slc22a7"),alpha = 1,images = c("Control_1")) + color2) + # PT-S3
  (SpatialFeaturePlot(ST,features = c("Slc12a3"),alpha = 1,images = c("Control_1")) + color2) + # DCT
  (SpatialFeaturePlot(ST,features = c("Aqp1"),alpha = 1,images = c("Control_1")) + color2) +  # CD-PC
  (SpatialFeaturePlot(ST,features = c("Aqp2"),alpha = 1,images = c("Control_1")) + color2)  # CD-PC
p1

p2 <- 
  (SpatialFeaturePlot(ST,features = c("Nphs1"),alpha = 1,images = c("Control_2.1")) + color2) + # Glom 
  (SpatialFeaturePlot(ST,features = c("Slc22a6"),alpha = 1,images = c("Control_2.1")) + color2) + # PT-S2
  (SpatialFeaturePlot(ST,features = c("Slc22a7"),alpha = 1,images = c("Control_2.1")) + color2) + # PT-S3
  (SpatialFeaturePlot(ST,features = c("Slc12a3"),alpha = 1,images = c("Control_2.1")) + color2) + # DCT
  (SpatialFeaturePlot(ST,features = c("Aqp1"),alpha = 1,images = c("Control_2.1")) + color2) +  # CD-PC
  (SpatialFeaturePlot(ST,features = c("Aqp2"),alpha = 1,images = c("Control_2.1")) + color2)  # CD-PC
p2

p3 <- 
  (SpatialFeaturePlot(ST,features = c("Nphs1"),alpha = 1,images = c("Control_5.4")) + color2) + # Glom 
  (SpatialFeaturePlot(ST,features = c("Slc22a6"),alpha = 1,images = c("Control_5.4")) + color2) + # PT-S2
  (SpatialFeaturePlot(ST,features = c("Slc22a7"),alpha = 1,images = c("Control_5.4")) + color2) + # PT-S3
  (SpatialFeaturePlot(ST,features = c("Slc12a3"),alpha = 1,images = c("Control_5.4")) + color2) + # DCT
  (SpatialFeaturePlot(ST,features = c("Aqp1"),alpha = 1,images = c("Control_5.4")) + color2) +  # CD-PC
  (SpatialFeaturePlot(ST,features = c("Aqp2"),alpha = 1,images = c("Control_5.4")) + color2)  # CD-PC
p3

p4 <- 
  (SpatialFeaturePlot(ST,features = c("Nphs1"),alpha = 1,images = c("Control_6.5")) + color2) + # Glom 
  (SpatialFeaturePlot(ST,features = c("Slc22a6"),alpha = 1,images = c("Control_6.5")) + color2) + # PT-S2
  (SpatialFeaturePlot(ST,features = c("Slc22a7"),alpha = 1,images = c("Control_6.5")) + color2) + # PT-S3
  (SpatialFeaturePlot(ST,features = c("Slc12a3"),alpha = 1,images = c("Control_6.5")) + color2) + # DCT
  (SpatialFeaturePlot(ST,features = c("Aqp1"),alpha = 1,images = c("Control_6.5")) + color2) +  # CD-PC
  (SpatialFeaturePlot(ST,features = c("Aqp2"),alpha = 1,images = c("Control_6.5")) + color2)  # CD-PC
p4

p2 / p3 / p4 # 20*20

# 5. Rename each cluster

ST2 <- ST
new.cluster.ids <- c("ALH", 
                     "PT-S2", "DCT", "ALH","PT-S3","PT-S3",
                     "PT-S1", "DLH", "CD-TC","CD-PC","PT-S2",
                     "PT-S1","CD-IC","CD-PC","Adipo","Uro",
                     "PT-S2","PT-S1","Glom","Inter_Immune","CD-IC")
table(new.cluster.ids)
names(new.cluster.ids) <- levels(ST2)
ST2 <- RenameIdents(ST2, new.cluster.ids)

ST2$spot_type <- ST2@active.ident
ST2$spot_type <- factor(ST2$spot_type,levels = c("Glom","PT-S1","PT-S2","PT-S3",
                                                 "DLH","ALH","DCT","CD-IC","CD-PC","CD-TC",
                                                 "Inter_Immune","Uro","Adipo"))
# ST2$type <- factor(ST2$type,levels = c("Control","AKI","CKI"))

# 6. visualization of cluster annotation 

cols = c("Glom" = '#E63863',
         "PT-S1" = '#E4C755',
         "PT-S2" = '#E59CC4',
         "PT-S3" = '#00BFC4',
         # "PT-Inj" = '#E95C59',
         "DLH" = '#53A85F',
         "ALH" = '#F1BB72',
         "DCT" = '#F3B1A0',
         "CD-IC" = '#D6E7A3',
         "CD-PC" = '#57C3F3',
         "CD-TC" = '#AB3282',
         "Inter_Immune" = "#23452F",
         "Inter" = '#8C549C',
         "Adipo" = "grey",
         "Uro" = '#3C5488B2')

p4 <- DimPlot(ST2,pt.size = 0.5,label = F,label.size = 5,
        group.by = "spot_type",cols = cols)  # 6*5
p4 # 6*5

# 7. spatial visualization of cluster annotation 

ST2@images$Control_1
ST2@images$Control_2.1
ST2@images$Control_3.2
ST2@images$Control_4.3
ST2@images$Control_5.4 
ST2@images$Control_6.5

ST2@active.ident <- factor(ST2$spot_type)

SpatialDimPlot(ST2,images = c("Control_1","Control_2.1",
                               "Control_3.2","Control_4.3",
                               "Control_5.4","Control_6.5"),
                ncol = 3,cols = cols,alpha = 1) + NoLegend()
# 12*5

# 8.Marker expression level of each spot type

p5 <- DotPlot(ST2,group.by = "spot_type",features = rev(c(
  "Nphs1","Nphs2", # Glom
  "Lrp2",#"Slc34a1", # PT
  "Slc17a3", # proximal tubule(PT)
  "Slc5a2","Slc5a12", # PT-S1
  "Slc13a3","Slc22a6",# "Slc17a3", # PT-S2
  "Slc22a7","Slc22a13", # "Slc7a13","Bcat1",# PT-S3
  # "Havcr1",#"Lcn2","Vim", # PT-injured
  "Aqp1","Bst1", # Descending loop of Henle (DLH)
  "Slc12a1","Umod",#"Cldn8","Krt18","Krt8",  # Ascending loop of Henle(ALH)
  "Slc12a3",  # Distal convoluted tubule(DCT)
  "Atp6v0d2","Atp6v1g3",#"Slc4a1","Aqp6","Slc26a4","Hmx2", # Collecting duct intercalated cell(CD-IC)
  "Aqp2","Hsd11b2", #  Collecting duct principal / epithelial cell (CD-PC)
  # "Insrr",
  "Stmn1", # Collecting duct transitional cell (CD-TC)
  "C1qa", # Macro (Immune)
  "Cd3d", # T lymph (Immune)
  # "S100a8","S100a9", # Neutrophil (Neutro)
  "Col1a1","Col1a2", # Interstitium （Inter）
  "Ly6d","Upk1b", # Urothelium (Uro）
  "Adipoq","Lpl" # Adipocytes (Adipo)
    ))) + scale_color_gradientn(colours = viridis::viridis(20), guide = guide_colorbar(ticks.colour = "black",
    frame.colour = "black"), name = "Average \n expression") + 
  coord_flip() + xlab(NULL) + ylab(NULL) + 
  annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+
  annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5) + 
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) # 5*8
p5 # 7*6

# 9.Regional markers expression

color2 <- scale_fill_gradientn(colours = colorRampPalette(c("royalblue","white","firebrick1"))(10),
                               guide = guide_colorbar(ticks.colour = "black",frame.colour = "black"))

p1 <- (SpatialFeaturePlot(ST,features = c("Lrp2"),alpha = 1,images = c("Control_1")) + color2) /
  (SpatialFeaturePlot(ST,features = c("Aqp1"),alpha = 1,images = c("Control_1")) + color2) /
  (SpatialFeaturePlot(ST,features = c("Aqp2"),alpha = 1,images = c("Control_1")) + color2)

p2 <- (SpatialFeaturePlot(ST,features = c("Lrp2"),alpha = 1,images = c("Control_5.4")) + color2) /
  (SpatialFeaturePlot(ST,features = c("Aqp1"),alpha = 1,images = c("Control_5.4")) + color2) /
  (SpatialFeaturePlot(ST,features = c("Aqp2"),alpha = 1,images = c("Control_5.4")) + color2)

p3 <- (SpatialFeaturePlot(ST,features = c("Lrp2"),alpha = 1,images = c("Control_6.5")) + color2) /
  (SpatialFeaturePlot(ST,features = c("Aqp1"),alpha = 1,images = c("Control_6.5")) + color2) /
  (SpatialFeaturePlot(ST,features = c("Aqp2"),alpha = 1,images = c("Control_6.5")) + color2)

p1 | p2 | p3 # 10*10

# 10. Save annotated ST-seq dataset
ST <- ST2
saveRDS(ST,"ST_anno.rds")

# ST <- readRDS("ST_anno.rds")
# table(ST$spot_type)

# 11.loadng ST-seq and scRNA-seq datasets

scRNA = readRDS("scRNA_anno.rds")
table(scRNA$cell_type)

# 12.Add module scores to the ST-seq dataset

all_DEGS <- FindAllMarkers(scRNA,only.pos = T,logfc.threshold = 0.5,min.pct = 0.5)
write.csv(all_DEGS,"scRNA_all_DEGS.csv")
table(all_DEGS$cluster)

PT_genes <- list(all_DEGS[all_DEGS$cluster %in% "PT",]$gene)
DLH_genes <- list(all_DEGS[all_DEGS$cluster %in% "DLH",]$gene)
ALH_genes <- list(all_DEGS[all_DEGS$cluster %in% "ALH",]$gene)
DCT_genes <- list(all_DEGS[all_DEGS$cluster %in% "DCT",]$gene)
CD_IC_genes <- list(all_DEGS[all_DEGS$cluster %in% "CD-IC",]$gene)
CD_PC_genes <- list(all_DEGS[all_DEGS$cluster %in% "CD-PC",]$gene)
Endo_genes <- list(all_DEGS[all_DEGS$cluster %in% "Endo",]$gene)
Fibro_genes <- list(all_DEGS[all_DEGS$cluster %in% "Fibro",]$gene)
Podo_genes <- list(all_DEGS[all_DEGS$cluster %in% "Podo",]$gene)
Peri_genes <- list(all_DEGS[all_DEGS$cluster %in% "Peri",]$gene)
Macro_genes <- list(all_DEGS[all_DEGS$cluster %in% "Macro",]$gene)
Neutro_genes <- list(all_DEGS[all_DEGS$cluster %in% "Neutro",]$gene)
B_genes <- list(all_DEGS[all_DEGS$cluster %in% "B lymph",]$gene)
T_NK_genes <- list(all_DEGS[all_DEGS$cluster %in% "T lymph/NK",]$gene)

ST <- AddModuleScore(object = ST,features = c(PT_genes,DLH_genes,ALH_genes,DCT_genes,
                                              CD_IC_genes,CD_PC_genes,Endo_genes,Fibro_genes,
                                              Podo_genes,Peri_genes,Macro_genes,Neutro_genes,B_genes,T_NK_genes),
                     assay="SCT",
                      name=c("PT","DLH","ALH","DCT","CD-IC","CD-PC","Endo","Fibro",
                             "Podo","Peri","Macro","Neutro","B cell","T/NK"))

colnames(head(ST@meta.data))

p6 <- DotPlot(ST,group.by = "spot_type",
              features = rev(c("Podo9","Endo7","PT1","DLH2","ALH3","DCT4","CD.IC5","CD.PC6",
                               "Fibro8","Peri10","Macro11","Neutro12","B.cell13","T.NK14")))  +
  coord_flip() + xlab(NULL) + ylab(NULL) + scale_color_gradientn(colours = magma(20), 
                                                                 guide = guide_colorbar(ticks.colour = "black",frame.colour = "black"), name = "Module socre") +
  annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1) +
  annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5) + 
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) # 5*8
 
p6

p5 | p6 # 12*7

meta_module_score <- ST@meta.data
write.csv(meta_module_score,"meta_module_score.csv")

# 13. Spot numbers of each sample and spot type
data <- as.data.frame(table(ST$spot_type))
data <- data[order(data$Freq,decreasing=TRUE),] 
data$Var1 <- factor(data$Var1,levels=unique(data$Var1)) 
data

ggplot(data,aes(x=Var1,y=Freq)) + geom_bar(aes(fill=Var1),stat="identity")+
  geom_text(aes(label=Freq),vjust=-1,size=3) + theme_test()+
  labs(x="",y="Cell number of spot types") + guides(fill="none")+
  theme(text=element_text(size=16),axis.text.x=element_text(angle=45,hjust=1)) + 
  scale_fill_manual(values = cols)


data2 <- as.data.frame(table(ST$sample))
ggplot(data2,aes(x=Var1,y=Freq)) + geom_bar(aes(fill=Var1),stat="identity")+
  geom_text(aes(label=Freq),vjust=-1,size=3) + theme_test()+
  labs(x="",y="Spot number of sample") + guides(fill="none")+
  theme(text=element_text(size=16),axis.text.x=element_text(angle=45,hjust=1)) 

####

ST2 <- subset(ST,subset =  sample %in% c("Control_1", "Control_2","Control_3","Control_4"))
data <- as.data.frame(table(ST2$spot_type))
data <- data[order(data$Freq,decreasing=TRUE),] 
data$Var1 <- factor(data$Var1,levels=unique(data$Var1)) 
data

DimPlot(ST2)

ggplot(data,aes(x=Var1,y=Freq)) + geom_bar(aes(fill=Var1),stat="identity")+
  geom_text(aes(label=Freq),vjust=-1,size=3) + theme_test()+
  labs(x="",y="Cell number of spot types") + guides(fill="none")+
  theme(text=element_text(size=16),axis.text.x=element_text(angle=45,hjust=1)) + 
  scale_fill_manual(values = cols)

dfsam <- as.data.frame(table(ST2$spot_type))

ggplot(dfsam, aes(x = "", y = Freq,fill = Var1)) + 
  geom_bar(stat = "identity")+
  coord_polar(theta = "y") + 
  theme_bw() + 
  labs(x = "", y = "", title = "")+
  theme(axis.ticks = element_blank())+
  theme(axis.text.x = element_blank()) + 
  theme(panel.grid=element_blank()) +    ## 去掉白色圆框和中间的坐标线
  theme(panel.border=element_blank())+
  scale_fill_manual(values = cols,labels = myLabel) # 6x6

myLabel = as.vector(dfsam$Var1)
myLabel = paste(myLabel, "(", round(dfsam$Freq / sum(dfsam$Freq) * 100, 2), "%)"
                , sep = "")

# 14. Spot proportion change

ST$slice <- factor(ST$sample,levels = rev(c("Control_1","Control_2","Control_3","Control_4","Control_5","Control_6")))
# ST$spot_type <- factor(ST$spot_type,levels = rev(levels(ST$spot_type)))

data <- data.frame(table(ST$spot_type,ST$slice))
head(data)
#  
table(ST$sample)
# Control_1 Control_2 Control_3 Control_4 Control_5 Control_6 
# 2768      3259      3099      3305      3124      1438

data$proportion <- c(data[1:14,"Freq"]/2768,data[15:28,"Freq"]/3259,
                     data[29:42,"Freq"]/3099,data[43:56,"Freq"]/3305,
                     data[57:70,"Freq"]/3124,data[71:84,"Freq"]/1438) * 100
data
 
colnames(data) <- c("spot_type","sample","Freq","proportion")

ggplot(data = data,aes(x = sample,y=Freq,fill = spot_type)) +
  geom_bar(stat="identity",position = "fill",width = 0.7)+
  scale_y_continuous(labels = scales::percent)+
  scale_fill_manual(values = cols, name = "" ) +
  theme_classic() + 
  labs(y = 'Fraction of different type',x="") + 
  coord_flip()  +
  # NoLegend() + 
  annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1) +
  annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

# 15. MIA

sc.markers <- FindAllMarkers(scRNA,only.pos = T,logfc.threshold = 0.25,min.pct = 0.25)
sc.markers$d = sc.markers$pct.1 - sc.markers$pct.2
sc.main.marker = subset(sc.markers , avg_log2FC > 0.25 & p_val_adj < 0.05 & d > 0.2)
sc.main.marker = sc.main.marker %>% arrange(cluster,desc(avg_log2FC))
sc.main.marker = as.data.frame(sc.main.marker)
sc.main.marker$cluster = paste('sc',sc.main.marker$cluster,sep = '_')
head(sc.main.marker)
table(sc.main.marker$cluster)
write.csv(sc.main.marker,"sc.main.marker.csv")
sc.main.marker <- read.csv("../3.stRNA/2.Integration/sc.main.marker.csv",row.names = 1)
  
region_marker <- FindAllMarkers(ST,only.pos = T,logfc.threshold = 0.25,min.pct = 0.25)
region_marker$d <- region_marker$pct.1 - region_marker$pct.2
region_main_marker <-  subset(region_marker,avg_log2FC > 0.25 & p_val_adj < 0.05 & d > 0.1)
region_main_marker <-  region_main_marker %>% arrange(cluster,desc(avg_log2FC))
region_main_marker <-  as.data.frame(region_main_marker)
region_main_marker$cluster = paste('spatial',region_main_marker$cluster,sep = '_')
head(region_main_marker)
table(region_main_marker$cluster)
write.csv(region_main_marker,"region.main_marker.csv")

celltype_specific = sc.main.marker[,c("cluster","gene")]
colnames(celltype_specific)[1]="celltype"

region_specific = region_main_marker[,c("cluster","gene")]
colnames(region_specific)[1]="region"
head(region_specific)

source("../3.stRNA/MIA.R")

Result = zhao_MIA(region_specific,celltype_specific,"Cis","../3.stRNA/2.Integration/")

MIA_data <- read.csv("../3.stRNA/2.Integration/plot.data.csv",row.names = 1)
pheatmap::pheatmap(t(MIA_data),
                   cluster_rows = F,cluster_cols = F,
                   border_color = "black",
                   main = "cellular module score (MIA)", 
                   angle_col = 45,
                   cellwidth = 15,cellheight = 15,
                   color = colorRampPalette(c("royalblue","white","firebrick1"))(10),
                   # gaps_row = c(5,10,15,20,25,30,35,40,45,50,55,60,65),
                   scale = "row",
                   show_rownames = T,
                   treeheight_row = 10,treeheight_col = 10) 
