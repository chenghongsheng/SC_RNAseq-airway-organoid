library(dplyr)
library(Seurat)
library(patchwork)
library(SeuratWrappers)
library(ggplot2)
library(egg)
library(tidyr)

SCBO1.data<-Read10X(data.dir="C:/Users/LKC/OneDrive - Nanyang Technological University/Project/Archives/Lung organoids/Rawdat/SCBO1_filtered_feature_bc_matrix")
SCBO2.data<-Read10X(data.dir="C:/Users/LKC/OneDrive - Nanyang Technological University/Project/Archives/Lung organoids/Rawdat/SCBO2_filtered_feature_bc_matrix")
SCNPO1.data<-Read10X(data.dir="C:/Users/LKC/OneDrive - Nanyang Technological University/Project/Archives/Lung organoids/Rawdat/SCNPO1_filtered_feature_bc_matrix")
SCNPO2.data<-Read10X(data.dir="C:/Users/LKC/OneDrive - Nanyang Technological University/Project/Archives/Lung organoids/Rawdat/SCNPO2_filtered_feature_bc_matrix")

SCBO1<-CreateSeuratObject(counts=SCBO1.data,project="SCBO1",min.cells = 3,min.features = 200)
SCBO2<-CreateSeuratObject(counts=SCBO2.data,project="SCBO2",min.cells = 3,min.features = 200)
SCNPO1<-CreateSeuratObject(counts=SCNPO1.data,project="SCNPO1",min.cells = 3,min.features = 200)
SCNPO2<-CreateSeuratObject(counts=SCNPO2.data,project="SCNPO2",min.cells = 3,min.features = 200)

SCBO1[["percent.mt"]]<-PercentageFeatureSet(SCBO1,pattern="^MT-")
SCBO2[["percent.mt"]]<-PercentageFeatureSet(SCBO2,pattern="^MT-")
SCNPO1[["percent.mt"]]<-PercentageFeatureSet(SCNPO1,pattern="^MT-")
SCNPO2[["percent.mt"]]<-PercentageFeatureSet(SCNPO2,pattern="^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(SCBO1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol=3)
VlnPlot(SCBO2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(SCNPO1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(SCNPO2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

tiff('SCBO1_QC.tiff',width=2000,height=2000,units='px',res=300,compression='lzw')
VlnPlot(SCBO1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol=3)
dev.off()

tiff('SCBO2_QC.tiff',width=2000,height=2000,units='px',res=300,compression='lzw')
VlnPlot(SCBO2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol=3)
dev.off()

tiff('SCNPO1_QC.tiff',width=2000,height=2000,units='px',res=300,compression='lzw')
VlnPlot(SCNPO1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol=3)
dev.off()

tiff('SCNPO2_QC.tiff',width=2000,height=2000,units='px',res=300,compression='lzw')
VlnPlot(SCNPO2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol=3)
dev.off()

#filter away cells with high mitochondrial gene percentage (>20%)
SCBO1<-subset(SCBO1,subset = percent.mt<20)
SCBO2<-subset(SCBO2,subset = percent.mt<20)
SCNPO1<-subset(SCNPO1,subset = percent.mt<20)
SCNPO2<-subset(SCNPO2,subset = percent.mt<20)

#normalization
lung<-merge(SCNPO1,y=c(SCNPO2,SCBO1,SCBO2),add.cell.ids = c("SCNPO1",'SCNPO2','SCBO1','SCBO2'),project = "LungOrganoids")
lung$orig.ident<-factor(lung$orig.ident,levels=c("SCNPO1",'SCNPO2','SCBO1','SCBO2'))
lung.list<-SplitObject(lung,split.by = "orig.ident")
lung.list <- lapply(X = lung.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = lung.list, nfeatures = 3000)
lung.list <- PrepSCTIntegration(object.list = lung.list, anchor.features = features)
lung.anchors <- FindIntegrationAnchors(object.list = lung.list, normalization.method = "SCT",
                                       anchor.features = features)
lung.combined.sct <- IntegrateData(anchorset = lung.anchors, normalization.method = "SCT")

DefaultAssay(lung.combined.sct) <- "integrated"

lung.combined.sct <- RunPCA(lung.combined.sct, verbose = FALSE)
DimHeatmap(lung.combined.sct, dims = 1:20, cells = 500, balanced = TRUE)
ElbowPlot(lung.combined.sct,ndims=20)
lung.combined.sct<- FindNeighbors(lung.combined.sct, dims = 1:20)
lung.combined.sct <- FindClusters(lung.combined.sct, resolution = 0.3)
lung.combined.sct <- RunUMAP(lung.combined.sct, reduction = "pca", dims = 1:20)

saveRDS(lung.combined.sct,file='lung.combined.sct.RDS')

#collapse replicates from healthy and COPD samples
lung.combined.sct<-readRDS("./lung.combined.sct.RDS")

lung.combined.sct$status<-as.character(lung.combined.sct$orig.ident)

lung.combined.sct$status<-replace(lung.combined.sct$status, lung.combined.sct$status=="SCNPO1","Healthy")
lung.combined.sct$status<-replace(lung.combined.sct$status, lung.combined.sct$status=="SCBO1","Healthy")
lung.combined.sct$status<-replace(lung.combined.sct$status, lung.combined.sct$status=="SCNPO2","COPD")
lung.combined.sct$status<-replace(lung.combined.sct$status, lung.combined.sct$status=="SCBO2","COPD")

lung.combined.sct$status<-factor(lung.combined.sct$status,levels=c("Healthy","COPD"))

lung.combined.sct <- FindClusters(lung.combined.sct, resolution = 0.3)
DefaultAssay(lung.combined.sct)<-"integrated"
Idents(lung.combined.sct)<-lung.combined.sct$integrated_snn_res.0.3

dev.new()
tiff('UMAP_res0.3.tiff',width=2000,height=2000,units='px',res=300,compression='lzw')
DimPlot(lung.combined.sct, reduction = "umap", label = TRUE,repel = TRUE)
dev.off()

dev.new()
tiff('UMAP_res0.3_split.tiff',width=2000,height=4000,units='px',res=300,compression='lzw')
DimPlot(lung.combined.sct, reduction = "umap", label = TRUE,repel = TRUE,split.by = "status",ncol = 1)
dev.off()

DefaultAssay(lung.combined.sct)<-"RNA"
lung.combined.marker_res.0.3 <- FindAllMarkers(lung.combined.sct, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
lung.combined.allmarker_res0.3<-FindAllMarkers(lung.combined.sct,only.pos = F,min.pct = 0.25,logfc.threshold = 0.25)

lung.combined.marker_res.0.3 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10

write.csv(lung.combined.marker_res.0.3,file="lung.combined.marker_res.0.3.csv")
write.csv(lung.combined.allmarker_res0.3,file="lung.combined.allmarker_res0.3.csv")
write.csv(top10,file="top10_res0.3.csv")

dev.new()
tiff('HM-cluster_res0.3.tiff',width=3000,height=3000,units='px',res=300,compression='lzw')
DoHeatmap(subset(lung.combined.sct,downsample=100), features = top10$gene)
dev.off()

saveRDS(lung.combined.sct,file='lung.combined.sct.RDS')

#cell type annotation
library(celldex)
library(SingleR)
library(pheatmap)
library(multtest)
library(metap)

ref_BlueEncode <- BlueprintEncodeData()
ref_HPCA <-HumanPrimaryCellAtlasData()

lung.sce<-lung.combined.sct
DefaultAssay(lung.sce)<-"RNA"

lung.sce[['SCT']]<-NULL
lung.sce[['integrated']]<-NULL
lung.sce<-as.SingleCellExperiment(lung.sce)

pred_BlueEncode <- SingleR(test=lung.sce, ref=ref_BlueEncode, labels=ref_BlueEncode$label.main)
table(pred_BlueEncode$labels)
plotScoreHeatmap(pred_BlueEncode)

pred_HPCA <- SingleR(test=lung.sce, ref=ref_HPCA, labels=ref_HPCA$label.main)
table(pred_HPCA$labels)
plotScoreHeatmap(pred_HPCA)

tab_pred_BlueEncode <- table(Assigned=pred_BlueEncode$pruned.labels, Cluster=lung.combined.sct$seurat_clusters) #resolution=0.1
tab_pred_HPCA <- table(Assigned=pred_HPCA$pruned.labels, Cluster=lung.combined.sct$seurat_clusters) #resolution=0.1

tiff('HM-BlueEncode.tiff',width=3000,height=3000,units='px',res=300,compression='lzw')
plotScoreHeatmap(pred_BlueEncode)
dev.off()

tiff('HM-HPCA.tiff',width=3000,height=3000,units='px',res=300,compression='lzw')
plotScoreHeatmap(pred_HPCA)
dev.off()

write.csv(tab_pred_BlueEncode,"tab_BlueEncode.csv")
write.csv(tab_pred_HPCA,"tab_HPCA.csv")
#BlueprintEncode & HPCA does not provide good annotation

#obtain marker genes of each cluster
C0<-lung.combined.allmarker_res0.3[which(lung.combined.allmarker_res0.3$cluster==0),]
C1<-lung.combined.allmarker_res0.3[which(lung.combined.allmarker_res0.3$cluster==1),]
C2<-lung.combined.allmarker_res0.3[which(lung.combined.allmarker_res0.3$cluster==2),]
C3<-lung.combined.allmarker_res0.3[which(lung.combined.allmarker_res0.3$cluster==3),]
C4<-lung.combined.allmarker_res0.3[which(lung.combined.allmarker_res0.3$cluster==4),]
C5<-lung.combined.allmarker_res0.3[which(lung.combined.allmarker_res0.3$cluster==5),]
C6<-lung.combined.allmarker_res0.3[which(lung.combined.allmarker_res0.3$cluster==6),]
C7<-lung.combined.allmarker_res0.3[which(lung.combined.allmarker_res0.3$cluster==7),]
C8<-lung.combined.allmarker_res0.3[which(lung.combined.allmarker_res0.3$cluster==8),]
C9<-lung.combined.allmarker_res0.3[which(lung.combined.allmarker_res0.3$cluster==9),]

rownames(C0)<-C0$gene
C0<-C0[,c(2),drop=F]
C0<-C0[order(C0$avg_log2FC,decreasing = T),,drop=F]

rownames(C1)<-C1$gene
C1<-C1[,2,drop=F]
C1<-C1[order(C1$avg_log2FC,decreasing = T),,drop=F]

rownames(C2)<-C2$gene
C2<-C2[,2,drop=F]
C2<-C2[order(C2$avg_log2FC,decreasing = T),,drop=F]

rownames(C3)<-C3$gene
C3<-C3[,2,drop=F]
C3<-C3[order(C3$avg_log2FC,decreasing = T),,drop=F]

rownames(C4)<-C4$gene
C4<-C4[,2,drop=F]
C4<-C4[order(C4$avg_log2FC,decreasing = T),,drop=F]

rownames(C5)<-C5$gene
C5<-C5[,2,drop=F]
C5<-C5[order(C5$avg_log2FC,decreasing = T),,drop=F]

rownames(C6)<-C6$gene
C6<-C6[,2,drop=F]
C6<-C6[order(C6$avg_log2FC,decreasing = T),,drop=F]

rownames(C7)<-C7$gene
C7<-C7[,2,drop=F]
C7<-C7[order(C7$avg_log2FC,decreasing = T),,drop=F]

rownames(C8)<-C8$gene
C8<-C8[,2,drop=F]
C8<-C8[order(C8$avg_log2FC,decreasing = T),,drop=F]

rownames(C9)<-C9$gene
C9<-C9[,2,drop=F]
C9<-C9[order(C9$avg_log2FC,decreasing = T),,drop=F]

#export for GSEA cell annotation
write.csv(C0,file="Cluster 0.csv")
write.csv(C1,file="Cluster 1.csv")
write.csv(C2,file="Cluster 2.csv")
write.csv(C3,file="Cluster 3.csv")
write.csv(C4,file="Cluster 4.csv")
write.csv(C5,file="Cluster 5.csv")
write.csv(C6,file="Cluster 6.csv")
write.csv(C7,file="Cluster 7.csv")
write.csv(C8,file="Cluster 8.csv")
write.csv(C9,file="Cluster 9.csv")

#newID
Idents(lung.combined.sct)<-lung.combined.sct$ClusterNames_0.3
lung.combined.sct$ClusterNames_0.3<-lung.combined.sct$integrated_snn_res.0.3
new.ID<-c("Basal cells 2",'Club cells 1',"Basal cells 1","Cycling basal cells","Club cells 2",'Goblet cells 1','Goblet cells 2','Basal cells 1','Ciliated epithelial cells','Epithelial cells')
names(new.ID) <- levels(lung.combined.sct)
lung.combined.sct <- RenameIdents(lung.combined.sct, new.ID)
lung.combined.sct$annotated_res0.3<-lung.combined.sct@active.ident

lung.combined.sct$annotated_res0.3<-factor(lung.combined.sct$annotated_res0.3,levels=c("Epithelial cells","Cycling basal cells","Basal cells 1",'Basal cells 2',
                                                                                       "Club cells 1",'Club cells 2',"Goblet cells 1",'Goblet cells 2','Ciliated epithelial cells'))

Idents(lung.combined.sct)<-lung.combined.sct$annotated_res0.3

#change the color of UMAP clusters
library(RColorBrewer)
col<-c('#680956',"#F03B20","#F07620","#DA8F00","#FFFA2F","#DAD400","#18B146","#007122","#273FA3")
pie(rep(1,9),col = col)

UMAP<-DimPlot(lung.combined.sct, reduction = "umap", label = F,repel = TRUE,pt.size = 0.8)
UMAP<-UMAP+scale_color_manual(values = col)

dev.new()
tiff('UMAP_res0.3_unlabel.tiff',width=2500,height=2000,units='px',res=300,compression='lzw')
UMAP
dev.off()

UMAP_split<-DimPlot(lung.combined.sct, reduction = "umap", label = F,repel = TRUE,split.by = "status",pt.size = 1)
UMAP_split<-UMAP_split+scale_color_manual(values = col)

dev.new()
tiff('UMAP_res0.3_split_unlabel.tiff',width=4000,height=2000,units='px',res=300,compression='lzw')
UMAP_split
dev.off()

DefaultAssay(lung.combined.sct)<-"integrated"
organoid.feature <- c('MKI67','TP63',"KRT5",'KRT13','KRT4','SCGB1A1','MUC5AC','FOXJ1')
combined_dot<-DotPlot(lung.combined.sct, features = organoid.feature,scale = T,dot.scale = 8)+RotatedAxis() 
combined_dot<-combined_dot+scale_color_gradient(low='lightgrey',high = 'red',limits=c(-1.5,2.5),breaks=c(-1,0,1,2))
combined_dot

dev.new()
tiff('combined_dotplot.tiff',width=2500,height=2000,units='px',res=300,compression='lzw')
combined_dot
dev.off()

#identify percentage of each subpopulation
lung.combined.sct_split<-SplitObject(lung.combined.sct,split.by = "status")

Healthy_count<-as.data.frame(summary(lung.combined.sct_split$Healthy@active.ident))
COPD_count<-as.data.frame(summary(lung.combined.sct_split$COPD@active.ident))

colnames(Healthy_count)<-c("Healthy")
colnames(COPD_count)<-c("COPD")

Healthy_count$Healthy_percentage<-c(round(Healthy_count$Healthy*100/sum(Healthy_count$Healthy),2))
COPD_count$COPD_percentage<-c(round(COPD_count$COPD*100/sum(COPD_count$COPD),2))

comb_count<-cbind(Healthy_count,COPD_count)

comb_count$cluster<-rownames(comb_count)
comb_count<-comb_count[,c(2,4,5)]
colnames(comb_count)<-c("Healthy","COPD",'cluster')
comb_count$cluster<-factor(comb_count$cluster,levels=c("Epithelial cells","Cycling basal cells","Basal cells 1",'Basal cells 2',
                                                       "Club cells 1",'Club cells 2',"Goblet cells 1",'Goblet cells 2','Ciliated epithelial cells'))
comb_count<-comb_count%>%gather(Group,Percentage,1:2)

comb_count$Group<-factor(comb_count$Group,levels=c("Healthy","COPD"))

write.csv(comb_count,file='comb_count.csv')

#barplot
abundance<-ggplot(comb_count,aes(fill=cluster,y=Percentage,x=Group))+
  geom_bar(position=position_fill(reverse = T), stat = 'identity',size=0.1)+
  scale_fill_manual(values=c(col))+
  guides(fill=guide_legend(reverse=T))+
  theme(text = element_text(size=15),axis.text.x = element_text(angle=45,hjust=1),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line=element_line(size=1,color="black"))+
  xlab("Groups")+
  ylab("Abundance (%)")

tiff('abundance.tiff',width=2000,height=2000,units='px',res=300,compression='lzw')
abundance
dev.off()

saveRDS(lung.combined.sct,file='lung.combined.sct.RDS')

#collapse trajectory
lung.combined.sct<-readRDS("./lung.combined.sct.RDS")

lung.combined.sct_split<-SplitObject(lung.combined.sct,split.by = "status")
Healthy.cds <- as.cell_data_set(lung.combined.sct_split$Healthy)
COPD.cds <- as.cell_data_set(lung.combined.sct_split$COPD)
lung.cds<-as.cell_data_set(lung.combined.sct)

Healthy.cds<-cluster_cells(cds = Healthy.cds,reduction_method = "UMAP")
COPD.cds<-cluster_cells(cds = COPD.cds,reduction_method = "UMAP")
lung.cds<-cluster_cells(cds=lung.cds,reduction_method = "UMAP")

Healthy.cds <- learn_graph(Healthy.cds, use_partition = F,close_loop = F)
COPD.cds <- learn_graph(COPD.cds, use_partition = F,close_loop = F,learn_graph_control = list(minimal_branch_len=8))
lung.cds<-learn_graph(lung.cds, use_partition=F, close_loop=F)

Healthy.cds<-order_cells(Healthy.cds)
COPD.cds<-order_cells(COPD.cds)
lung.cds<-order_cells(lung.cds)

Healthy.cds<-order_cells(Healthy.cds,root_pr_nodes = c('Y_40','Y_70','Y_24'))
COPD.cds<-order_cells(COPD.cds,root_pr_nodes = c('Y_20','Y_26','Y_112'))
lung.cds<-order_cells(lung.cds,root_pr_nodes = c('Y_85','Y_15','Y_106','Y_49'))


lung.Pseudotime<-plot_cells(lung.cds,
                            color_cells_by = "pseudotime",
                            trajectory_graph_segment_size = 1,
                            label_cell_groups=F,
                            label_leaves=F,
                            label_roots=F,
                            label_branch_points=F,
                            graph_label_size=3, label_principal_points = F)

Healthy.Pseudotime<-plot_cells(Healthy.cds,
                               color_cells_by = "pseudotime",
                               trajectory_graph_segment_size = 1,
                               label_cell_groups=F,
                               label_leaves=F,
                               label_roots=F,
                               label_branch_points=F,
                               graph_label_size=3, label_principal_points = F)


COPD.Pseudotime<-plot_cells(COPD.cds,
                            color_cells_by = "pseudotime",
                            trajectory_graph_segment_size = 1,
                            label_cell_groups=F,
                            label_leaves=F,
                            label_roots =F,
                            label_branch_points=F,
                            graph_label_size=3, label_principal_points = F)
dev.new()
tiff('Healthy.pseudotime.tiff',width=2500,height=2000,units='px',res=300,compression='lzw')
Healthy.Pseudotime
dev.off()

dev.new()
tiff('COPD.pseudotime.tiff',width=2500,height=2000,units='px',res=300,compression='lzw')
COPD.Pseudotime
dev.off()

dev.new()
tiff('lung.pseudotime.tiff',width=2500,height=2000,units='px',res=300,compression='lzw')
lung.Pseudotime
dev.off()

lung.combined.sct_split$Healthy <- AddMetaData(
  object = lung.combined.sct_split$Healthy,
  metadata = Healthy.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "Pseudotime")

lung.combined.sct_split$COPD<- AddMetaData(
  object = lung.combined.sct_split$COPD,
  metadata = COPD.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "Pseudotime")

lung.combined.sct<- AddMetaData(
  object = lung.combined.sct,
  metadata = lung.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "Pseudotime")

FeaturePlot(lung.combined.sct_split$Healthy, c("Pseudotime"), pt.size = 0.1) & scale_color_viridis(option="plasma")
FeaturePlot(lung.combined.sct_split$COPD, c("Pseudotime"), pt.size = 0.1) & scale_color_viridis(option="plasma")
dev.new()
FeaturePlot(lung.combined.sct, c("Pseudotime"), pt.size = 0.1) & scale_color_viridis(option="plasma")

saveRDS(lung.combined.sct,file='lung.combined.sct.RDS')


library('gplots')
library('RColorBrewer')
library("genefilter")
library('pheatmap')
library('viridis')
library('circlize')

DefaultAssay(lung.combined.sct_split$Healthy)<-"integrated"
lung_healthy_marker<-FindAllMarkers(lung.combined.sct_split$Healthy, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

lung_healthy_marker %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC) -> top50_healthy
write.csv(top50_healthy,"top50.healthy.csv")

HM_genemarker <- read.csv("C:/Users/LKC/OneDrive - Nanyang Technological University/Project/Lung organoids/HM_genemarker.csv")

#Healthy 
Healthy.sub<-subset(lung.combined.sct_split$Healthy,features=HM_genemarker$Gene)

Healthy.sub<-subset(Healthy.sub,idents=c("Cycling basal cells","Basal cells 1",
                                         'Basal cells 2',"Club cells 1",'Club cells 2',
                                         "Goblet cells 1",'Goblet cells 2','Ciliated epithelial cells'))
Healthy.mat<-as.matrix(Healthy.sub@assays$integrated@scale.data)

Healthy_col<-as.data.frame(cbind(Healthy.sub@meta.data[,12,drop=F],Idents(Healthy.sub)))
colnames(Healthy_col)<-c("Pseudotime","ID")

Healthy_col<-Healthy_col[!c(Healthy_col$ID=="Cycling basal cells" & Healthy_col$Pseudotime>20),]
Healthy_col<-Healthy_col[!c(Healthy_col$ID=="Basal cells 2" & Healthy_col$Pseudotime<20),]


Healthy_col<-Healthy_col[order(Healthy_col$ID,Healthy_col$Pseudotime,decreasing = F),]

Healthy.mat<-Healthy.mat[HM_genemarker$Gene,rownames(Healthy_col)]

my_color_annotation<-list(ID = c("Cycling basal cells"="#F03B20","Basal cells 1"="#F07620",
                                 'Basal cells 2'="#DA8F00","Club cells 1"="#FFFA2F",'Club cells 2'="#DAD400",
                                 "Goblet cells 1"="#18B146",'Goblet cells 2'="#007122",'Ciliated epithelial cells'="#273FA3"),
                          Pseudotime = viridis(3,option="plasma"), status=c("Healthy"="#FFD9A7","COPD"="#F0AD54"))

hm_color<- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

HM_Healthy<-pheatmap(Healthy.mat, scale = 'row', color = hm_color,show_rownames = T,
                     cluster_rows = F, cluster_cols = F,
                     annotation_col =Healthy_col,annotation_colors = my_color_annotation,
                     clustering_distance_rows = "correlation",
                     clustering_distance_cols= "correlation",
                     labels_col = F, show_colnames = F,gaps_row = c(5,10,15,20))

tiff('HM_Healthy.tiff',width=2500,height=2000,units='px',res=300,compression='lzw')
HM_Healthy
dev.off()

#COPD 
COPD.sub<-subset(lung.combined.sct_split$COPD,features=HM_genemarker$Gene)

COPD.sub<-subset(COPD.sub,idents=c("Cycling basal cells","Basal cells 1",
                                   'Basal cells 2',"Club cells 1",'Club cells 2',
                                   "Goblet cells 1",'Goblet cells 2','Ciliated epithelial cells'))
COPD.mat<-as.matrix(COPD.sub@assays$integrated@scale.data)

COPD_col<-as.data.frame(cbind(COPD.sub@meta.data[,12,drop=F],Idents(COPD.sub)))
colnames(COPD_col)<-c("Pseudotime","ID")

COPD_col<-COPD_col[order(COPD_col$ID,COPD_col$Pseudotime,decreasing = F),]

COPD.mat<-COPD.mat[HM_genemarker$Gene,rownames(COPD_col)]

HM_COPD<-pheatmap(COPD.mat, scale = 'row', color = hm_color,show_rownames = T,
                  cluster_rows = F, cluster_cols = F,
                  annotation_col =COPD_col,annotation_colors = my_color_annotation,
                  clustering_distance_rows = "correlation",
                  clustering_distance_cols= "correlation",
                  labels_col = F, show_colnames = F,gaps_row = c(5,10,15,20))

tiff('HM_COPD.tiff',width=2500,height=2000,units='px',res=300,compression='lzw')
HM_COPD
dev.off()

#combine Healthy & COPD
lung.sub<-subset(lung.combined.sct,features=HM_genemarker$Gene)

lung.sub<-subset(lung.sub,idents=c("Cycling basal cells","Basal cells 1",
                                   'Basal cells 2',"Club cells 1",'Club cells 2',
                                   "Goblet cells 1",'Goblet cells 2','Ciliated epithelial cells'))


lung.mat<-as.matrix(lung.sub@assays$SCT@scale.data)

lung_col<-as.data.frame(lung.sub@meta.data[,c(7,11,12)])
colnames(lung_col)<-c('status',"ID","Pseudotime")
lung_col<-lung_col[,c(2:3,1)]
lung_col$ID<-factor(lung_col$ID,levels=c("Cycling basal cells","Basal cells 1",
                                         'Basal cells 2',"Club cells 1",'Club cells 2',
                                         "Goblet cells 1",'Goblet cells 2','Ciliated epithelial cells'))

lung_col<-lung_col[order(lung_col$status,lung_col$ID,lung_col$Pseudotime,decreasing = F),]
lung_col<-lung_col[c(8403:1,8404:19601),]

set.seed(1234)
size=1000
healthy_col<-lung_col[sample(rownames(lung_col[c(1:8403),]),size=1000),]
copd_col<-lung_col[sample(rownames(lung_col[c(8404:19601),]),size=1000),]
comb_col<-as.data.frame(rbind(healthy_col,copd_col))

comb_col<-comb_col[order(comb_col$status,comb_col$ID,comb_col$Pseudotime,decreasing = F),]
comb_col<-comb_col[c(1000:1,1001:2000),]

lung.mat<-lung.mat[HM_genemarker$Gene,rownames(comb_col)]

mybreak = seq(-3,3,length.out = 101)
lung.mat_scaled <- t(scale(t(lung.mat)))

HM_lung<-pheatmap(lung.mat_scaled, scale = 'none', color = hm_color,show_rownames = T,
                  cluster_rows = F, cluster_cols = F,breaks = mybreak,cellheight = 10,
                  annotation_col =comb_col,annotation_colors = my_color_annotation,
                  clustering_distance_rows = "correlation",
                  clustering_distance_cols= "correlation",
                  labels_col = F, show_colnames = F,gaps_row = c(5,10,15,20),gaps_col=c(1000))
dev.new()
tiff('HM_lung.tiff',width=2500,height=2000,units='px',res=300,compression='lzw')
HM_lung
dev.off()

ggplot(data = df, aes(x = session, y = score, group = Group, shape = Group)) +
  geom_line(aes(linetype = Group), size = 1, stat="summary", fun.y=mean) + 
  geom_point(size = 5, fill = "white", stat="summary", fun.y=mean)

DefaultAssay(lung.sub)<-"SCT"
lung.sub_split<-SplitObject(lung.sub,split.by = "status")

plain <- function(x,...) {
  format(x, ..., scientific = FALSE, trim = TRUE)
}

Healthy_SCGB1A1<-VlnPlot(lung.sub_split$Healthy,features = c("SCGB1A1"),log=T,pt.size=0) 
Healthy_SCGB1A1<-Healthy_SCGB1A1+
  geom_point(size = 5, fill = "black", stat="summary", fun=median)+
  scale_fill_manual(values=c("#F03B20","#F07620","#DA8F00","#FFFA2F","#DAD400","#18B146","#007122","#273FA3"))+
  scale_y_log10(limits=c(1,10),labels=plain)

COPD_SCGB1A1<-VlnPlot(lung.sub_split$COPD,features = c("SCGB1A1"),log=T,pt.size=0) 
COPD_SCGB1A1<-COPD_SCGB1A1+
  geom_point(size = 5, fill = "black", stat="summary", fun=median)+
  scale_fill_manual(values=c("#F03B20","#F07620","#DA8F00","#FFFA2F","#DAD400","#18B146","#007122","#273FA3"))+
  scale_y_log10(limits=c(1,10),labels=plain)

Healthy_MUC5AC<-VlnPlot(lung.sub_split$Healthy,features = c("MUC5AC"),log=T,pt.size=0) 
Healthy_MUC5AC<-Healthy_MUC5AC+
  stat_summary(fun = median, geom='point', size = 5, colour = "black")+
  scale_fill_manual(values=c("#F03B20","#F07620","#DA8F00","#FFFA2F","#DAD400","#18B146","#007122","#273FA3"))+
  scale_y_log10(limits=c(1,10),labels=plain)

COPD_MUC5AC<-VlnPlot(lung.sub_split$COPD,features = c("MUC5AC"),log=T,pt.size=0) 
COPD_MUC5AC<-COPD_MUC5AC+
  stat_summary(fun = median, geom='point', size = 5, colour = "black")+
  scale_fill_manual(values=c("#F03B20","#F07620","#DA8F00","#FFFA2F","#DAD400","#18B146","#007122","#273FA3"))+
  scale_y_log10(limits=c(1,10),labels=plain)

Healthy_TP63<-VlnPlot(lung.sub_split$Healthy,features = c("TP63"),log=T,pt.size=0) 
Healthy_TP63<-Healthy_TP63+
  stat_summary(fun = median, geom='point', size = 5, colour = "black")+
  scale_fill_manual(values=c("#F03B20","#F07620","#DA8F00","#FFFA2F","#DAD400","#18B146","#007122","#273FA3"))+
  scale_y_log10(limits=c(1,10),labels=plain)

COPD_TP63<-VlnPlot(lung.sub_split$COPD,features = c("TP63"),log=T,pt.size=0) 
COPD_TP63<-COPD_TP63+
  stat_summary(fun = median, geom='point', size = 5, colour = "black")+
  scale_fill_manual(values=c("#F03B20","#F07620","#DA8F00","#FFFA2F","#DAD400","#18B146","#007122","#273FA3"))+
  scale_y_log10(limits=c(1,10),labels=plain)

Healthy_FOXJ1<-VlnPlot(lung.sub_split$Healthy,features = c("FOXJ1"),log=T,pt.size=0) 
Healthy_FOXJ1<-Healthy_FOXJ1+
  stat_summary(fun = median, geom='point', size = 5, colour = "black")+
  scale_fill_manual(values=c("#F03B20","#F07620","#DA8F00","#FFFA2F","#DAD400","#18B146","#007122","#273FA3"))+
  scale_y_log10(limits=c(1,10),labels=plain)

COPD_FOXJ1<-VlnPlot(lung.sub_split$COPD,features = c("FOXJ1"),log=T,pt.size=0) 
COPD_FOXJ1<-COPD_FOXJ1+
  stat_summary(fun = median, geom='point', size = 5, colour = "black")+
  scale_fill_manual(values=c("#F03B20","#F07620","#DA8F00","#FFFA2F","#DAD400","#18B146","#007122","#273FA3"))+
  scale_y_log10(limits=c(1,10),labels=plain)


tiff('Healthy_SCGB1A1.tiff',width=2000,height=2000,units='px',res=300,compression='lzw')
Healthy_SCGB1A1
dev.off()

tiff('COPD_SCGB1A1.tiff',width=2000,height=2000,units='px',res=300,compression='lzw')
COPD_SCGB1A1
dev.off()

tiff('Healthy_MUC5AC.tiff',width=2000,height=2000,units='px',res=300,compression='lzw')
Healthy_MUC5AC
dev.off()

tiff('COPD_MUC5AC.tiff',width=2000,height=2000,units='px',res=300,compression='lzw')
COPD_MUC5AC
dev.off()

tiff('Healthy_TP63.tiff',width=2000,height=2000,units='px',res=300,compression='lzw')
Healthy_TP63
dev.off()

tiff('COPD_TP63.tiff',width=2000,height=2000,units='px',res=300,compression='lzw')
COPD_TP63
dev.off()

tiff('Healthy_FOXJ1.tiff',width=2000,height=2000,units='px',res=300,compression='lzw')
Healthy_FOXJ1
dev.off()

tiff('COPD_FOXJ1.tiff',width=2000,height=2000,units='px',res=300,compression='lzw')
COPD_FOXJ1
dev.off()

#get DEGs

lung.combined.sct$cluster.group<-paste(Idents(lung.combined.sct),lung.combined.sct$status,sep="_")

lung.combined.sct$cluster.group<-factor(lung.combined.sct$cluster.group,
                                        level=c("Epithelial cells_Healthy","Cycling basal cells_Healthy","Basal cells 1_Healthy",'Basal cells 2_Healthy',
                                                "Club cells 1_Healthy",'Club cells 2_Healthy',"Goblet cells 1_Healthy",'Goblet cells 2_Healthy','Ciliated epithelial cells_Healthy',
                                                "Epithelial cells_COPD","Cycling basal cells_COPD","Basal cells 1_COPD",'Basal cells 2_COPD',
                                                "Club cells 1_COPD",'Club cells 2_COPD',"Goblet cells 1_COPD",'Goblet cells 2_COPD','Ciliated epithelial cells_COPD'))

Idents(lung.combined.sct)<-"cluster.group"

DefaultAssay(lung.combined.sct)<-"RNA"
CBC_COPDvH <- FindMarkers(lung.combined.sct, ident.1 = "Cycling basal cells_COPD", ident.2 = "Cycling basal cells_Healthy", verbose = T)

BC_COPDvH <- FindMarkers(lung.combined.sct, ident.1 = c("Basal cells 1_COPD","Basal cells 2_COPD"), 
                         ident.2 = c("Basal cells 1_Healthy","Basal cells 2_Healthy"), verbose = T)

Club_COPDvH <- FindMarkers(lung.combined.sct, ident.1 = c("Club cells 1_COPD","Club cells 2_COPD"), 
                           ident.2 = c("Club cells 1_Healthy","Club cells 2_Healthy"), verbose = T)

Goblet_COPDvH <- FindMarkers(lung.combined.sct, ident.1 = c("Goblet cells 1_COPD","Goblet cells 2_COPD"), 
                             ident.2 = c("Goblet cells 1_Healthy","Goblet cells 2_Healthy"), verbose = T)

Cilia_COPDvH <- FindMarkers(lung.combined.sct, ident.1 = "Ciliated epithelial cells_COPD", ident.2 = "Ciliated epithelial cells_Healthy", verbose = T)

Idents(lung.combined.sct)<-"status"
COPDvH <- FindMarkers(lung.combined.sct, ident.1 = "COPD", ident.2 = "Healthy", verbose = T)

CBC_COPDvH.sig<-CBC_COPDvH[(CBC_COPDvH$p_val_adj<0.05 & abs(CBC_COPDvH$avg_log2FC)>1),]
BC_COPDvH.sig<-BC_COPDvH[(BC_COPDvH$p_val_adj<0.05 & abs(BC_COPDvH$avg_log2FC)>1),]
Club_COPDvH.sig<-Club_COPDvH[(Club_COPDvH$p_val_adj<0.05 & abs(Club_COPDvH$avg_log2FC)>1),]
Goblet_COPDvH.sig<-Goblet_COPDvH[(Goblet_COPDvH$p_val_adj<0.05 & abs(Goblet_COPDvH$avg_log2FC)>1),]
Cilia_COPDvH.sig<-Cilia_COPDvH[(Cilia_COPDvH$p_val_adj<0.05 & abs(Cilia_COPDvH$avg_log2FC)>1),]
COPDvH.sig<-COPDvH[(COPDvH$p_val_adj<0.05 & abs(COPDvH$avg_log2FC)>1),]

write.csv(CBC_COPDvH,file="CBC_COPDvH.csv")
write.csv(BC_COPDvH,file="BC_COPDvH.csv")
write.csv(Club_COPDvH,file="Club_COPDvH.csv")
write.csv(Goblet_COPDvH,file="Goblet_COPDvH.csv")
write.csv(Cilia_COPDvH,file="Cilia_COPDvH.csv")
write.csv(COPDvH,file="COPDvH.csv")
write.csv(CBC_COPDvH.sig,file="CBC_COPDvH.sig.csv")
write.csv(BC_COPDvH.sig,file="BC_COPDvH.sig.csv")
write.csv(Club_COPDvH.sig,file="Club_COPDvH.sig.csv")
write.csv(Goblet_COPDvH.sig,file="Goblet_COPDvH.sig.csv")
write.csv(Cilia_COPDvH.sig,file="Cilia_COPDvH.sig.csv")
write.csv(COPDvH.sig,file="COPDvH.sig.csv")

saveRDS(lung.combined.sct,file='lung.combined.sct.RDS')