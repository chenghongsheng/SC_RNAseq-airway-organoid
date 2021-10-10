library(GEOquery)
library(Biobase)
library(tidyverse)
library(DESeq2)
library(ggplot2)

setwd("C:/Users/LKC/OneDrive - Nanyang Technological University/Project/Lung organoids/RNAseq")

#GSE124180
GSE124180 <- getGEO('GSE124180')
supp_GSE124180 = getGEOSuppFiles('GSE124180')
fnames_GSE124180 = rownames(supp_GSE124180)

GSE124180_mat<-read.delim(fnames_GSE124180[1],header=TRUE,stringsAsFactors = F)
GSE124180_mat<-GSE124180_mat[-c(1),]
rownames(GSE124180_mat)<-GSE124180_mat$ENSEMBL_GENEID
GSE124180_mat$ENSEMBL_GENEID<-NULL
GSE124180_mat$X<-NULL

GSE124180_coldat2<-GSE124180$GSE124180_series_matrix.txt.gz@phenoData@data
GSE124180_coldat<-GSE124180_coldat[GSE124180_coldat$source_name_ch1=="large airway",]
GSE124180_mat<-GSE124180_mat[,GSE124180_coldat$title]

i<-c(1:21)

GSE124180_mat[ , i] <- apply(GSE124180_mat[ , i], 2,            # Specify own function within apply
                             function(x) as.numeric(as.character(x)))
sapply(GSE124180_mat, class)

#GSE162154
GSE162154 <- getGEO('GSE162154')
GSE162154_coldat <- GSE162154$GSE162154_series_matrix.txt.gz@phenoData@data
GSE162154_PRJNA680768_mat <- read.delim("C:/Users/LKC/OneDrive - Nanyang Technological University/Project/Lung organoids/RNAseq/PRJNA680768.txt", comment.char="#")
GSE162154_PRJNA680768_mat<-GSE162154_PRJNA680768_mat[,-c(2:6)]
rownames(GSE162154_PRJNA680768_mat)<-GSE162154_PRJNA680768_mat$Geneid
GSE162154_PRJNA680768_mat$Geneid<-NULL

GSE162154_PRJNA680768_coldat <- read.csv("C:/Users/LKC/OneDrive - Nanyang Technological University/Project/Lung organoids/RNAseq/PRJNA680768_coldat.csv")
GSE162154_PRJNA680768_coldat$Group<-GSE162154_coldat$`sample type:ch1`

#GSE146532
GSE146532 <- getGEO('GSE146532')
GSE146532_coldat<-GSE146532$GSE146532_series_matrix.txt.gz@phenoData@data
GSE146532_PRJNA610819_mat <- read.delim("C:/Users/LKC/OneDrive - Nanyang Technological University/Project/Lung organoids/RNAseq/PRJNA610819.txt", comment.char="#")
GSE146532_PRJNA610819_mat<-GSE146532_PRJNA610819_mat[,-c(2:6)]
rownames(GSE146532_PRJNA610819_mat)<-GSE146532_PRJNA610819_mat$Geneid
GSE146532_PRJNA610819_mat$Geneid<-NULL

GSE146532_PRJNA610819_coldat <- read.csv("C:/Users/LKC/OneDrive - Nanyang Technological University/Project/Lung organoids/RNAseq/PRJNA610819_coldat.csv")
GSE146532_PRJNA610819_coldat<-GSE146532_PRJNA610819_coldat[c(1:40),]

GSE146532_coldat<-GSE146532_coldat[GSE146532_coldat$`treatment:ch1`=="Uninfected" & !GSE146532_coldat$`diagnosis:ch1`=="Asthma",]

GSE146532_dds<-DESeqDataSetFromMatrix(countData = GSE146532_PRJNA610819_mat, colData = GSE146532_PRJNA610819_coldat,design = ~Group)

GSE146532.coll_dds<-collapseReplicates(GSE146532_dds,GSE146532_dds$Subject)
GSE146532_coll_mat<-assay(GSE146532.coll_dds)

GSE146532.coll_coldat<-as.data.frame(colData(GSE146532.coll_dds))

#Ivan's data
Ivan_mat <- read.csv("C:/Users/LKC/OneDrive - Nanyang Technological University/Project/Lung organoids/RNAseq/Ivan_FC.csv", header = T)
Ivan_mat$X<-NULL
rownames(Ivan_mat)<-Ivan_mat$ensembl_gene_id
Ivan_mat$ensembl_gene_id<-NULL
Ivan_mat<-Ivan_mat[,c(1:6)]

Ivan_coldat <- read.csv("C:/Users/LKC/OneDrive - Nanyang Technological University/Project/Lung organoids/RNAseq/Ivan_FC_coldat.csv")
rownames(Ivan_coldat)<-Ivan_coldat$SampleID
Ivan_coldat$SampleID<-NULL
Ivan_coldat<-Ivan_coldat[,1,drop=F]


#common geneID
common<-intersect(rownames(Ivan_mat),rownames(GSE124180_mat))
common<-intersect(common,rownames(GSE162154_PRJNA680768_mat))
common<-intersect(common,rownames(GSE146532_PRJNA610819_mat))

#combined matrix
Ivan_mat<-Ivan_mat[common,]
GSE124180_mat<-GSE124180_mat[common,]
GSE162154_PRJNA680768_mat<-GSE162154_PRJNA680768_mat[common,]
GSE146532_PRJNA610819_mat<-GSE146532_PRJNA610819_mat[common,]

comb_mat<-cbind(Ivan_mat,GSE124180_mat)
comb_mat<-cbind(comb_mat,GSE162154_PRJNA680768_mat)
comb_mat<-cbind(comb_mat,GSE146532_coll_mat)

#combined coldata
Ivan_coldat$Group<-c(rep("COPD",3),rep("Healthy",3))
Ivan_coldat$batch<-c(rep("Ivan",6))
Ivan_coldat$Sample_Details<-NULL

GSE124180_coldat<-GSE124180_coldat[,c(1,66)]
GSE124180_coldat$Group<-ifelse(GSE124180_coldat$`copd:ch1`=='case',"COPD","Healthy")
GSE124180_coldat$batch<-rep("GSE124180",nrow(GSE124180_coldat))
rownames(GSE124180_coldat)<-GSE124180_coldat$title
GSE124180_coldat<-GSE124180_coldat[,c(3:4)]

GSE162154_PRJNA680768_coldat$Group<-ifelse(GSE162154_PRJNA680768_coldat$Group=="COPD","COPD","Healthy")
GSE162154_PRJNA680768_coldat$batch<-rep('GSE162154',nrow(GSE162154_PRJNA680768_coldat))
rownames(GSE162154_PRJNA680768_coldat)<-GSE162154_PRJNA680768_coldat$ID
GSE162154_PRJNA680768_coldat$ID<-NULL
GSE162154_PRJNA680768_coldat<-GSE162154_PRJNA680768_coldat[,c(3:4)]

GSE146532.coll_coldat$batch<-rep("GSE146532",nrow(GSE146532.coll_coldat))
GSE146532.coll_coldat<-GSE146532.coll_coldat[,c(5,7)]

comb_coldat<-rbind(Ivan_coldat,GSE124180_coldat)
comb_coldat<-rbind(comb_coldat,GSE162154_PRJNA680768_coldat)
comb_coldat<-rbind(comb_coldat,GSE146532.coll_coldat)

comb_coldat$Group<-factor(comb_coldat$Group,levels=c("Healthy",'COPD'))
comb_coldat$batch<-factor(comb_coldat$batch)

#Deseq2
library(DESeq2)
dds<-DESeqDataSetFromMatrix(countData = comb_mat,colData = comb_coldat,design=~batch+Group)
dds<-DESeq(dds)
vsd<-varianceStabilizingTransformation(dds, blind=TRUE)

plotPCA(vsd,intgroup="Group")
plotPCA(vsd,intgroup="batch")

#remove batch effect
mat <- assay(vsd)
mm <- model.matrix(~Group, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$batch, design=mm)
assay(vsd) <- mat
PCA1<-plotPCA(vsd,intgroup="Group")
PCA2<-plotPCA(vsd,intgroup="batch")
PCA1+PCA2

vsdtable<-as.data.frame(assay(vsd))

#DEG
#obtain Deseq2 results
res<-results(dds, contrast =c('Group',"COPD","Healthy"),alpha = 0.05 )
res<-res[order(res$padj),]
summary(res)
##batch correction does not resolve the inter-study differences. Proceed to independent analysis

res$symbols<-genesymbols[rownames(res)]
sig.res<-as.matrix(res[which(res$padj<0.05),])

dev.new()
tiff('PCA_group.tiff',width=4000,height=2000,units='px',res=300,compression='lzw')
PCA1+PCA2
dev.off()

write.csv(sig.res,file="sig.res.csv")

#GSE124180
GSE124180_mat<-read.delim(fnames_GSE124180[1],header=TRUE,stringsAsFactors = F)
GSE124180_mat<-GSE124180_mat[-c(1),]
rownames(GSE124180_mat)<-GSE124180_mat$ENSEMBL_GENEID
GSE124180_mat$ENSEMBL_GENEID<-NULL
GSE124180_mat$X<-NULL
GSE124180_mat<-GSE124180_mat[,rownames(GSE124180_coldat)]

i<-c(1:21)

GSE124180_mat[ , i] <- apply(GSE124180_mat[ , i], 2,            # Specify own function within apply
                             function(x) as.numeric(as.character(x)))
sapply(GSE124180_mat, class)

GSE124180_coldat$Group<-factor(GSE124180_coldat$Group,levels=c("Healthy","COPD"))

dds_GSE124180<-DESeqDataSetFromMatrix(countData = GSE124180_mat,colData = GSE124180_coldat,design=~Group)
dds_GSE124180<-DESeq(dds_GSE124180)
vsd_GSE124180<-varianceStabilizingTransformation(dds_GSE124180, blind=TRUE)

plotPCA(vsd_GSE124180,intgroup="Group")+geom_text(aes(label=name),vjust=-0.1,col="black",size=3)

res_GSE124180<-results(dds_GSE124180, contrast =c('Group',"COPD","Healthy"),alpha = 0.05 )
res_GSE124180<-res_GSE124180[order(res_GSE124180$padj),]
summary(res_GSE124180)

#out of 36917 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 9, 0.024%
#LFC < 0 (down)     : 6, 0.016%
#outliers [1]       : 58, 0.16%
#low counts [2]     : 10472, 28%
#(mean count < 0)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results

#GSE162154
GSE162154_PRJNA680768_mat <- read.delim("C:/Users/LKC/OneDrive - Nanyang Technological University/Project/Lung organoids/RNAseq/PRJNA680768.txt", comment.char="#")
GSE162154_PRJNA680768_mat<-GSE162154_PRJNA680768_mat[,-c(2:6)]
rownames(GSE162154_PRJNA680768_mat)<-GSE162154_PRJNA680768_mat$Geneid
GSE162154_PRJNA680768_mat$Geneid<-NULL

GSE162154_PRJNA680768_coldat <- read.csv("C:/Users/LKC/OneDrive - Nanyang Technological University/Project/Lung organoids/RNAseq/PRJNA680768_coldat.csv")
GSE162154_PRJNA680768_coldat$Group<-factor(GSE162154_PRJNA680768_coldat$Group,levels=c("Healthy",'COPD'))

dds_GSE162154<-DESeqDataSetFromMatrix(countData = GSE162154_PRJNA680768_mat,colData = GSE162154_PRJNA680768_coldat ,design=~Group)
dds_GSE162154<-DESeq(dds_GSE162154)
vsd_GSE162154<-varianceStabilizingTransformation(dds_GSE162154, blind=TRUE)

plotPCA(vsd_GSE162154,intgroup="Group")+geom_text(aes(label=name),vjust=-0.1,col="black",size=3) #good separation of clusters based on disease status, selected for further comparison with sc-RNAseq data

res_GSE162154<-results(dds_GSE162154, contrast =c('Group',"COPD","Healthy"),alpha = 0.05 )
res_GSE162154<-res_GSE162154[order(res_GSE162154$padj),]
summary(res_GSE162154)

res_GSE162154<-as.data.frame(res_GSE162154)
res_GSE162154$symbols<-genesymbols[rownames(res_GSE162154)]

sig.res_GSE162154<-res_GSE162154[which((res_GSE162154$padj<0.05)& (abs(res_GSE162154$log2FoldChange)>1)),]

write.csv(res_GSE162154,"res_GSE162154.csv")
write.csv(sig.res_GSE162154,"sig.res_GSE162154.csv")

#GSE146532
GSE146532_coldat<-GSE146532$GSE146532_series_matrix.txt.gz@phenoData@data
GSE146532_PRJNA610819_mat <- read.delim("C:/Users/LKC/OneDrive - Nanyang Technological University/Project/Lung organoids/RNAseq/PRJNA610819.txt", comment.char="#")
GSE146532_PRJNA610819_mat<-GSE146532_PRJNA610819_mat[,-c(2:6)]
rownames(GSE146532_PRJNA610819_mat)<-GSE146532_PRJNA610819_mat$Geneid
GSE146532_PRJNA610819_mat$Geneid<-NULL

GSE146532_PRJNA610819_coldat <- read.csv("C:/Users/LKC/OneDrive - Nanyang Technological University/Project/Lung organoids/RNAseq/PRJNA610819_coldat.csv")
GSE146532_PRJNA610819_coldat<-GSE146532_PRJNA610819_coldat[c(1:40),]

GSE146532_coldat<-GSE146532_coldat[GSE146532_coldat$`treatment:ch1`=="Uninfected" & !GSE146532_coldat$`diagnosis:ch1`=="Asthma",]

GSE146532_dds<-DESeqDataSetFromMatrix(countData = GSE146532_PRJNA610819_mat, colData = GSE146532_PRJNA610819_coldat,design = ~Group)

GSE146532.coll_dds<-collapseReplicates(GSE146532_dds,GSE146532_dds$Subject)
GSE146532_coll_mat<-assay(GSE146532.coll_dds)

GSE146532.coll_coldat<-as.data.frame(colData(GSE146532.coll_dds))
GSE146532.coll_coldat$Group<-factor(GSE146532.coll_coldat$Group,levels = c("Healthy","COPD"))

dds_GSE146532<-DESeqDataSetFromMatrix(countData = GSE146532_coll_mat,colData = GSE146532.coll_coldat ,design=~Group)
dds_GSE146532<-DESeq(dds_GSE146532)
vsd_GSE146532<-varianceStabilizingTransformation(dds_GSE146532, blind=TRUE)

plotPCA(vsd_GSE146532,intgroup="Group")+geom_text(aes(label=name),vjust=-0.1,col="black",size=3)

res_GSE146532<-results(dds_GSE146532, contrast =c('Group',"COPD","Healthy"),alpha = 0.05 )
res_GSE146532<-res_GSE146532[order(res_GSE146532$padj),]
summary(res_GSE146532)

#In-house data
Ivan_mat <- read.csv("C:/Users/LKC/OneDrive - Nanyang Technological University/Project/Lung organoids/RNAseq/Ivan_FC.csv", header = T)
Ivan_mat$X<-NULL
rownames(Ivan_mat)<-Ivan_mat$ensembl_gene_id
Ivan_mat$ensembl_gene_id<-NULL
Ivan_mat<-Ivan_mat[,c(1:6)]

Ivan_coldat$Group <-factor(Ivan_coldat$Group,levels=c("Healthy","COPD")) 

dds_Ivan<-DESeqDataSetFromMatrix(countData = Ivan_mat,colData = Ivan_coldat ,design=~Group)
dds_Ivan<-DESeq(dds_Ivan)
vsd_Ivan<-varianceStabilizingTransformation(dds_Ivan, blind=TRUE)

plotPCA(vsd_Ivan,intgroup="Group")+geom_text(aes(label=name),vjust=-0.1,col="black",size=3)

res_Ivan<-results(dds_Ivan, contrast =c('Group',"COPD","Healthy"),alpha = 0.05 )
res_Ivan<-res_Ivan[order(res_Ivan$padj),]
summary(res_Ivan)

res_Ivan<-as.data.frame(res_Ivan)
res_Ivan$symbols<-genesymbols[rownames(res_Ivan)]

sig.res_Ivan<-res_Ivan[which((res_Ivan$padj<0.05)& (abs(res_Ivan$log2FoldChange)>1)),]

write.csv(res_Ivan,"res_Ivan.csv")
write.csv(sig.res_Ivan,"sig.res_Ivan.csv")

#map ensembl id to gene symbols
library(biomaRt)

Hsa.dataset<-useMart(dataset='hsapiens_gene_ensembl',biomart=("ensembl"))
query.gene<-union(rownames(Ivan_mat),rownames(GSE124180_mat))
query.gene<-union(query.gene,rownames(GSE162154_PRJNA680768_mat))
query.gene<-union(query.gene,rownames(GSE146532_PRJNA610819_mat))

Genemap<-getBM(attributes = c('ensembl_gene_id','hgnc_symbol'), 
               filters="ensembl_gene_id",
               values=query.gene,mart=Hsa.dataset)

genesymbols <- tapply(Genemap$hgnc_symbol, 
                      Genemap$ensembl_gene_id, paste, collapse="; ")


library(VennDiagram)
vennplot<-venn.diagram(x=list(rownames(sig.res_GSE162154),rownames(sig.res_Ivan)),euler.d = F,scaled=F,filename=NULL,
                       category.names = c("GSE162154","Ivan"), fill=c("#F8766D","#A3A500"), 
                       alpha=c(0.5,0.5), height=3000, width=3000, resolution=300,
                       imagetype='tiff',units='px', lwd='3', cex='1',
                       cat.cex='0.5', cat.fontface='bold',cat.pos=c(10,350))
grid.draw(vennplot)

#heatmap of comparison
library('gplots')
library('RColorBrewer')
library("genefilter")
library('pheatmap')
library(RColorBrewer)

HM_comparison_mat <- read.delim("C:/Users/LKC/OneDrive - Nanyang Technological University/Project/Lung organoids/RNAseq/HM_comparison.txt")

rownames(HM_comparison_mat)<-HM_comparison_mat$Entity.Name
HM_comparison_mat$Entity.Name<-NULL
HM_comparison_mat<-HM_comparison_mat[,c(1,3,7,4,2,5,6)]
colnames(HM_comparison_mat)<-c("GSE162154",'All cells','Basal cells','Club cells','Goblet cells','Ciliated epithelial cells','Cycling basal cells')


annotation_row <- read.csv("C:/Users/LKC/OneDrive - Nanyang Technological University/Project/Lung organoids/RNAseq/annotation_row.csv", row.names=1)
annotation_row<-annotation_row[c(!annotation_row$Function=="DNA damage and repair"),,drop=F]
annotation_row$Function<-factor(annotation_row$Function)


data_type<-c("GEO dataset",rep("Lung organoid SCRNAseq",6))
annotation_col <- as.data.frame(data_type, drop=F)
rownames(annotation_col)<-colnames(HM_comparison_mat)
annotation_col$data_type<-factor(annotation_col$data_type)
annotation_col<-annotation_col[c(1,3:7),,drop=F]

HM_comparison_mat<-HM_comparison_mat[rownames(annotation_row),rownames(annotation_col)]

hm_color<- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

HM_comparison<-pheatmap(HM_comparison_mat, scale = 'row', color = hm_color,show_rownames = F,
                        cluster_rows = T, cluster_cols = F,border_color = NA,
                        annotation_row = annotation_row,
                        annotation_col = annotation_col,
                        clustering_distance_rows = "correlation",
                        show_colnames = T)

HM_comparison_mat<-HM_comparison_mat[HM_comparison$tree_row$order,]

HM_comparison<-pheatmap(HM_comparison_mat, scale = 'row', color = hm_color,show_rownames = F,
                        cluster_rows =F, cluster_cols = F,border_color = NA,
                        annotation_row = annotation_row,
                        annotation_col = annotation_col,
                        clustering_distance_rows = "correlation",
                        show_colnames = T)

dev.new()
tiff('HM_comparison.tiff',width=1800,height=1800,units='px',res=300,compression='lzw')
HM_comparison
dev.off()
