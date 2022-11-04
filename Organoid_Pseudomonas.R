library(DESeq2)

counts <- read.csv("./counts.csv",check.names = F,row.names = 1)

coldata <- read.csv("./bulkRNAseq_metadata.csv",check.names = F,row.names = 1)
coldata$Treatment<-factor(coldata$Treatment,levels = c("HC_ctrl","HC_Pseudo", "COPD_ctrl","COPD_Pseudo"))

colnames(counts)==rownames(coldata)  ####sanity check

dds<-DESeqDataSetFromMatrix(countData=counts,colData = coldata, design=~Treatment)
dds<-DESeq(dds)

vsd<-varianceStabilizingTransformation(dds, blind=TRUE)

###PCA plot####
pca<-plotPCA(vsd,intgroup="Treatment",returnData=T)

library(ggplot2)

percentVar<-round(100*attr(pca,"percentVar"))

pca<-ggplot(pca, aes(PC1, PC2, color=Treatment)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme(aspect.ratio = 1)+
  coord_fixed()

tiff('PCA.tiff',width=2000,height=2000,units='px',res=300,compression='lzw')
pca
dev.off()

##### DESEq results #########
res_COPDvsHC_ctrl <-as.data.frame(results(dds,contrast = c("Treatment",'COPD_ctrl','HC_ctrl')))
sig_COPDvsHC_ctrl<-res_COPDvsHC_ctrl[which((res_COPDvsHC_ctrl$padj<0.05) & (abs(res_COPDvsHC_ctrl$log2FoldChange)>1)),]

res_COPDvsHC_Pseudo<-as.data.frame(results(dds,contrast = c("Treatment",'COPD_Pseudo','HC_Pseudo')))
sig_COPDvsHC_Pseudo<-res_COPDvsHC_Pseudo[which((res_COPDvsHC_Pseudo$padj<0.05) & (abs(res_COPDvsHC_Pseudo$log2FoldChange)>1)),]

res_PseudovsCtrl_HC<-as.data.frame(results(dds,contrast = c("Treatment",'HC_Pseudo','HC_ctrl')))
sig_PseudovsCtrl_HC<-res_PseudovsCtrl_HC[which((res_PseudovsCtrl_HC$padj<0.05) & (abs(res_PseudovsCtrl_HC$log2FoldChange)>1)),]

res_PseudovsCtrl_COPD<-as.data.frame(results(dds,contrast = c("Treatment",'COPD_Pseudo','COPD_ctrl')))
sig_PseudovsCtrl_COPD<-res_PseudovsCtrl_COPD[which((res_PseudovsCtrl_COPD$padj<0.05) & (abs(res_PseudovsCtrl_COPD$log2FoldChange)>1)),]

allDEG<-unique(c(rownames(sig_COPDvsHC_ctrl),
         rownames(sig_COPDvsHC_Pseudo),
         rownames(sig_PseudovsCtrl_COPD),
         rownames(sig_PseudovsCtrl_HC)))

##### add gene information #######
library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ENSEMBL_MART_ENSEMBL",host="www.ensembl.org"))
#mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

gene.info<-getBM(
  filters= "ensembl_gene_id",
  attributes= c("ensembl_gene_id", "entrezgene_id", "entrezgene_accession", "hgnc_symbol", "description"),
  values= allDEG,
  mart= mart, useCache = FALSE
)

gene.info<- gene.info[-which(duplicated(gene.info$ensembl_gene_id)), ]
rownames(gene.info) <- gene.info$ensembl_gene_id

sig_COPDvsHC_ctrl <- cbind(sig_COPDvsHC_ctrl, gene.info[rownames(sig_COPDvsHC_ctrl), ])
sig_COPDvsHC_Pseudo <- cbind(sig_COPDvsHC_Pseudo, gene.info[rownames(sig_COPDvsHC_Pseudo), ])
sig_PseudovsCtrl_COPD <- cbind(sig_PseudovsCtrl_COPD, gene.info[rownames(sig_PseudovsCtrl_COPD), ])
sig_PseudovsCtrl_HC <- cbind(sig_PseudovsCtrl_HC, gene.info[rownames(sig_PseudovsCtrl_HC), ])

write.csv(sig_COPDvsHC_ctrl,file="sig_COPDvsHC_ctrl.csv")
write.csv(sig_COPDvsHC_Pseudo,file="sig_COPDvsHC_Pseudo.csv")
write.csv(sig_PseudovsCtrl_COPD,file="sig_PseudovsCtrl_COPD.csv")
write.csv(sig_PseudovsCtrl_HC,file="sig_PseudovsCtrl_HC.csv")


##### heatmap ######
library('RColorBrewer')
library('pheatmap')

vsdtable<-as.data.frame(assay(vsd))
mat<-vsdtable
mat<-vsdtable[allDEG,]
mat<-mat[,c(4:6,10:12,1:3,7:9)]
mat<-as.matrix(mat)

anno_col<-coldata[,2,drop=F]
annotation_colour<-list(Treatment=c("HC_ctrl"="#F8766D",'HC_Pseudo'='#7CAE00',
                                    "COPD_ctrl"="#00BFC4",'COPD_Pseudo'='#C77CFF'))
  


hm_color<- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

hm<-pheatmap(mat, scale = 'row', color = hm_color,show_rownames = F,
         cluster_rows = T, cluster_cols =F,
         annotation_col=anno_col,annotation_colors = annotation_colour,
         clustering_distance_rows = "correlation",
         clustering_distance_cols= "correlation",
         labels_row = F,labels_col = F,show_colnames = F)


tiff('heatmap.tiff',width=3000,height=3000,units='px',res=300,compression='lzw')
hm
dev.off()

vsdmat<-mat[hm$tree_row$order,]
write.csv(vsdmat,file="vsdmat.csv")

vsdmat<-read.csv("C:/Users/hscheng/OneDrive - Nanyang Technological University/Project/Archives/Lung organoids/Lung_fig/Rebuttal/Pseudomonas RNAseq/vsdmat.csv", 
                 row.names=1)
clus1<-rownames(vsdmat[vsdmat$Clus==1,])
clus2<-rownames(vsdmat[vsdmat$Clus==2,])
clus3<-rownames(vsdmat[vsdmat$Clus==3,])
clus4<-rownames(vsdmat[vsdmat$Clus==4,])
clus5<-rownames(vsdmat[vsdmat$Clus==5,])

write.table(clus1,"clus1.txt")
write.table(clus2,"clus2.txt")
write.table(clus3,"clus3.txt")
write.table(clus4,"clus4.txt")
write.table(clus5,"clus5.txt")


library(ViSEAGO)
expressed_genes<-rownames(GSEA_COPDvsHC_ctrl)
write.table(expressed_genes,'expressed_genes.txt')

background<-scan("expressed_genes.txt",
                 quiet=TRUE,
                 what=""
)

selection1<-scan(
  "clus1.txt",
  quiet=TRUE,
  what=""
)

selection2<-scan(
  "clus2.txt",
  quiet=TRUE,
  what=""
)
selection3<-scan(
  "clus3.txt",
  quiet=TRUE,
  what=""
)
selection4<-scan(
  "clus4.txt",
  quiet=TRUE,
  what=""
)
selection5<-scan(
  "clus5.txt",
  quiet=TRUE,
  what=""
)


Ensembl<-ViSEAGO::Ensembl2GO()
ViSEAGO::available_organisms(Ensembl)

myGENE2GO<-ViSEAGO::annotate(
  "hsapiens_gene_ensembl",
  Ensembl
)
#############cluster 1############
Clus1_GO<-ViSEAGO::create_topGOdata(
  geneSel=selection1,
  allGenes=background,
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=5
)

classic1<-topGO::runTest(
  Clus1_GO,
  algorithm ="classic",
  statistic = "fisher"
)

Clus1_sResults<-ViSEAGO::merge_enrich_terms(
  Input=list(
    condition=c("Clus1_GO","classic1")
  )
)

Clus1_GO<-as.data.frame(Clus1_sResults@data)
write.table(Clus1_GO,"Clus1_GO.txt")

#############cluster 2############
Clus2_GO<-ViSEAGO::create_topGOdata(
  geneSel=selection2,
  allGenes=background,
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=5
)

classic2<-topGO::runTest(
  Clus2_GO,
  algorithm ="classic",
  statistic = "fisher"
)

Clus2_sResults<-ViSEAGO::merge_enrich_terms(
  Input=list(
    condition=c("Clus2_GO","classic2")
  )
)

Clus2_GO<-as.data.frame(Clus2_sResults@data)
write.table(Clus2_GO,"Clus2_GO.txt")

#############cluster 3############
Clus3_GO<-ViSEAGO::create_topGOdata(
  geneSel=selection3,
  allGenes=background,
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=5
)

classic3<-topGO::runTest(
  Clus3_GO,
  algorithm ="classic",
  statistic = "fisher"
)

Clus3_sResults<-ViSEAGO::merge_enrich_terms(
  Input=list(
    condition=c("Clus3_GO","classic3")
  )
)

Clus3_GO<-as.data.frame(Clus3_sResults@data)
write.table(Clus3_GO,"Clus3_GO.txt")

#############cluster 4############
Clus4_GO<-ViSEAGO::create_topGOdata(
  geneSel=selection4,
  allGenes=background,
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=5
)

classic4<-topGO::runTest(
  Clus4_GO,
  algorithm ="classic",
  statistic = "fisher"
)

Clus4_sResults<-ViSEAGO::merge_enrich_terms(
  Input=list(
    condition=c("Clus4_GO","classic4")
  )
)

Clus4_GO<-as.data.frame(Clus4_sResults@data)
write.table(Clus4_GO,"Clus4_GO.txt")

#############cluster 5############
Clus5_GO<-ViSEAGO::create_topGOdata(
  geneSel=selection5,
  allGenes=background,
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=5
)

classic5<-topGO::runTest(
  Clus5_GO,
  algorithm ="classic",
  statistic = "fisher"
)

Clus5_sResults<-ViSEAGO::merge_enrich_terms(
  Input=list(
    condition=c("Clus5_GO","classic5")
  )
)

Clus5_GO<-as.data.frame(Clus5_sResults@data)
write.table(Clus5_GO,"Clus5_GO.txt")

#######GSEA rnk file############
GSEA_COPDvsHC_ctrl<-res_COPDvsHC_ctrl[order(res_COPDvsHC_ctrl$log2FoldChange,decreasing = T),]
GSEA_COPDvsHC_ctrl<-GSEA_COPDvsHC_ctrl[,2,drop=F]
GSEA_COPDvsHC_ctrl<-GSEA_COPDvsHC_ctrl[complete.cases(GSEA_COPDvsHC_ctrl),,drop=F]
write.table(GSEA_COPDvsHC_ctrl,file="GSEA_COPDvsHC_ctrl.txt")

GSEA_COPDvsHC_Pseudo<-res_COPDvsHC_Pseudo[order(res_COPDvsHC_Pseudo$log2FoldChange,decreasing = T),]
GSEA_COPDvsHC_Pseudo<-GSEA_COPDvsHC_Pseudo[,2,drop=F]
GSEA_COPDvsHC_Pseudo<-GSEA_COPDvsHC_Pseudo[complete.cases(GSEA_COPDvsHC_Pseudo),,drop=F]
write.table(GSEA_COPDvsHC_Pseudo,file="GSEA_COPDvsHC_Pseudo.txt")

GSEA_PseudovsCtrl_HC<-res_PseudovsCtrl_HC[order(res_PseudovsCtrl_HC$log2FoldChange,decreasing = T),]
GSEA_PseudovsCtrl_HC<-GSEA_PseudovsCtrl_HC[,2,drop=F]
GSEA_PseudovsCtrl_HC<-GSEA_PseudovsCtrl_HC[complete.cases(GSEA_PseudovsCtrl_HC),,drop=F]
write.table(GSEA_PseudovsCtrl_HC,file="GSEA_PseudovsCtrl_HC.txt")

GSEA_PseudovsCtrl_COPD<-res_PseudovsCtrl_COPD[order(res_PseudovsCtrl_COPD$log2FoldChange,decreasing = T),]
GSEA_PseudovsCtrl_COPD<-GSEA_PseudovsCtrl_COPD[,2,drop=F]
GSEA_PseudovsCtrl_COPD<-GSEA_PseudovsCtrl_COPD[complete.cases(GSEA_PseudovsCtrl_COPD),,drop=F]
write.table(GSEA_PseudovsCtrl_COPD,file="GSEA_PseudovsCtrl_COPD.txt")

