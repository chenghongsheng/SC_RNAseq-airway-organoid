library(Rsubread)

# define the reference genome fasta file
REF_GENOME <- "/gpfs1/scratch/fivan/databases/ref_data/hg38.fa"
# define the output directory for the Rsubread index
# (admin note: requires data/ref_data/download_hg19.sh to be run first)
RSUBREAD_INDEX_PATH <- "/gpfs1/scratch/fivan/databases/ref_data"
# define the basename for the index
RSUBREAD_INDEX_BASE <- "hg38"
# check what is in the reference directory
#dir(RSUBREAD_INDEX_PATH)
# build the index
buildindex(basename=file.path(RSUBREAD_INDEX_PATH,RSUBREAD_INDEX_BASE), reference=REF_GENOME)

# Perform fastqc check
## Set input and output directory
indir <- '/media/DP/COPDbulk/data/FASTQ'
outdir <- '/media/DP/COPDbulk/data/FastQC'
if (!dir.exists(outdir)) try(dir.create(outdir, recursive= TRUE), silent= TRUE)

# Sample IDs
files <- list.files(indir)
files <- files[grep('\\.fastq\\.gz', files)]

## Perform fastqc
for (f in files) {
  
  cmd <- paste0('fastqc --outdir=', outdir, ' ', file.path(indir, f))
  system(cmd)
  
}


library(Rsubread)

args <- commandArgs(trailingOnly=TRUE)
lwr_idx <- args[1]
upp_idx <- args[2]

#REF_GENOME <- "/home/fivan/Databases/ref_data/hg38.fa"
#RSUBREAD_INDEX_PATH <- "/home/fivan/Databases/ref_data"
#RSUBREAD_INDEX_BASE <- "hg38"

#wkdir <- "/mediaDP/COPDbulk/data"

REF_GENOME <- "/gpfs1/scratch/fivan/databases/ref_data/hg38.fa"
RSUBREAD_INDEX_PATH <- "/gpfs1/scratch/fivan/databases/ref_data"
RSUBREAD_INDEX_BASE <- "hg38"

wkdir <- "/gpfs1/scratch/fivan/COPDbulk/data"

# define shared directory for RNAseq data
RNAseqDATADIR <- file.path(wkdir, "FASTQ")
outdir <- file.path(wkdir, "BAM")
if (!dir.exists(outdir)) try(dir.create(outdir, recursive= TRUE), silent= TRUE)

RNAseqIDs <- unique(gsub("\\.fastq\\.gz", "", list.files(RNAseqDATADIR)))

for (i in lwr_idx:upp_idx) {
  
  ID <- RNAseqIDs[i]
  
  # define the fastq file with forward reads
  inputfile <- file.path(RNAseqDATADIR,paste0(ID, ".fastq.gz"))
  
  # run the align command to map the reads
  align(index=file.path(RSUBREAD_INDEX_PATH,RSUBREAD_INDEX_BASE), readfile1=inputfile,
        output_file=file.path(outdir, paste0(ID, ".bam")))
  
}


library(Rsubread)

args <- commandArgs(trailingOnly=TRUE)
nfiles <- args[1]

REF_GENOME <- "/gpfs1/scratch/fivan/databases/ref_data/hg38.fa"
RSUBREAD_INDEX_PATH <- "/gpfs1/scratch/fivan/databases/ref_data"
RSUBREAD_INDEX_BASE <- "hg38"

wkdir <- "/gpfs1/scratch/fivan/COPDbulk/data"

# define shared directory for RNAseq data
BAMDIR <- file.path(wkdir, "BAM")
outdir <- wkdir

BAMIDs <- list.files(BAMDIR)
BAMfiles <- BAMIDs[grep("\\.bam$", BAMIDs)]

outputBAMfile <- file.path(BAMDIR, BAMfiles[1])
df <- propmapped(outputBAMfile)

for (i in 2:nfiles) {
  
  outputBAMfile <- file.path(BAMDIR, BAMfiles[i])
  df <- rbind(df, propmapped(outputBAMfile))
  
}

write.csv(df, file.path(outdir, "propmapped.csv"))


library(Rsubread)

args <- commandArgs(trailingOnly=TRUE)
lwridx <- args[1]
uppidx <- args[2]

REF_GENOME <- "/gpfs1/scratch/fivan/databases/ref_data/hg38.fa"
RSUBREAD_INDEX_PATH <- "/gpfs1/scratch/fivan/databases/ref_data"
RSUBREAD_INDEX_BASE <- "hg38"

wkdir <- "/gpfs1/scratch/fivan/COPDbulk/data"

# define shared directory for RNAseq data

BAMDIR <- file.path(wkdir, "BAM")
outdir <- file.path(wkdir, "featureCounts")
if (!dir.exists(outdir)) try(dir.create(outdir, recursive= TRUE), silent= TRUE)

BAMIDs <- list.files(BAMDIR)
BAMfiles <- BAMIDs[grep("\\.bam$", BAMIDs)]

for (i in lwridx:uppidx) {
  
  F <- BAMfiles[i]
  outputBAMfile <- file.path(BAMDIR, F)
  mycounts <- featureCounts(outputBAMfile, annot.ext=file.path(RSUBREAD_INDEX_PATH,"hg38.ensGene.gtf"),
                            isGTFAnnotationFile=TRUE, isPairedEnd=FALSE)
  saveRDS(mycounts, file.path(outdir, paste0(gsub("\\.bam", "", F), ".RData")))
  
}


# get table for feature counts

#library(Rsubread)
#wkdir <- "/gpfs1/scratch/fivan/COPDbulk/data"
#indir <- file.path(wkdir, "featureCounts")
#
#outdir <- file.path(wkdir, "featureCountsTables")
#if (!dir.exists(outdir)) try(dir.create(outdir, recursive= TRUE), silent= TRUE)
#
#FILES <- list.files(indir)
#
#for (F in FILES) {
#
#  mycounts <- readRDS(file.path(indir, F))
#  write.table(mycounts$counts, file= file.path(outdir, paste0(gsub("_.*", "", F), ".txt")), sep="\t", quote=FALSE, append=FALSE)
#
#}



# get table for feature counts

library(Rsubread)
wkdir <- "/gpfs1/scratch/fivan/COPDbulk/data"
indir <- file.path(wkdir, "featureCounts")

outdir <- file.path(wkdir, "featureCountsTables")
if (!dir.exists(outdir)) try(dir.create(outdir, recursive= TRUE), silent= TRUE)

FILES <- list.files(indir)

for (F in FILES) {
  
  mycounts <- readRDS(file.path(indir, F))
  write.table(mycounts$counts, file= file.path(outdir, paste0(gsub("_.*", "", F), ".txt")), sep="\t", quote=FALSE, append=FALSE)
  
}

# merge tables 
library(Rsubread)
wkdir <- "/gpfs1/scratch/fivan/COPDbulk/data"
indir <- file.path(wkdir, "featureCountsTables")

outdir <- wkdir

FILES <- list.files(indir)
counts <- read.table(file.path(indir, FILES[1]))
for (s in 2:length(FILES)){
  counts <- cbind(counts, read.table(file.path(indir, FILES[s])))
}

write.csv(counts, file= file.path(outdir, "featureCountsTables.csv"))


# Experimental data

wkdir <- "/gpfs1/scratch/fivan/COPDbulk/data"
wkdir <- "Y:/FXI-BACKUPS-SHOTGUN/COPDbulk/data"
outdir <- wkdir

# Experimental data file

df0 <- read.csv(file.path(wkdir, "bulkRNAseq_metadata.csv"), header= TRUE)
df1 <- read.csv(file.path(wkdir, "GSMbulkRNAseq.csv"), header= TRUE)
df2 <- read.csv(file.path(wkdir, "SraRunTableGSE146532.txt"), header= TRUE)
df3 <- read.csv(file.path(wkdir, "SraRunTableGSE162154.txt"), header= TRUE)
df4 <- read.csv(file.path(wkdir, "GSE124180.txt"), header= TRUE)

# Combine df2 and df3 to df
df3$time_point <- "0 hours"
df2$SAMPLE_TYPE <- paste0(df2$treatment, " ", df2$diagnosis, " ", df2$SUBJECT_ID, " ", df2$time_point)
cols <- intersect(colnames(df3), colnames(df2))
df <- rbind(df2[, cols], df3[, cols])
df$Subject <- gsub(" .*", "", gsub(".*ubject_|.*ubject ", "", df$SAMPLE_TYPE))
df$Subject[121:126] <- c("01", "02", "03", "04", "05", "06")

# Combine df and df1 to df
df <- df[df$Sample.Name %in% df1$GSM,]
rownames(df) <- df$Sample.Name
rownames(df1) <- df1$GSM
df <- cbind(df, df1[rownames(df), ])
df <- df[, c("Run","GSE","GSM","Instrument","SAMPLE_TYPE","Metadata","source_name","time_point","Subject")]
rownames(df) <- df$Run

df$Experiment0 <- gsub(" ", ".", gsub(" hours", "h", gsub("Uninfected Healthy Subject_", "Healthy",
                                                          gsub("Uninfected COPD Subject_", "COPD", df$SAMPLE_TYPE))))
df$Experiment <- gsub("Healthy.*|never.smoked", "Healthy", gsub("COPD.*", "COPD", df$Experiment0))
df <- df[, c("Run","GSM","GSE","Instrument","source_name","time_point","Experiment","Subject")]

# Combine df and df4 to df
df4$Subject <- c(paste0("0", 1:9), 10:21)
df4$Instrument <- "Illumina HiSeq 2500"
df4$source_name <- df4$Sample_source
df4$time_point <- "0 hours"
df4$Experiment <- df4$Group
df4$Run <- df4$GSM
df4 <- df4[, c("Run","GSM","GSE","Instrument","source_name","time_point", "Experiment","Subject")]
df <- rbind(df, df4)

# Combine df and df0 to df
rownames(df0) <- df0$SampleID
df0$Subject <- c(paste0("0", 1:9), 10:12)
df0$GSE <- "SHRlab"
df0$Instrument <- "Illumina HiSeq 2500"
df0$source_name <- "Nasal pharyngal"
df0$time_point <- "6 hours"
df0$Experiment <- gsub("NPO-HC", "Healthy", gsub("NPO-COPD", "COPD", df0$Sample_Details))
df0$GSM <- df0$SampleID
df0$Run <- df0$SampleID
df0 <- df0[, c("Run","GSM","GSE","Instrument","source_name","time_point","Experiment","Subject")]
df0 <- df0[grep("NRMK", rownames(df0)), ]
df <- rbind(df, df0)

df$Subject <- paste0("S",df$Subject)

# Save df as experiment_design.csv
write.csv(df, file= file.path(outdir, "experiment_design.csv"))

#experiment_design <- df
#samples <- as.character(experiment_design$SampleID)
#group <- factor(experiment_design$Experiment)
#group


# QC and stats

library(Rsubread)
wkdir <- "/gpfs1/scratch/fivan/COPDbulk/data"
wkdir <- "Y:/FXI-BACKUPS-SHOTGUN/COPDbulk/data"
indir <- file.path(wkdir, "featureCountsTables")
outdir <- wkdir


# read df as experiment_design.csv

experiment_design <- read.csv(file.path(wkdir, "experiment_design.csv"), row.names= 1)

experiment_design$source_name <- gsub("Nasal.*", "NP", gsub("large.*", "LA", gsub("Primary.*", "SA", gsub("Bronch.*", "BE", experiment_design$source_name))))
experiment_design$time_point <- gsub(" hours", "", experiment_design$time_point)
experiment_design$Experiment <- gsub("COPD", "CO", gsub("Healthy", "HE", experiment_design$Experiment))
experiment_design$group <- paste0(experiment_design$source_name, ".", experiment_design$Experiment, gsub("S", "", experiment_design$Subject), ".", experiment_design$time_point)
experiment_design <- experiment_design[sort(as.vector(experiment_design$group), index.return= TRUE)$ix, ]
samples <- experiment_design$group
samples
group <- factor(experiment_design$group)
group


# FEATURE COUNTS

# Rsubread counts
df <- read.csv(file.path(wkdir, "featureCountsTables.csv"))
colnames(df) <- gsub(".bam", "", colnames(df))
# GSE124180
df1 <- read.csv(file.path(wkdir, "GSMbulkRNAseq.csv"), header= TRUE)
df4 <- read.csv(file.path(wkdir, "GSE124180.txt"), header= TRUE)
df4 <- df4[df4$Metadata %in% df1$Metadata, ]
GSE124180 <- read.csv(file.path(wkdir, "NCBIGeneCount", "GSE124180_gene_count_table.tsv"), header= TRUE, sep= "\t", row.names= 1)
GSE124180 <- GSE124180[-1,]
myrows <- rownames(GSE124180)
GSE124180 <- as.data.frame(apply(GSE124180, 2, as.numeric))
rownames(GSE124180) <- myrows
GSE124180 <- GSE124180[, as.character(df4$Metadata)]
GSE124180 <- data.frame(X=myrows, GSE124180)
head(GSE124180)
colnames(GSE124180) <- c("X", rownames(df4))
head(GSE124180)
# Merge
df <- merge(x= df, y= GSE124180, by.x= "X", by.y= "X", all= TRUE)
df[is.na(df)] <- 0 # -------------------------------------------------------- MAKE IT 0 IF NA!!!


rownames(df) <- df$X
df <- df[,-1]
df <- df[, rownames(experiment_design)]
colnames(df) <- as.vector(group)

my.colors <- as.numeric(experiment_design$GSE)
my.colors <- sapply( as.vector(experiment_design$GSE), function(x) { if (x=="GSE146532") "black" else if (x=="GSE124180") "blue" else if (x=="SHRlab") "red" else "green"} )

rownames(experiment_design) <- experiment_design$group



### UNNORMALIZED
### --------------------------------------------------------------------------

# Density plot of raw read counts (log10)

png(file=file.path(outdir, "Raw_read_counts_per_gene.density.png"))
counts <- df[,1]
logcounts <- log(counts,10) 
d <- density(logcounts)
plot(d,xlim=c(1,8),main="",ylim=c(0,.2),xlab="Raw read counts per gene (log10)", ylab="Density", col= my.colors[1])
for (s in 2:ncol(df)){
  counts <- df[,s]
  logcounts <- log(counts,10) 
  d <- density(logcounts)
  lines(d, col= my.colors[s])
}
legend(5, 0.15, legend=c("GSE146532", "GSE124180","SHRlab","GSE162154"),
       col=c("black", "blue", "red", "green"), lty=1, cex=0.8)
dev.off()

# Box plots of raw read counts (log10)

counts <- df
png(file=file.path(outdir, "Raw_read_counts_per_gene.boxplot.png"), width= 2000, height= 800)
logcounts <- log(counts+1,10)
boxplot(logcounts, main="", xlab="", ylab="Raw read counts per gene (log10)",axes=FALSE)
axis(2)
axis(1,at=c(1:ncol(df)),labels=colnames(logcounts),las=2,cex.axis=0.8)
dev.off()


# Heatmap

## select data for the 100 most highly expressed genes
select = order(rowMeans(counts), decreasing=TRUE)[1:1000]
highexprgenes_counts <- as.matrix(counts[select,])

## heatmap with sample name on X-axis
png(file=file.path(outdir, "High_expr_genes.heatmap.png"), width= 2000, height= 800)
heatmap(highexprgenes_counts, col=topo.colors(50), margin=c(10,6), Colv= NA)
dev.off()

## heatmap with condition group as labels
highexprgenes_counts <- highexprgenes_counts[, samples]
colnames(highexprgenes_counts)<- group
png(file=file.path(outdir, "High_expr_genes.heatmap.samples.png"), width= 2000, height= 800)
heatmap(highexprgenes_counts, col = topo.colors(50), margin=c(10,6))
dev.off()


# PCA

group <- factor(experiment_design$Experiment)
group

# select data for the 1000 most highly expressed genes
select = order(rowMeans(counts), decreasing=TRUE)[1:100]
highexprgenes_counts <- counts[select,]
# annotate the data with condition group as labels
colnames(highexprgenes_counts)<- group
# transpose the data to have variables (genes) as columns
data_for_PCA <- t(highexprgenes_counts)
dim(data_for_PCA)

## calculate MDS (matrix of dissimilarities)
mds <- cmdscale(dist(data_for_PCA), k=3, eig=TRUE)  
# k = the maximum dimension of the space which the data are to be represented in
# eig = indicates whether eigenvalues should be returned

mds$eig

# transform the Eigen values into percentage
eig_pc <- mds$eig * 100 / sum(mds$eig)
# plot the PCA
png(file=file.path(outdir, "PCA_PropExplainedVariance.png"))
barplot(eig_pc,
        las=1,
        xlab="Dimensions", 
        ylab="Proportion of explained variance (%)", y.axis=NULL,
        col="darkgrey")
dev.off()

## calculate MDS
mds <- cmdscale(dist(data_for_PCA)) # Performs MDS analysis 
#Samples representation
png(file=file.path(outdir, "PCA_Dim1vsDim2.png"))
plot(mds[,1], -mds[,2], type="n", xlab="Dimension 1", ylab="Dimension 2", main="")
text(mds[,1], -mds[,2], rownames(mds), cex=0.8) 
dev.off()


### NORMALIZED
### --------------------------------------------------------------------------

all( colnames(df) %in% rownames(experiment_design) )
all( colnames(df) == rownames(experiment_design) )

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = df, colData = experiment_design, design = ~ source_name+Experiment)
View(counts(dds))
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
df <- normalized_counts

# Density plot of raw read counts (log10)

png(file=file.path(outdir, "Raw_read_counts_per_gene.density.normalized.png"))
counts <- df[,1]
logcounts <- log(counts,10) 
d <- density(logcounts)
plot(d,xlim=c(1,8),main="",ylim=c(0,.2),xlab="Raw read counts per gene (log10)", ylab="Density", col= my.colors[1])
for (s in 2:ncol(df)){
  counts <- df[,s]
  logcounts <- log(counts,10) 
  d <- density(logcounts)
  lines(d, col= my.colors[s])
}
legend(5, 0.15, legend=c("GSE146532", "GSE124180","SHRlab","GSE162154"),
       col=c("black", "blue", "red", "green"), lty=1, cex=0.8)
dev.off()

# Box plots of raw read counts (log10)

counts <- df
png(file=file.path(outdir, "Raw_read_counts_per_gene.boxplot.normalized.png"), width= 2000, height= 800)
logcounts <- log(counts+1,10)
boxplot(logcounts, main="", xlab="", ylab="Raw read counts per gene (log10)",axes=FALSE)
axis(2)
axis(1,at=c(1:ncol(df)),labels=colnames(logcounts),las=2,cex.axis=0.8)
dev.off()

#DEG analysis

dds<-DESeq(dds)
vsd<-varianceStabilizingTransformation(dds, blind=TRUE)

res<-results(dds, contrast =c('Experiment',"COPD","Healthy"),alpha = 0.05 )
res<-res[order(res$padj),]
summary(res)

sig.res<-res[which((res$padj<0.05)& (abs(res$log2FoldChange)>1)),]




