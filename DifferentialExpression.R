setwd("~/Documents/piRNA_clusterAnnotation/piRNA_cluster_Annotation_Shuo/not_miRNAnormalized/Mismatch_1/pdens0.1")
#reads data tables into R
library(tidyverse)
library(tidyr)

RIL21147.R1 <- read.table("21147-R1_piRNAcluster_TotalDensity.txt",header=FALSE,col.names = c("piRNAcluster","RIL21147.R1"))
RIL21147.R2 <- read.table("21147-R2_piRNAcluster_TotalDensity.txt",header=FALSE,col.names = c("piRNAcluster","RIL21147.R2"))
RIL21147.R3 <- read.table("21147-R3_piRNAcluster_TotalDensity.txt",header=FALSE,col.names = c("piRNAcluster","RIL21147.R3"))
RIL21183.R1 <- read.table("21183-R1_piRNAcluster_TotalDensity.txt",header=FALSE,col.names = c("piRNAcluster","RIL21183.R1"))
RIL21183.R2 <- read.table("21183-R2_piRNAcluster_TotalDensity.txt",header=FALSE,col.names = c("piRNAcluster","RIL21183.R2"))
RIL21183.R3 <- read.table("21183-R3_piRNAcluster_TotalDensity.txt",header=FALSE,col.names = c("piRNAcluster","RIL21183.R3"))
RIL21188.R1 <- read.table("21188-R1_piRNAcluster_TotalDensity.txt",header=FALSE,col.names = c("piRNAcluster","RIL21188.R1"))
RIL21188.R2 <- read.table("21188-R2_piRNAcluster_TotalDensity.txt",header=FALSE,col.names = c("piRNAcluster","RIL21188.R2"))
RIL21188.R3 <- read.table("21188-R3_piRNAcluster_TotalDensity.txt",header=FALSE,col.names = c("piRNAcluster","RIL21188.R3"))
RIL21213.R1 <- read.table("21213-R1_piRNAcluster_TotalDensity.txt",header=FALSE,col.names = c("piRNAcluster","RIL21213.R1"))
RIL21213.R2 <- read.table("21213-R2_piRNAcluster_TotalDensity.txt",header=FALSE,col.names = c("piRNAcluster","RIL21213.R2"))
RIL21213.R3 <- read.table("21213-R3_piRNAcluster_TotalDensity.txt",header=FALSE,col.names = c("piRNAcluster","RIL21213.R3"))
RIL21291.R1 <- read.table("21291-R1_piRNAcluster_TotalDensity.txt",header=FALSE,col.names = c("piRNAcluster","RIL21291.R1"))
RIL21291.R2 <- read.table("21291-R2_piRNAcluster_TotalDensity.txt",header=FALSE,col.names = c("piRNAcluster","RIL21291.R2"))
RIL21291.R3 <- read.table("21291-R3_piRNAcluster_TotalDensity.txt",header=FALSE,col.names = c("piRNAcluster","RIL21291.R3"))
RIL21346.R1 <- read.table("21346-R1_piRNAcluster_TotalDensity.txt",header=FALSE,col.names = c("piRNAcluster","RIL21346.R1"))
RIL21346.R2 <- read.table("21346-R2_piRNAcluster_TotalDensity.txt",header=FALSE,col.names = c("piRNAcluster","RIL21346.R2"))
RIL21346.R3 <- read.table("21346-R3_piRNAcluster_TotalDensity.txt",header=FALSE,col.names = c("piRNAcluster","RIL21346.R3"))

#makes a big table, each column corresponds to read counts from 1 library.
all.counts <- merge(RIL21147.R1,RIL21147.R2,by="piRNAcluster",all=TRUE)
all.counts <- merge(all.counts,RIL21147.R3,by="piRNAcluster",all=TRUE)
all.counts <- merge(all.counts,RIL21183.R1,by="piRNAcluster",all=TRUE)
all.counts <- merge(all.counts,RIL21183.R2,by="piRNAcluster",all=TRUE)
all.counts <- merge(all.counts,RIL21183.R3,by="piRNAcluster",all=TRUE)
all.counts <- merge(all.counts,RIL21188.R1,by="piRNAcluster",all=TRUE)
all.counts <- merge(all.counts,RIL21188.R2,by="piRNAcluster",all=TRUE)
all.counts <- merge(all.counts,RIL21188.R3,by="piRNAcluster",all=TRUE)
all.counts <- merge(all.counts,RIL21213.R1,by="piRNAcluster",all=TRUE)
all.counts <- merge(all.counts,RIL21213.R2,by="piRNAcluster",all=TRUE)
all.counts <- merge(all.counts,RIL21213.R3,by="piRNAcluster",all=TRUE)
all.counts <- merge(all.counts,RIL21291.R1,by="piRNAcluster",all=TRUE)
all.counts <- merge(all.counts,RIL21291.R2,by="piRNAcluster",all=TRUE)
all.counts <- merge(all.counts,RIL21291.R3,by="piRNAcluster",all=TRUE)
all.counts <- merge(all.counts,RIL21346.R1,by="piRNAcluster",all=TRUE)
all.counts <- merge(all.counts,RIL21346.R2,by="piRNAcluster",all=TRUE)
all.counts <- merge(all.counts,RIL21346.R3,by="piRNAcluster",all=TRUE)
all.counts[is.na(all.counts)] <- 0
rownames(all.counts) <- all.counts$piRNAcluster
all.counts[,-1] <- all.counts[,-1] %>% round

#setwd("~/Documents/piRNA_clusterAnnotation/piRNA_cluster_Annotation_Shuo/not_miRNAnormalized/Mismatch_1/pdens0.1/42AB")
#write.csv(all.counts, "all.counts.csv")

#makes an array indicating strain
RIL <- c(rep("21147",3),rep("21183",3),rep("21188",3),rep("21213",3),rep("21291",3),rep("21346",3))
#makes an array indicating genotype
genotype <- c(rep("BB6",3),rep("BB4",3),rep("BB4",3),rep("BB6",3),rep("BB6",3),rep("BB4",3))
#makes an array indicating which pair
pair <- c(rep("A",3),rep("B",3),rep("C",3),rep("B",3),rep("C",3),rep("A",3))
#makes a table combining strain, genotype, and pair
exp.Design <- data.frame(genotype,pair,rownames=names(all.counts[-c(1)]))

#Normalizing using micro RNA counts
microRNA_norm <- read.delim("miRNA_coverage.txt", header = F, col.names = c("RIL", "norm"))
microRNA_norm$norm.Factor <- microRNA_norm$norm/1000000
write.csv(microRNA_norm, "miRNA_coverage.csv")

library(DESeq2)
#run a DEseq analysis estimating abundance of different repeats in BB6 (sensitive genotype) vs BB4 (tolerant genotype)
dds <- DESeqDataSetFromMatrix(as.matrix(all.counts[,2:19]), colData = exp.Design, design =~ pair + genotype)
sizeFactors(dds) <- microRNA_norm$norm.Factor
dds <- DESeq(dds)

#check the normalized counts. does it match manual normalization by microRNA counts.
normalized_counts <- counts(dds, normalized=TRUE) 
normalized_counts <- data.frame(normalized_counts)

#check which sizeFactor/ normalization factor they've used.
sizeFactors(dds)

#Perform differential expression between BB6/sensitive an BB4/tolerant genotypes
res <- results(dds, contrast=c("genotype","BB6","BB4")) 
res1 <- res %>% data.frame()
write.csv(normalized_counts, "normalized_counts_piRNAcluster.csv")
write.csv(all.counts, "rawcounts_piRNAcluster.csv")

#################################################################################################################

#density plot
res1$cluster <- rownames(res1) 
res1 <- res1[,c(7,1:6)]
res1$cluster <- gsub(":", ",", res1$cluster)
res1$cluster <- gsub("\\.", " ", res1$cluster, perl = T)
res1$cluster <- gsub("  ", ",", res1$cluster, perl = T)
res1 <- res1 %>%  separate(cluster, c("chr", "start", "stop"), sep = ",")

#convert all NAs to 0 so that we can calculate the mean.
res1[is.na(res1)] <- 0
data <- c(mean(res1$log2FoldChange))
data <- data.frame(data)
res2 <- res1[,c(4,5,9)]
write.csv(res2, "DiffExp_piRNAcluster.csv")

################################################################################
################################################################################

#Make PCA plot
rld <- rlog(dds)
plotPCA(rld, intgroup = c("pair" , "genotype"), ntop =497) 
#plotPCA(rld, intgroup =  "genotype", ntop = 497)
PCA_graph <- plotPCA(rld, intgroup = "genotype", ntop = 497, returnData = TRUE)
PCA_graph2 <- plotPCA(rld, intgroup = c("pair" , "genotype"), ntop =497, returnData = TRUE) 

setwd("~/Documents/Paper/Supplementals")
write.csv(PCA_graph2, "PCA_piRNAcluster.csv")

mynamestheme <- theme(plot.title = element_text(family = "Helvetica", face = "bold", size = (20), hjust = 0.5), 
                      legend.title = element_text(colour = "grey38",  face = "bold", family = "Helvetica", size = (18)), 
                      legend.text = element_text(face = "bold", colour="grey38",family = "Helvetica", size = (15)), 
                      axis.title = element_text(family = "Helvetica", size = (18), colour = "grey38",  face = "bold"),
                      axis.text = element_text(family = "Courier", colour = "slategray", size = (15),  face = "bold"))


#ggsave("PCplot_new.pdf",width=7,height=5)

ggplot(PCA_graph2, aes(PC1, PC2, group)) +  geom_point( aes(color=group, shape = group),size = 5) + 
  scale_fill_manual(name = "Strains", labels = c("A:B6:21147","A:B4:21346", "B:B4:21183","B:B6:21213", "C:B4:21188","C:B6:21291"))+
  scale_colour_manual(values=c("turquoise4", "olivedrab3","turquoise4", "olivedrab3","turquoise4", "olivedrab3"))+ 
  scale_shape_manual(values = c(19,19,3,3,17,17)) +labs(y="PC2: 22% variance", x = "PC1: 48% variance") + 
  theme_bw()+ mynamestheme + theme(legend.position="top")

ggsave("PCplot_shapes_new.pdf",width=7,height=5)

#######################################################

#In order to access rld values
rld<-assay(rlog(dds))
#write.csv(rld, "rld.csv")
#To get top genes contributing to PCA
prcomp.DESeqTransform <- function(x, ntop = 497, genes = NULL, ...){
  require(magrittr)
  # Get sample info
  sample_info <- as.data.frame(SummarizedExperiment::colData(x))
  # Get counts
  x <- SummarizedExperiment::assay(x)
  if(!is.null(genes)){
    message("Only using ", genes, " genes as requested.")
    if(!all(genes %in% rownames(x))) stop("Not all provided genes are in the gene count matrix.")
    selected_genes <- which(genes %in% rownames(x))
  } else if(is.numeric(ntop) & ntop < nrow(x)){
    ntop <- round(ntop)
    message("Only using ", ntop, " most variable genes.")
    # calculate the variance for each gene
    rv <- genefilter::rowVars(x)
    # select the ntop genes by variance
    selected_genes <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  } else {
    message("Using all ", nrow(x), " genes.")
    selected_genes <- 1:nrow(x)
  }
  # Get the data for those genes
  selected_expr <- x[selected_genes, ]
  # perform a PCA on the data in assay(x) for the selected genes
  ## Need to transpose the matrix as prcomp clusters by rows
  pca <- prcomp(t(selected_expr), ...)
  #### Eigen scores table ####
  # Get sample information from DESeq x
  # and bind the PC scores
  pc_scores <- sample_info %>%
    dplyr::bind_cols(as.data.frame(pca$x))
  #### Eigen values table ####
  eigen_values <- data.frame(PC = colnames(pca$x), stdev = pca$sdev) %>%
    dplyr::mutate(var = stdev^2,
                  var_pct = var/sum(var),
                  cum_var = cumsum(var_pct),
                  PC = forcats::fct_inorder(PC))
  #### Factor loadings table ####
  factor_loadings <- pca$rotation %>%
    as.data.frame() %>%
    dplyr::mutate(gene = row.names(.)) %>%
    dplyr::select(gene, dplyr::everything())
  #### Convert the original data to a data.frame ####
  selected_expr <- selected_expr %>%
    as.data.frame() %>%
    dplyr::rename_all(dplyr::funs(paste0("sample", .))) %>%
    dplyr::mutate(gene = rownames(.)) %>%
    dplyr::select(gene, dplyr::everything())
  # Return a list with each of these xs
  return(list(pc_scores = pc_scores,
              eigen_values = eigen_values,
              loadings = factor_loadings,
              original = selected_expr))
}

require(magrittr)

PCA1<- prcomp.DESeqTransform(rld, intgroup=c("pair" , "genotype"))
genes<- as.data.frame(PCA1$loadings)

PC_1<- subset(genes[,1:2])
PC_2<- subset(genes[,c(1,3)])

################################################################################
################################################################################

setwd("~/Documents/piRNA_clusterAnnotation/piRNA_cluster_Annotation_Shuo/not_miRNAnormalized/Mismatch_1/pdens0.1/StrainsSeparatedAnalysis/")

TopPC2clusters <- read.csv("StrainPairsAtrophyPhenotype.csv", header=T ,sep=",")
clusterInfo <- read.csv("piRNAclusterName-strandedness.csv", header=T, sep= ",")

#Top_PC2 <- merge(TopPC2clusters, clusterInfo, by="cluster")

res1$cluster <- rownames(res1)

res1 <- subset(res1, select=c("cluster", "baseMean", "log2FoldChange", "padj"))
#res1_merge <- merge(res1,clusterInfo, by="cluster")
Top_PC2 <- merge(TopPC2clusters,res1, by ="cluster")
write.csv(Top_PC2, "Top_PC2_all.csv")

################################################################################
################################################################################
#42AB : flybase: 6,255,432-6,499,291/ kelleher- 6256844-6499214
# shuos 42AB1 - chr2R:6281001..6383027
#shuos 42AB2 - chr2R:6393001..6464026
# in my data: - chr2R:6,254,004..6,496,026

cluster <- subset(res1, chr == "chr2R")
cluster <- subset(cluster, start >= 6254004)

#for all clusters

#6255432-6254004 #1.428 kb
#6499291-6496026 #3.265 kb

library(ggplot2)
mynamestheme <- theme(plot.title = element_text( size = (20), hjust = 0.5), 
                      legend.title = element_blank(), 
                      legend.text = element_text( colour="black", size = (18)), 
                      axis.title = element_text(size = (20), colour = "black"),
                      axis.text.x = element_text(colour = "black", size = (16)),
                      axis.text.y = element_text(colour = "black", size = (16))
)



ggplot(res1, aes(x=log2FoldChange)) + geom_density(alpha=.2,fill="#FF6666")  +
  theme_bw() + geom_vline(data = data,  aes(xintercept= data), linetype="dashed", size=1) + 
  geom_vline( alpha=.5,aes(xintercept= 0), linetype="solid", size=1)+ mynamestheme  + xlim(-2.5,2.5) +
  labs(x="\nLog2FoldChange (Sensitive/Tolerant)", y = "Density\n") 
  annotate("text", x = 2.5, y= 0.9, label=expression(paste(italic("p"), "= 1.804e-11")), size= 7) 
  
ggplot(res1, aes(x=log2FoldChange)) + geom_density(alpha=.2,fill="grey42")  +
  theme_bw() + geom_vline(data = data,  aes(xintercept= data), linetype="dashed", size=1) + 
  geom_vline( alpha=.5,aes(xintercept= 0), linetype="solid", size=1)+ mynamestheme  + xlim(-2.5,2.5) +
  labs(title="",x="\nLog2FoldChange (Sensitive/Tolerant)", y = "Density\n") +
  annotate("text", x = -0.5, y= 1.5, label= "***", size= 10) 


library("BSDA")
#usage is sign.test(vector of values, hypothesized median, alternative = "two.sided"/"greater"/"less",conf.level)
#our hypothesized median is zero which is the default
#this is a two sided test (default)
#conf.level = 0.95 (default)
#497
qqnorm(res1$log2FoldChange) 
shapiro.test(res1$log2FoldChange)
#W = 0.72176, p-value < 2.2e-16 Data is not normal. Left skewed 
#s = 145, p-value = 1.804e-11
#One-sample Sign-Test
SIGN.test(pull(res1,log2FoldChange))
#s = 64, p-value < 2.2e-16 ***
t.test(res1$log2FoldChange, mu = 0)


#t = 0.022648, df = 164, p-value = 0.982, mean 0.001075575



res1$start <- res1$start %>% as.numeric
#extracting the pericentromeric clusters based on the pericentromeric boundaries.
Count_2L <- subset(res1, chr == "chr2L") #55 genes
Count_2R <- subset(res1, chr == "chr2R") #91 genes
Count_2 <- rbind(Count_2L,Count_2R)

#pericentromere of 2nd chromosome in both the left and right arms.
Count_QTL2_2L <- subset(Count_2L, start >= 20820000) # 15
Count_QTL2_2R <- subset(Count_2R, start <= 6942495) #63

Count_QTL1_2L <- subset(Count_2L, start >= 19010000 & start <= 20000000) #4
Count_QTL3_2R <- subset(Count_2R, start >= 6142495  &  start <=6422495)
Count_QTL <- rbind(Count_QTL2_2L,Count_QTL2_2R,Count_QTL1_2L,Count_QTL3_2R)
#write.csv(Count_QTL2L, "piRNAcluster_QTL2L.csv")
#write.csv(Count_QTL2R, "piRNAcluster_QTL2R.csv")

Count_QTL2 <- rbind(Count_QTL2_2L, Count_QTL2_2R)
Count_QTL2$arm <- "pericentromere"

#euchromatin
Count_euchromatin_2L <- subset(Count_2L, start <= 20550000) #26 ///38 clusters
Count_euchromatin_2L <- subset(Count_euchromatin_2L, start >=5041) #37 clusters
Count_euchromatin_2R <- subset(Count_2R, start >= 5662495) #30 genes
Count_euchromatin_2R <- subset(Count_euchromatin_2R, start <= 25258060)
Count_euchromatin_2 <- rbind(Count_euchromatin_2L, Count_euchromatin_2R)

Count_euchromatin_2$arm <- "euchromatin"
Count_merge2 <- rbind(Count_QTL2, Count_euchromatin_2)

#Telomere associated region (TAS)

Count_telomere_2L <- subset(Count_2L, start <= 5041) #6
Count_telomere_2R <- subset(Count_2R, start >= 25258060 & start <= 25261551)
Count_telomere_2 <- rbind(Count_telomere_2L, Count_telomere_2R)
Count_telomere_2$arm <- "telomere"

Count_merge2 <- rbind(Count_merge2, Count_telomere_2)


data <- c(mean(Count_QTL2$log2FoldChange), mean(Count_euchromatin_2$log2FoldChange), mean(Count_telomere_2$log2FoldChange))
data <- data.frame(data)
colnames(data) <- "grp_mean"
data$arm <- c("pericentromere", "euchromatin", "telomere")
mu <- data
colnames(mu) <- c("grp_mean", "arm")


mynamestheme <- theme(plot.title = element_text( face = "bold", size = (20), hjust = 0.5), 
                      legend.title = element_text(colour = "black",  face = "bold", size = (18)), 
                      legend.text = element_text(face = "bold", colour="black", size = (15)), 
                      axis.title = element_text(size = (20), colour = "black",  face = "bold"),
                      axis.text = element_text(colour = "black", size = (16),  face = "bold"))

#peri and euchromatin arm
ggplot(Count_merge2, aes(x=log2FoldChange, fill=arm)) + geom_density(alpha=.2)  +
  theme_bw() + geom_vline(data = mu,  aes(xintercept= grp_mean,
                                                          color=arm), linetype="dashed", size=1) + 
  geom_vline( alpha=.5,aes(xintercept= 0), linetype="solid", size=1)+ mynamestheme  +
  labs(title="2nd Chromosome\n",x="\nLog2FoldChange (Sensitive/Tolerant)", y = "Density\n")  + 
  annotate("text", x = -0.5, y= 0.8, label= "***", size= 7)  +
  annotate("text", x = 0.6, y= 0.8, label= "ns", size= 7) + xlim(-2,2)

qqnorm(Count_euchromatin_2$log2FoldChange)  
#left skewed. positively skewed
shapiro.test(Count_euchromatin_2$log2FoldChange)
#W = 0.91501, p-value = 0.0002222. Data is not normal as p-value is significant.

SIGN.test(pull(Count_euchromatin_2,log2FoldChange))
#s = 12, p-value = 0.0003059 **


t.test(Count_euchromatin_2$log2FoldChange, mu = 0)
# t = -4.5436, df = 66, p-value = 2.416e-05 mu -0.4771526 


qqnorm(Count_QTL2$log2FoldChange)
#+ qqline(Count_QTL2$log2FoldChange) 
shapiro.test(Count_QTL2$log2FoldChange)
#W = 0.69305, p-value = 1.759e-11 Data is not normally distributed.
SIGN.test(pull(Count_QTL2,log2FoldChange))
#s = 8, p-value = 1.808e-12 ***


t.test(Count_QTL2$log2FoldChange, mu = 0) 
#t = -5.2671, df = 77, p-value = 1.222e-06, mu =-0.8025526 
library(car)
leveneTest(Count_merge2$log2FoldChange, Count_merge2$arm,center=mean)
var.test(Count_euchromatin_2$log2FoldChange,Count_QTL2$log2FoldChange)
#F = 0.40802, num df = 66, denom df = 77, p-value = 0.0002527
#variance varies between samples so we cannot perform two-sample t-test.
t.test(Count_euchromatin_2$log2FoldChange,Count_QTL2$log2FoldChange, var.equal = T)
#t = 1.7028, df = 143, p-value = 0.09078
library("coin")
Count_merge2_1 <- Count_merge2[-146,]
oneway_test(Count_merge2_1$log2FoldChange~as.factor(Count_merge2_1$arm))
#Z = 1.6917, p-value = 0.09071


#whole 2nd chromosome
ggplot(data=Count_2, aes(x=log2FoldChange)) +
  geom_histogram(aes(y=..density..),color="black", fill="white",binwidth=0.5)+ 
  geom_density(alpha=.2, fill="#FF6666") +
  theme_bw() + geom_vline(aes(xintercept=mean(log2FoldChange)),
                                          color="blue", linetype="dashed", size=1) +
  geom_vline( alpha=.5,aes(xintercept= 0), linetype="solid", size=1)+
  annotate("text", x = 1, y= 1, label=expression(paste(italic("p"), "= 3.178e-05")), size= 7) +
  mynamestheme + labs(title="2nd chromosome\n",x="\nLog2FoldChange", y = "Density\n")


# 3rd chromosome peri/centromere
#2l	21900000-23507622
#2R	1-6,612,495
#3L	22880000-28095748
#3R	1-4,552,934
#X	21913910-23473918

#2L :22001009 (cytogenomic)/ 22160000 (epigenomic)
#2R :5,398,184(c)/5692495 (e)
#3L	:22962476 (c)/22926900 (e)
#3R	:4,552,934 (c)/ -- (e)
#X	:22628490(c)/ 22838164 (e)

#22880000

Count_3L <- subset(res1, chr == "chr3L") #106
Count_3R <- subset(res1, chr == "chr3R") #81
Count_3 <- rbind(Count_3L,Count_3R)

Count_QTL3L <- subset(Count_3L, start >= 22926900) #28 clusters///68
Count_QTL3R <- subset(Count_3R, start <= 4552934) #12 clusters///39

Count_QTL3 <- rbind(Count_QTL3L, Count_QTL3R) #40//107
Count_QTL3$arm <- "pericentromere"

#euchromatin
Count_euchromatin_3L <- subset(Count_3L, start <= 22926900) #8///#38 clusters
Count_euchromatin_3L <- subset(Count_euchromatin_3L, start >= 19608) #37
Count_euchromatin_3R <- subset(Count_3R, start >= 4552934) #42
Count_euchromatin_3R <- subset(Count_euchromatin_3R, start <= 31173015) #40

Count_euchromatin_3 <- rbind(Count_euchromatin_3L, Count_euchromatin_3R)
Count_euchromatin_3$arm <- "euchromatin"

Count_arm3 <- rbind(Count_QTL3, Count_euchromatin_3)


#telomere associated clusters
 
Count_telomere_3L <- subset(Count_3L, start <= 19608)  #1
Count_telomere_3R <- subset(Count_3R, start >= 31173015 ) # 2 clusters in telomere #& start <= 31179331

Count_telomere_3 <- rbind(Count_telomere_3L, Count_telomere_3R)
Count_telomere_3$arm <- "telomere"
Count_arm3 <- rbind(Count_arm3, Count_telomere_3)

data3 <- c(mean(Count_QTL3$log2FoldChange), mean(Count_euchromatin_3$log2FoldChange),mean(Count_telomere_3$log2FoldChange))
data3 <- data.frame(data3)
colnames(data3) <- "grp_mean"
data3$arm <- c("pericentromere", "euchromatin", "telomere")
mu <- data3
colnames(mu) <- c("grp_mean", "arm")

#peri and euchromatin arm

ggplot(Count_arm3, aes(x=log2FoldChange, fill=arm)) + geom_density(alpha=.2)  +
  theme_bw() + geom_vline(data = mu,  aes(xintercept= grp_mean,
                                                          color=arm), linetype="dashed", size=1) + 
  mynamestheme + geom_vline( alpha=.5,aes(xintercept= 0), linetype="solid", size=1)+
  labs(title="3rd Chromosome\n",x="\nLog2FoldChange (Tolerant/Sensitive)", y = "Density\n")  + 
  annotate("text", x = -0.5, y= 1.5, label= "**", size= 7)  +
  annotate("text", x = 0.8, y= 1.5, label= "***", size= 7)
  
qqnorm(Count_QTL3$log2FoldChange) 
#+ qqline(Count_QTL3$log2FoldChange,datax=TRUE) 
shapiro.test(Count_QTL3$log2FoldChange)
#W = 0.51207, p-value < 2.2e-16
SIGN.test(pull(Count_QTL3,log2FoldChange))
#s = 78, p-value = 0.007047 **
#up s = 35, p-value = 0.01035 *

library("BSDA")
#usage is sign.test(vector of values, hypothesized median, alternative = "two.sided"/"greater"/"less",conf.level)
#our hypothesized median is zero which is the default
#this is a two sided test (default)
#conf.level = 0.95 (default)
SIGN.test(pull(Count_QTL3,log2FoldChange))
#s = 78, p-value = 0.007047 **
t.test(Count_QTL3$log2FoldChange, mu = 0) 
#t = 1.7847, df = 146, p-value = 0.07639 , 0.1057082  

qqnorm(Count_euchromatin_3$log2FoldChange)
#+ qqline(Count_QTL3$log2FoldChange,datax=TRUE)  left skewed. up.
shapiro.test(Count_euchromatin_3$log2FoldChange)
#W = 0.70924, p-value = 1.385e-07
SIGN.test(pull(Count_euchromatin_3,log2FoldChange))
#s = 27, p-value = 0.0001131 ***


t.test(Count_euchromatin_3$log2FoldChange, mu = 0)
#t = 1.9847, df = 39, p-value = 0.05424, p-value 0.2548511 

leveneTest(Count_arm3$log2FoldChange, Count_arm3$arm,center=mean)
var.test(Count_euchromatin_3$log2FoldChange,Count_QTL3$log2FoldChange)
#variance donot vary between samples so we can perform two-sample t-test.

t.test(Count_QTL3$log2FoldChange, Count_euchromatin_3$log2FoldChange, var.equal = T)
#t = -1.1318, df = 185, p-value = 0.2592



#X chromosome

Count_X <- subset(res1, chr == "chrX")  #55

Count_QTLX <- subset(Count_X, start >= 22838164) #14 ///4 genes
Count_QTLX$arm <- "pericentromere"
#euchromatin
Count_euchromatin_X <- subset(Count_X, start <= 22838164 )
Count_euchromatin_X$arm <- "euchromatin"

Count_telomere_X <- subset(res1,chr=="chrX_TAS") #1
Count_telomere_X$arm <- "telomere"
Count_armX <- rbind(Count_QTLX, Count_euchromatin_X, Count_telomere_X)
#X-TAS- LO3284

data4 <- c(mean(Count_QTLX$log2FoldChange), mean(Count_euchromatin_X$log2FoldChange), mean(Count_telomere_X$log2FoldChange))
data4 <- data.frame(data4)
colnames(data4) <- "grp_mean"
data4$arm <- c("pericentromere", "euchromatin", "telomere")
mu <- data4
colnames(mu) <- c("grp_mean", "arm")

#euchorm and pericentromere 
ggplot(Count_armX, aes(x=log2FoldChange, fill=arm)) + geom_density(alpha=.2)  +
  theme_bw() + geom_vline(data = mu,  aes(xintercept= grp_mean,
                                                          color=arm), linetype="dashed", size=1) +
  mynamestheme + geom_vline( alpha=.5,aes(xintercept= 0), linetype="solid", size=1) +
  labs(title="X Chromosome\n",x="\nLog2FoldChange (Tolerant/Sensitive)", y = "Density\n") +
  annotate("text", x = -0.3, y= 1, label= "ns", size= 5)  +
  annotate("text", x = 0.5, y= 0.96, label= "*", size= 7)

qqnorm(Count_euchromatin_X$log2FoldChange) 
#+ qqline(Count_QTL3$log2FoldChange,datax=TRUE) # right-skewed 
shapiro.test(Count_euchromatin_X$log2FoldChange)
#W = 0.75054, p-value = 7.101e-06 #data is not normally distributed
SIGN.test(pull(Count_euchromatin_X,log2FoldChange))
#s = 20, p-value = 0.01916 *
# updated s = 14, p-value = 0.02263 *

t.test(Count_euchromatin_X$log2FoldChange, mu = 0)
#t = 2.1723, df = 30, p-value = 0.03787, m = 0.2532891
#updated : t = -1.7407, df = 50, p-value = 0.08789, m = -0.1908661 

qqnorm(Count_QTLX$log2FoldChange) 
#+ qqline(Count_QTL3$log2FoldChange,datax=TRUE) # right-skewed 
shapiro.test(Count_QTLX$log2FoldChange)
#W = 0.92632, p-value = 0.08069
# The null hypothesis in shapiro test is that data is normally distributed.
#If p-value significant then we reject the null,data is not normally distributed.
#Here the p-value is not significant so we can use t-test.

#SIGN.test(pull(Count_QTLX,log2FoldChange))
#s = 12, p-value = 0.6636 *

t.test(Count_QTLX$log2FoldChange, mu = 0)
#t = 0.36704, df = 23, p-value = 0.7169, mu = 0.06712891
#updated t = 0.28695, df = 3, p-value = 0.7928 mu=0.06777926 

leveneTest(Count_armX$log2FoldChange, Count_armX$arm,center=mean)
var.test(Count_euchromatin_X$log2FoldChange,Count_QTLX$log2FoldChange)
#VAR is not different

t.test(Count_euchromatin_X$log2FoldChange,Count_QTLX$log2FoldChange, var.equal = T)
#t = 0.89369, df = 53, p-value = 0.3755



#4th chrom
#all.count <- read.csv("Diffexp_histogram.csv", sep=",", header = T)


#all.count<- subset(all.count, baseMean >=100)
Count_4 <- subset(res1, chr == "chr4") #7 genes

ggplot(data=Count_4, aes(x=log2FoldChange)) +
  
  geom_density(alpha=.2, fill="cyan") +
  theme_bw() + geom_vline(aes(xintercept=mean(log2FoldChange)),
                          color="blue", linetype="dashed", size=1) + 
  geom_vline( alpha=.5,aes(xintercept= 0), linetype="solid", size=1) + xlim(-1.5,1.5) +
  annotate("text", x = -1, y= 0.7, label=expression(paste(italic("p"), "= 0.05074")), size= 7) +
  mynamestheme + labs(title="4 chromosome\n",x="\nLog2FC (Tolerant/Sensitive)", y = "Density\n")

qqnorm(Count_4$log2FoldChange) 
#+ qqline(Count_QTL3$log2FoldChange,datax=TRUE) # left-skewed 
shapiro.test(Count_4$log2FoldChange)
#W = 0.94074, p-value = 0.6454.Data is  normal.
#SIGN.test(pull(Count_4,log2FoldChange))
#s = 32, p-value = 0.0293 *


t.test(Count_4$log2FoldChange, mu = 0) 
#t = -6.0567, df = 6, p-value = 0.000918 mean of x -1.121693 



#all pericentromere genes
Count_4$arm <- "pericentromere"
all_peri_eu_tel <- rbind(Count_arm3, Count_armX, Count_merge2, Count_4) #2375 #396
all_peri <- subset(all_peri_eu_tel, arm=="pericentromere") #150 #196
all_eu <- subset(all_peri_eu_tel, arm == "euchromatin") #2225 #195
all_tel <- subset(all_peri_eu_tel, arm == "telomere") #5
all_peri_eu_tel1 <- rbind(all_peri, all_eu,all_tel)
all_peri_eu <- rbind(all_peri, all_eu)


data <- c(mean(all_peri$log2FoldChange), mean(all_eu$log2FoldChange),mean(all_tel$log2FoldChange))
data <- data.frame(data)
colnames(data) <- "grp_mean"
data$arm <- c("pericentromere", "euchromatin", "telomere")
mu <- data
colnames(mu) <- c("grp_mean", "arm")


mynamestheme <- theme(plot.title = element_text( size = (20), hjust = 0.5), 
                      legend.title = element_blank(), 
                      legend.text = element_text( colour="black", size = (18)), 
                      axis.title = element_text(size = (20), colour = "black"),
                      axis.text.x = element_text(colour = "black", size = (16)),
                      axis.text.y = element_text(colour = "black", size = (16))
)


library(RColorBrewer)
#peri and euchromatin arm 
library(ggplot2)
A <- ggplot(all_peri_eu_tel1, aes(x=log2FoldChange)) + geom_density(alpha=.3, aes(fill=arm))  + 
  theme_bw() + geom_vline(data = mu,  aes(xintercept= grp_mean,
                                                          color=arm), linetype="dashed", size=1) + 
  geom_vline( alpha=.5,aes(xintercept= 0), linetype="solid", size=1)+ mynamestheme + xlim(-3,2.3) +
  labs(title="",x="\nLog2FoldChange (Sensitive/tolerant)", y = "Density\n") + theme(legend.position = "top") +
  annotate("text", x = -0.35, y= 0.8, label= "***", size= 10, color= "blue") +
  annotate("text", x = 0.09, y= 0.8, label= "***", size= 10, color = "red3") +
  annotate("text", x = -0.75, y= 0.87, label= "ns", size= 7, color = "chartreuse3") +
  scale_fill_manual(values =c("deepskyblue", "lightcoral", "chartreuse3")) +
  scale_color_manual(values =c("blue", "red3", "chartreuse3")) 
A


SIGN.test(pull(all_eu,log2FoldChange))
#s = 32, p-value = 3.531e-13
SIGN.test(pull(all_tel,log2FoldChange))
#s = 3, p-value = 1
SIGN.test(pull(all_peri,log2FoldChange))
#s = 12, p-value < 2.2e-16

qqnorm(all_eu$log2FoldChange) 
#+ qqline(Count_QTL3$log2FoldChange,datax=TRUE) # left-skewed 
shapiro.test(all_eu$log2FoldChange)
#W = 0.78622, p-value = 2.358e-11.Data is not normal.
SIGN.test(pull(all_eu,log2FoldChange))
#s = 66, p-value = 1.41e-06 ***
#updated s = 47, p-value = 2.9e-06 ***
t.test(all_eu$log2FoldChange, mu = 0)
#t = 4.1043, df = 109, p-value = 7.853e-05, *** 0.2648571
#t = -2.9994, df = 194, p-value = 0.00306 ** mu -0.180519

qqnorm(all_tel$log2FoldChange) 
#+ qqline(Count_QTL3$log2FoldChange,datax=TRUE) # left-skewed 
shapiro.test(all_tel$log2FoldChange)
#W = 0.92747, p-value = 0.5792. Data is normal.
SIGN.test(pull(all_tel,log2FoldChange))
t.test(all_eu$log2FoldChange, mu = 0)
#t = 4.1043, df = 109, p-value = 7.853e-05, 0.2648571
t.test(all_tel$log2FoldChange, mu = 0)
#t = -2.751, df = 4, p-value = 0.05132



qqnorm(all_peri$log2FoldChange) 
#+ qqline(Count_QTL3$log2FoldChange,datax=TRUE) # left-skewed 
shapiro.test(all_peri$log2FoldChange)
#W = 0.66009, p-value < 2.2e-16.Data is not normal.
SIGN.test(pull(all_peri,log2FoldChange))
#s = 163, p-value = 1.765e-06 ***
#updated s = 61, p-value = 1.378e-05 ***

t.test(all_peri$log2FoldChange, mu = 0)
 #t = 2.9061, df = 284, p-value = 0.003948, 0.111702
# updated t = -4.0569, df = 195, p-value = 7.19e-05, -0.1170835

library("coin")
oneway_test(all_peri_eu$log2FoldChange~as.factor(all_peri_eu$arm))

#Z = 2.4111, p-value = 0.0159

#chi-squared = 3.4645, df = 2, p-value = 0.1769
t.test(all_eu$log2FoldChange,all_peri$log2FoldChange, var.equal = T)
#t = 2.4262, df = 389, p-value = 0.01571


