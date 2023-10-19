#Generating sliding window plots that depicts expression of two piRNA clusters named cluster1 and cluster2 in sensitive and tolerant genoypes

setwd("~/Documents/piRNA_clusterAnnotation/piRNA_cluster_Annotation_Shuo/not_miRNAnormalized/Mismatch_1/pdens0.1/slidingWindow")

library(dplyr)
library(ggplot2)

#For cluster 1 
RIL21183.p <- read.delim("21183.QTLcluster.p.bed", header = FALSE, col.names = c("chr", "pos", "coverage"))
RIL21213.p <-read.delim("21213.QTLcluster.p.bed", header = FALSE, col.names = c("chr", "pos", "coverage"))

RIL21183.n <- read.delim("21183.QTLcluster.n.bed", header = FALSE, col.names = c("chr", "pos", "coverage"))
RIL21213.n <-read.delim("21213.QTLcluster.n.bed", header = FALSE, col.names = c("chr", "pos", "coverage"))

#For cluster 2
RIL21183_cluster2.p <- read.delim("21183.QTLcluster2.p.bed", header = FALSE, col.names = c("chr", "pos", "coverage"))
RIL21213_cluster2.p <-read.delim("21213.QTLcluster2.p.bed", header = FALSE, col.names = c("chr", "pos", "coverage"))

RIL21183_cluster2.n <- read.delim("21183.QTLcluster2.m.bed", header = FALSE, col.names = c("chr", "pos", "coverage"))
RIL21213_cluster2.n <-read.delim("21213.QTLcluster2.m.bed", header = FALSE, col.names = c("chr", "pos", "coverage"))

###################################################################################################
#For cluster 1
RIL21183.p$strand <- "+"
RIL21183.n$strand <- "-"

RIL21183.p$coverage <- log(RIL21183.p$coverage + 1, base = 2) 
RIL21183.n$coverage <- -log(RIL21183.n$coverage + 1, base = 2) 
RIL21183 <- rbind(RIL21183.p,RIL21183.n)
RIL21183$genotype <- "Tolerant"

RIL21213.p$strand <- "+"
RIL21213.n$strand <- "-"
RIL21213.p$coverage <- log(RIL21213.p$coverage + 1, base = 2) 
RIL21213.n$coverage <- -log(RIL21213.n$coverage + 1, base = 2) 
RIL21213 <- rbind(RIL21213.p,RIL21213.n)
RIL21213$genotype <- "Sensitive"

df <- rbind(RIL21183,RIL21213)

###################################################################################################
#For cluster 2

RIL21183_cluster2.p$strand <- "+"
RIL21183_cluster2.n$strand <- "-"
RIL21183_cluster2.p$coverage <- log(RIL21183_cluster2.p$coverage + 1, base = 2) 
RIL21183_cluster2.n$coverage <- -log(RIL21183_cluster2.n$coverage + 1, base = 2) 
RIL21183_cluster2 <- rbind(RIL21183_cluster2.p,RIL21183_cluster2.n)
RIL21183_cluster2$genotype <- "Tolerant"

RIL21213_cluster2.p$strand <- "+"
RIL21213_cluster2.n$strand <- "-"
RIL21213_cluster2.p$coverage <- log(RIL21213_cluster2.p$coverage + 1, base = 2) 
RIL21213_cluster2.n$coverage <- -log(RIL21213_cluster2.n$coverage + 1, base = 2) 
RIL21213_cluster2 <- rbind(RIL21213_cluster2.p,RIL21213_cluster2.n)
RIL21213_cluster2$genotype <- "Sensitive"

df2 <- rbind(RIL21183_cluster2,RIL21213_cluster2)
###################################################################################################
#uploading csv file that has TE annotations within the piRNA cluster of interest
TE_cluster1 <- read.csv("TE-QTLcuster1.csv", sep = ",", header = TRUE)
TE_cluster1$genotype <- factor("extra",levels = c("extra","Sensitive","Tolerant"))

extra1<- subset(df, genotype=="Tolerant")
extra1$genotype <- "extra"
extra1$coverage <- 0
df <- rbind(df, extra1)

library(RColorBrewer)
# Define the number of colors you want
nb.cols <- 13
mycolors <- colorRampPalette(brewer.pal(9, "Set1"))(nb.cols)
mycolors <- c("#E41A1C","#BA5E6C", "#999999","#4DAF4A","#7E6E85","#0066FF","#3E8E93","#FF7F00","#FFD421","#E1C62F","#A65628","#DB728C","#D789B2")

#Plotting
df$genotype <- factor(df$genotype, levels = c("extra","Sensitive","Tolerant"))
a<-  ggplot(df, aes(x=pos,y=coverage)) + geom_bar(stat='identity', fill ="black") + 
  ylab("log2(RPM+1)\n") + xlab("\nNucleotide Position") + theme_classic()  +
  theme(axis.text=element_text(size=13, color = "black"),axis.title=element_text(size=24),strip.text = element_text(size=20))  +
  geom_rect(data=TE_cluster1, mapping = aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=Name), color= "white", inherit.aes=FALSE) +
  scale_fill_manual(values = mycolors) + ylim(-13,13) +
  facet_grid(genotype~., scales = "free_y")
a 

###################################################################################################

library(ggplot2)
TE_cluster2 <- read.csv("TE-QTLcluster2.csv", sep = ",", header = TRUE)
TE_cluster2$genotype <- factor("extra",levels = c("extra","Sensitive","Tolerant"))
#TE_cluster2$col <- ifelse(TE_cluster2$col, 'white', 'gray53')
#Plotting
extra<- subset(df2, genotype=="Tolerant")
extra$genotype <- "extra"
extra$coverage <- 0
df2 <- rbind(df2, extra)
df2$genotype <- factor(df2$genotype, levels = c("extra","Sensitive","Tolerant"))
b <-  ggplot(df2, aes(x=pos,y=coverage)) + geom_bar(stat='identity', fill ="black") + 
  ylab("log2(RPM+1)\n") + xlab("\nNucleotide Position") + theme_classic()  +
  theme(axis.text=element_text(size=15, color = "black"),axis.title=element_text(size=24),strip.text = element_text(size=20))  +
  geom_rect(data=TE_cluster2, mapping = aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=Name), color= "white", inherit.aes=FALSE) +
  scale_fill_manual(values = mycolors) + ylim(-13,13) +
  facet_grid(genotype~., scales = "free_y")
b
