setwd("~/Documents/Shuo_smallRNA_heatmap")

piRNAcluster <- read.csv("micor_normalized_clusterExp.csv", header = T, sep=",")
piRNAcluster <- read.csv("cluster_heatMap.csv", header = T, sep=",")

piRNAcluster <- piRNAcluster[,-c(1:4,11)]

#plot sample
piRNAcluster<- piRNAcluster[,c(1,6,4,2,5,3)]
colnames(piRNAcluster) <- c("Sensitive 1", "Tolerant 1", "Sensitive 2", "Tolerant 2", "Sensitive 3", "Tolerant 3")
library(RColorBrewer)
library(gplots)

piRNA.exp_matrix <- data.matrix(piRNAcluster)
normalized <- log(piRNA.exp_matrix + 1, base = 2) 

#For selected piRNA clusters
colfunc <- colorRampPalette(c("black", "white", "red"))
dev.off()

heatmap.2(normalized,trace="none",dendrogram = 'none',Rowv = NA, Colv = NA, #dendogram removal
          keysize =1.2,key.title = "A", key.xlab = expression("log"[2]*" (RPKM + 1)"), key.ylab = NA,key.ytickfun=NA,
          density.info = c('none'), #no frequency 
          #key.par=list(mar=c(4,1,2.5,1), cex=1, cex.lab=1.2, cex.axis=1.2), #determine the key size
          key.par=list(mgp=c(1.8, 0.7, 0),
                       mar=c(0.2, 2, 2, 2), cex.lab=1.5, cex.axis=1.5),
          col = colorRampPalette(c('black', 'red'))(25), mar = c(8,20),
          cexRow=0.5, cexCol =1.7, srtCol=45, #cexRow and cexCol regulate font size; srtCol determines the slop of xlabel
          adjCol = c(0.9, 0.5), lhei = c(0.88,6), adjRow=c(0,1), hclustfun = function(d) hclust(d, method = "average")
)

#ggsave("heatmap_selected.pdf", width=8.5, height =7)
