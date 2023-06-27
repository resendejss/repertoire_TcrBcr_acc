library(ggplot2)
library(ggpubr)
library(cowplot)

load("v.counts.RData")
load("d.counts.RData")
load("j.counts.RData")

load("metadata.RData")

identical(colnames(v.counts), colnames(d.counts))
identical(colnames(v.counts), colnames(j.counts))

identical(rownames(v.counts), rownames(d.counts))
identical(rownames(v.counts), rownames(j.counts))

data_vdj <- v.counts + d.counts + j.counts

identical(rownames(data_vdj), metadata$sample_id)
idx <- match(rownames(data_vdj), metadata$sample_id)
identical(rownames(data_vdj), metadata$sample_id[idx])

data_vdj$steroid <- metadata$steroid[idx]

IGK <- data_vdj[,c("IGK","steroid")]
IGL <- data_vdj[,c("IGL","steroid")]
IGH <- data_vdj[,c("IGH","steroid")]
TRA <- data_vdj[,c("TRA","steroid")]
TRB <- data_vdj[,c("TRB","steroid")]
TRG <- data_vdj[,c("TRG","steroid")]
TRD <- data_vdj[,c("TRD","steroid")]

IGK$receptor <- rep("IGK",76)
IGL$receptor <- rep("IGL",76)
IGH$receptor <- rep("IGH",76)
TRA$receptor <- rep("TRA",76)
TRB$receptor <- rep("TRB",76)
TRG$receptor <- rep("TRG",76)
TRD$receptor <- rep("TRD",76)

colnames(IGK)[1] <- "count"
colnames(IGL)[1] <- "count"
colnames(IGH)[1] <- "count"
colnames(TRA)[1] <- "count"
colnames(TRB)[1] <- "count"
colnames(TRD)[1] <- "count"
colnames(TRG)[1] <- "count"

vdj <- rbind(IGK,IGL,IGH,TRA,TRB,TRD,TRG)

vdj_log2 <- vdj
vdj$count <- log2(vdj$count + 1)

vdj <- na.omit(vdj)

g1 <- ggplot(vdj, aes(x=receptor, y=count, fill=steroid))+
  geom_boxplot(width=c(0.4), lwd=0.2, outlier.size = 0.5)+
  scale_y_continuous("Counts(log2)")+
  scale_x_discrete("")+
  scale_fill_manual(values=c("#F8766D","#00BFC4"), 
                    labels=c("HSP (n=45)","LSP (n=30)"))+
  theme(legend.position = "bottom", 
        legend.title = element_text(face="bold", size=10),
        legend.key=element_rect(size=10, color=NA), 
        legend.key.size=unit(8,"mm"),
        legend.text=element_text(size=10), 
        legend.direction = "horizontal",
        legend.box = "horizontal" )+
  stat_compare_means(method="wilcox.test", 
                     label.x = c(1.3), 
                     label.y = 17,
                     label="p.format",
                     size=3)

vdj_norm <- data_vdj
vdj_norm[,1:7] <- log2(vdj_norm[,1:7]+1)
identical(rownames(vdj_norm),metadata$sample_id)
idx <- match(metadata$sample_id, rownames(vdj_norm))
vdj_norm[,1:7] <- vdj_norm[,1:7]/metadata$reads[idx]

IGK <- vdj_norm[,c("IGK","steroid")]
IGL <- vdj_norm[,c("IGL","steroid")]
IGH <- vdj_norm[,c("IGH","steroid")]
TRA <- vdj_norm[,c("TRA","steroid")]
TRB <- vdj_norm[,c("TRB","steroid")]
TRG <- vdj_norm[,c("TRG","steroid")]
TRD <- vdj_norm[,c("TRD","steroid")]

IGK$receptor <- rep("IGK",76)
IGL$receptor <- rep("IGL",76)
IGH$receptor <- rep("IGH",76)
TRA$receptor <- rep("TRA",76)
TRB$receptor <- rep("TRB",76)
TRG$receptor <- rep("TRG",76)
TRD$receptor <- rep("TRD",76)

colnames(IGK)[1] <- "count"
colnames(IGL)[1] <- "count"
colnames(IGH)[1] <- "count"
colnames(TRA)[1] <- "count"
colnames(TRB)[1] <- "count"
colnames(TRD)[1] <- "count"
colnames(TRG)[1] <- "count"

vdj <- rbind(IGK,IGL,IGH,TRA,TRB,TRD,TRG)

#vdj_log2 <- vdj
#vdj$count <- log2(vdj$count + 1)

vdj <- na.omit(vdj)

g2 <- ggplot(vdj, aes(x=receptor, y=count, fill=steroid))+
  geom_boxplot(width=c(0.4), lwd=0.2, outlier.size = 0.5)+
  scale_y_continuous("Expression(log2)")+
  scale_x_discrete("")+
  scale_fill_manual(values=c("#F8766D","#00BFC4"), 
                    labels=c("HSP (n=45)","LSP (n=30)"))+
  theme(legend.position = "bottom", 
        legend.title = element_text(face="bold", size=10),
        legend.key=element_rect(size=10, color=NA), 
        legend.key.size=unit(8,"mm"),
        legend.text=element_text(size=10), 
        legend.direction = "horizontal",
        legend.box = "horizontal" )+
  stat_compare_means(method="wilcox.test", 
                     label.x = c(1.3), 
                     #label.y = 17,
                     label="p.format",
                     size=3)

plot_grid(g1,g2, nrow = 2, labels="AUTO")

#data_vdj_BT <- data_vdj
#data_vdj_BT$b <- data_vdj_BT$IGK + data_vdj_BT$IGL + data_vdj_BT$IGH
#data_vdj_BT$t <- data_vdj_BT$TRB + data_vdj_BT$TRA + data_vdj_BT$TRG + data_vdj_BT$TRD

#data_vdj_BT$IGK[1] + data_vdj_BT$IGL[1] + data_vdj_BT$IGH[1]
#data_vdj_BT$b[1]
#data_vdj_BT$TRB[1] + data_vdj_BT$TRA[1] + data_vdj_BT$TRG[1] + data_vdj_BT$TRD[1]
#data_vdj_BT$t[1]

#data_vdj_BT$bnorm <- data_vdj_BT$b/metadata$reads[idx]
#data_vdj_BT$tnorm <- data_vdj_BT$t/metadata$reads[idx]








###############################################################################