
library(ggplot2)
library(ggpubr)
library(cowplot)

getwd()
dir()

load("data.entropy.RData")
load("metadata.RData")

head(data.entropy)
gsub("_report","", rownames(data.entropy)) %in%  metadata$sample_id

idx <- match(gsub("_report","", rownames(data.entropy)), metadata$sample_id)

data <- data.entropy
data$steroid <- metadata$steroid[idx]

IGK <- data[,c("IGK","steroid")]
IGL <- data[,c("IGL","steroid")]
IGH <- data[,c("IGH","steroid")]
TRA <- data[,c("TRA","steroid")]
TRB <- data[,c("TRB","steroid")]
TRD <- data[,c("TRD","steroid")]
TRG <- data[,c("TRG","steroid")]

IGK$receptor <- rep("IGK", 76)
IGL$receptor <- rep("IGL", 76)
IGH$receptor <- rep("IGH", 76)
TRA$receptor <- rep("TRA", 76)
TRB$receptor <- rep("TRB", 76)
TRD$receptor <- rep("TRD", 76)
TRG$receptor <- rep("TRG", 76)

colnames(IGK)[1] <- "entropy"
colnames(IGL)[1] <- "entropy"
colnames(IGH)[1] <- "entropy"
colnames(TRA)[1] <- "entropy"
colnames(TRB)[1] <- "entropy"
colnames(TRD)[1] <- "entropy"
colnames(TRG)[1] <- "entropy"

entropy <- rbind(IGK,IGL,IGH,TRA,TRB,TRD,TRG)

entropy <- na.omit(entropy)

g1 <- ggplot(entropy, aes(x=receptor, y=entropy, fill=steroid))+
  ggtitle("Shannon Entropy")+
  geom_boxplot(width=0.4, lwd=0.2, outlier.size = 0.5)+
  scale_y_continuous("Shannon Entropy")+
  scale_x_discrete("")+
  scale_fill_manual(values = c("#F8766D","#00BFC4"))+
  theme(legend.position = "bottom", 
        legend.title = element_text(face="bold", size=10),
        legend.key=element_rect(size=10, color=NA), 
        legend.key.size=unit(8,"mm"),
        legend.text=element_text(size=10), 
        legend.direction = "horizontal",
        legend.box = "horizontal" )+
  stat_compare_means(method = "wilcox.test",
                     #label.x = c(1.3),
                     label = "p.format",
                     size = 3)

#############
load("../2_expression_IG_TCR/v.counts.RData")
load("../2_expression_IG_TCR/d.counts.RData")
load("../2_expression_IG_TCR/j.counts.RData")

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

g2 <- ggplot(vdj, aes(x=receptor, y=count, fill=steroid))+
  ggtitle("TCR and BCR count")+
  geom_boxplot(width=c(0.4), lwd=0.2, outlier.size = 0.5)+
  scale_y_continuous("Counts(log2)")+
  scale_x_discrete("")+
  scale_fill_manual(values=c("#F8766D","#00BFC4"), 
                    labels=c("HSP (n=45)","LSP (n=30)"))+
  theme(legend.position = "none", 
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

g3 <- ggplot(vdj, aes(x=receptor, y=count, fill=steroid))+
  ggtitle("Expression of TCR and BCR")+
  geom_boxplot(width=c(0.4), lwd=0.2, outlier.size = 0.5)+
  scale_y_continuous("Expression(log2)")+
  scale_x_discrete("")+
  scale_fill_manual(values=c("#F8766D","#00BFC4"), 
                    labels=c("HSP (n=45)","LSP (n=30)"))+
  theme(legend.position = "none", 
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

plot_grid(g2,g3,g1, nrow = 3, labels="AUTO")

################################################################################

