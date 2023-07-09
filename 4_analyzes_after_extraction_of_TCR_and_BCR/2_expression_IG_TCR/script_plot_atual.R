library(ggplot2)
library(ggpubr)
library(cowplot)

# -- counts -- #################################################################
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

vdj_log10 <- vdj
vdj_log10$count <- log10(vdj_log10$count + 1)
vdj_log10 <- na.omit(vdj_log10)

g_count_log10 <- ggplot(vdj_log10, aes(x=receptor, y=count, fill=steroid))+
  geom_boxplot(width=c(0.4), lwd=0.2, outlier.size = 0.5)+
  #geom_jitter(color="black", size=0.4, alpha=0.9) +
  scale_y_continuous("Counts(log10)")+
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

data.frame(xtabs(~ receptor + steroid, data = vdj_log10))

# -- contagem normalizada -- ###################################################
vdj_norm <- data_vdj
vdj_norm[,1:7] <- log10(vdj_norm[,1:7]+1)
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

vdj_log10 <- rbind(IGK,IGL,IGH,TRA,TRB,TRD,TRG)

vdj_log10 <- na.omit(vdj_log10)

g_norm_log10 <- ggplot(vdj_log10, aes(x=receptor, y=count, fill=steroid))+
  geom_boxplot(width=c(0.4), lwd=0.2, outlier.size = 0.5)+
  #geom_jitter(color="black", size=0.4, alpha=0.9) +
  scale_y_continuous("Expression(log10)")+
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

data.frame(xtabs(~ receptor + steroid, data = vdj_log10))

# -- entropia -- ###############################################################
load("../afterExtraction_trust4/data.entropy.RData")

head(rownames(data.entropy))
head(metadata$sample_id)
head(gsub("_report","",rownames(data.entropy)))

idx <- match(gsub("_report","",rownames(data.entropy)),metadata$sample_id)

data.entropy$steroid <- metadata$steroid[idx]

IGK <- data.entropy[,c("IGK","steroid")]
IGL <- data.entropy[,c("IGL","steroid")]
IGH <- data.entropy[,c("IGH","steroid")]
TRA <- data.entropy[,c("TRA","steroid")]
TRB <- data.entropy[,c("TRB","steroid")]
TRG <- data.entropy[,c("TRG","steroid")]
TRD <- data.entropy[,c("TRD","steroid")]

IGK$receptor <- rep("IGK",76)
IGL$receptor <- rep("IGL",76)
IGH$receptor <- rep("IGH",76)
TRA$receptor <- rep("TRA",76)
TRB$receptor <- rep("TRB",76)
TRG$receptor <- rep("TRG",76)
TRD$receptor <- rep("TRD",76)

colnames(IGK)[1] <- "entropy"
colnames(IGL)[1] <- "entropy"
colnames(IGH)[1] <- "entropy"
colnames(TRA)[1] <- "entropy"
colnames(TRB)[1] <- "entropy"
colnames(TRD)[1] <- "entropy"
colnames(TRG)[1] <- "entropy"

vdj <- rbind(IGK,IGL,IGH,TRA,TRB,TRD,TRG)

vdj$entropy[is.na(vdj$entropy)] <- 0

vdj <- na.omit(vdj)

g_entropy <- ggplot(vdj, aes(x=receptor, y=entropy, fill=steroid))+
  geom_boxplot(width=c(0.4), lwd=0.2, outlier.size = 0.5)+
  scale_y_continuous("Shannon Entropy")+
  scale_x_discrete("")+
  scale_fill_manual(values=c("#F8766D","#00BFC4"), 
                    labels=c("HSP","LSP"))+
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

plot_grid(g_count_log10,g_norm_log10,g_entropy, nrow = 3, labels="AUTO")


data.frame(xtabs(~ receptor + steroid, data = vdj))
################################################################################
