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

T_LSP <- c(TRA$count[TRA$steroid=="LSP"],
           TRB$count[TRB$steroid=="LSP"],
           TRG$count[TRG$steroid=="LSP"],
           TRD$count[TRD$steroid=="LSP"])

T_HSP <- c(TRA$count[TRA$steroid=="HSP"],
           TRB$count[TRB$steroid=="HSP"],
           TRG$count[TRG$steroid=="HSP"],
           TRD$count[TRD$steroid=="HSP"])

B_LSP <- c(IGH$count[IGH$steroid=="LSP"],
           IGK$count[IGK$steroid=="LSP"],
           IGL$count[IGL$steroid=="LSP"])

B_HSP <- c(IGH$count[IGH$steroid=="HSP"],
           IGK$count[IGK$steroid=="HSP"],
           IGL$count[IGL$steroid=="HSP"])


vdj_norm <- data_vdj
#vdj_norm[,1:7] <- log2(vdj_norm[,1:7]+1)
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

T_LSP <- c(TRA$count[TRA$steroid=="LSP"],
           TRB$count[TRB$steroid=="LSP"],
           TRG$count[TRG$steroid=="LSP"],
           TRD$count[TRD$steroid=="LSP"])

T_HSP <- c(TRA$count[TRA$steroid=="HSP"],
           TRB$count[TRB$steroid=="HSP"],
           TRG$count[TRG$steroid=="HSP"],
           TRD$count[TRD$steroid=="HSP"])

B_LSP <- c(IGH$count[IGH$steroid=="LSP"],
           IGK$count[IGK$steroid=="LSP"],
           IGL$count[IGL$steroid=="LSP"])

B_HSP <- c(IGH$count[IGH$steroid=="HSP"],
           IGK$count[IGK$steroid=="HSP"],
           IGL$count[IGL$steroid=="HSP"])

T_LSP <- na.omit(T_LSP)
T_HSP <- na.omit(T_HSP)
B_LSP <- na.omit(B_LSP)
B_HSP <- na.omit(B_HSP)

wilcox.test(T_LSP, B_LSP)
mean(T_LSP)
mean(B_LSP)

wilcox.test(T_HSP, B_HSP)
mean(T_HSP)
mean(B_HSP)


wilcox.test(T_HSP, B_HSP)
mean(T_HSP)
mean(B_HSP)

wilcox.test(T_HSP, B_HSP)
mean(T_HSP)
mean(B_HSP)
