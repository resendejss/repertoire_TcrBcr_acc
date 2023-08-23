################################################################################
# Name: entropy_TCR_BCR_boxplot                                             #
#                                                                              #
# Author: Jean Resende (jean.s.s.resende@gmail.com)                            #
# Creation date: 2023/08/23                                                    #
# Last update date: 2023/08/23                                                 #
################################################################################
library(ggplot2)
library(ggpubr)
library(cowplot)


load("data.entropy.RData")
load("metadata.RData")

head(data.entropy)
head(metadata)

data <- data.entropy[,c(13,14,16,15,1,11,12)]

gsub("_report","", rownames(data)) %in%  metadata$sample_id
idx <- match(gsub("_report","", rownames(data)), metadata$sample_id)
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

# head(data.entropy)
# gsub("_report","", rownames(data.entropy)) %in%  metadata$sample_id

# idx <- match(gsub("_report","", rownames(data.entropy)), metadata$sample_id)

# data <- data.entropy
# data$steroid <- metadata$steroid[idx]

# IGK <- data[,c("IGK","steroid")]
# IGL <- data[,c("IGL","steroid")]
# IGH <- data[,c("IGH","steroid")]
# TRA <- data[,c("TRA","steroid")]
# TRB <- data[,c("TRB","steroid")]
# TRD <- data[,c("TRD","steroid")]
# TRG <- data[,c("TRG","steroid")]

# IGK$receptor <- rep("IGK", 76)
# IGL$receptor <- rep("IGL", 76)
# IGH$receptor <- rep("IGH", 76)
# TRA$receptor <- rep("TRA", 76)
# TRB$receptor <- rep("TRB", 76)
# TRD$receptor <- rep("TRD", 76)
# TRG$receptor <- rep("TRG", 76)

# colnames(IGK)[1] <- "entropy"
# colnames(IGL)[1] <- "entropy"
# colnames(IGH)[1] <- "entropy"
# colnames(TRA)[1] <- "entropy"
# colnames(TRB)[1] <- "entropy"
# colnames(TRD)[1] <- "entropy"
# colnames(TRG)[1] <- "entropy"

# entropy <- rbind(IGK,IGL,IGH,TRA,TRB,TRD,TRG)

# entropy <- na.omit(entropy)
# 
# g1 <- ggplot(entropy, aes(x=receptor, y=entropy, fill=steroid))+
#   ggtitle("Shannon Entropy")+
#   geom_boxplot(width=0.4, lwd=0.2, outlier.size = 0.5)+
#   scale_y_continuous("Shannon Entropy")+
#   scale_x_discrete("")+
#   scale_fill_manual(values = c("#F8766D","#00BFC4"))+
#   theme(legend.position = "bottom", 
#         legend.title = element_text(face="bold", size=10),
#         legend.key=element_rect(size=10, color=NA), 
#         legend.key.size=unit(8,"mm"),
#         legend.text=element_text(size=10), 
#         legend.direction = "horizontal",
#         legend.box = "horizontal" )+
#   stat_compare_means(method = "wilcox.test",
#                      #label.x = c(1.3),
#                      label = "p.format",
#                      size = 3)


