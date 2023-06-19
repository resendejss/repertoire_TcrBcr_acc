################################################################################
# Name: expression_TCR_BCR_boxplot                                             #
#                                                                              #
# Author: Jean Resende (jean.s.s.resende@gmail.com)                            #
# Creation date: 2023/06/19                                                    #
# Last update date: 2023/06/19                                                 #
################################################################################
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
data_vdj_log2 <- data_vdj
data_vdj_log2[,1:7] <- log2(data_vdj_log2[,1:7]+ 1)

data_vdj_log2 <- na.omit(data_vdj_log2)


g1 <- ggplot(data_vdj_log2, aes(x=steroid, y=IGK, fill=steroid))+
  ggtitle("IGK")+
  geom_boxplot(width=c(0.4), lwd=0.2, outlier.size = 0.5)+
  scale_y_continuous("IGL counts")+
  scale_x_discrete("")+
  scale_fill_manual(values=c("#F8766D","#00BFC4"), 
                    labels=c("HSP (n=45)","LSP (n=30)"))+
  geom_rug(data = data_vdj_log2[data_vdj_log2$steroid=="LSP",], 
           aes(x=NULL), 
           sides=  "r",
           colour= "#00BFC4")+
  geom_rug(data = data_vdj_log2[data_vdj_log2$steroid=="HSP",], 
         aes(x=NULL), 
         sides=  "l",
         colour="#F8766D")+
  theme(legend.position = "bottom", 
        legend.title = element_text(face="bold", size=10),
        legend.key=element_rect(size=10, color=NA), 
        legend.key.size=unit(8,"mm"),
        legend.text=element_text(size=10), 
        legend.direction = "horizontal",
        legend.box = "horizontal" )+
  stat_compare_means(method="wilcox.test", 
                     label.x = c(1.3), 
                     label.y = 16.1,
                     label="p.format",
                     size=3)

g2 <- ggplot(data_vdj_log2, aes(x=steroid, y=IGL, fill=steroid))+
  ggtitle("IGL")+
  geom_boxplot(width=c(0.4), lwd=0.2, outlier.size = 0.5)+
  scale_y_continuous("IGL counts")+
  scale_x_discrete("")+
  scale_fill_manual(values=c("#F8766D","#00BFC4"), 
                    labels=c("HSP (n=45)","LSP (n=30)"))+
  geom_rug(data = data_vdj_log2[data_vdj_log2$steroid=="LSP",], 
           aes(x=NULL), 
           sides=  "r",
           colour= "#00BFC4")+
  geom_rug(data = data_vdj_log2[data_vdj_log2$steroid=="HSP",], 
           aes(x=NULL), 
           sides=  "l",
           colour="#F8766D")+
  theme(legend.position = "bottom", 
        legend.title = element_text(face="bold", size=10),
        legend.key=element_rect(size=10, color=NA), 
        legend.key.size=unit(8,"mm"),
        legend.text=element_text(size=10), 
        legend.direction = "horizontal",
        legend.box = "horizontal" )+
  stat_compare_means(method="wilcox.test", 
                     label.x = c(1.3), 
                     label.y = 16.1,
                     label="p.format",
                     size=3)

g3 <- ggplot(data_vdj_log2, aes(x=steroid, y=IGH, fill=steroid))+
  ggtitle("IGH")+
  geom_boxplot(width=c(0.4), lwd=0.2, outlier.size = 0.5)+
  scale_y_continuous("IGH counts")+
  scale_x_discrete("")+
  scale_fill_manual(values=c("#F8766D","#00BFC4"), 
                    labels=c("HSP (n=45)","LSP (n=30)"))+
  geom_rug(data = data_vdj_log2[data_vdj_log2$steroid=="LSP",], 
           aes(x=NULL), 
           sides=  "r",
           colour= "#00BFC4")+
  geom_rug(data = data_vdj_log2[data_vdj_log2$steroid=="HSP",], 
           aes(x=NULL), 
           sides=  "l",
           colour="#F8766D")+
  theme(legend.position = "bottom", 
        legend.title = element_text(face="bold", size=10),
        legend.key=element_rect(size=10, color=NA), 
        legend.key.size=unit(8,"mm"),
        legend.text=element_text(size=10), 
        legend.direction = "horizontal",
        legend.box = "horizontal" )+
  stat_compare_means(method="wilcox.test", 
                     label.x = c(1.3), 
                     label.y = 16.1,
                     label="p.format",
                     size=3)

g4 <- ggplot(data_vdj_log2, aes(x=steroid, y=TRB, fill=steroid))+
  ggtitle("TRB")+
  geom_boxplot(width=c(0.4), lwd=0.2, outlier.size = 0.5)+
  scale_y_continuous("TRB counts")+
  scale_x_discrete("")+
  scale_fill_manual(values=c("#F8766D","#00BFC4"), 
                    labels=c("HSP (n=45)","LSP (n=30)"))+
  geom_rug(data = data_vdj_log2[data_vdj_log2$steroid=="LSP",], 
           aes(x=NULL), 
           sides=  "r",
           colour= "#00BFC4")+
  geom_rug(data = data_vdj_log2[data_vdj_log2$steroid=="HSP",], 
           aes(x=NULL), 
           sides=  "l",
           colour="#F8766D")+
  theme(legend.position = "bottom", 
        legend.title = element_text(face="bold", size=10),
        legend.key=element_rect(size=10, color=NA), 
        legend.key.size=unit(8,"mm"),
        legend.text=element_text(size=10), 
        legend.direction = "horizontal",
        legend.box = "horizontal" )+
  stat_compare_means(method="wilcox.test", 
                     label.x = c(1.3), 
                     label.y = 16.1,
                     label="p.format",
                     size=3)

g5 <- ggplot(data_vdj_log2, aes(x=steroid, y=TRA, fill=steroid))+
  ggtitle("TRA")+
  geom_boxplot(width=c(0.4), lwd=0.2, outlier.size = 0.5)+
  scale_y_continuous("TRA counts")+
  scale_x_discrete("")+
  scale_fill_manual(values=c("#F8766D","#00BFC4"), 
                    labels=c("HSP (n=45)","LSP (n=30)"))+
  geom_rug(data = data_vdj_log2[data_vdj_log2$steroid=="LSP",], 
           aes(x=NULL), 
           sides=  "r",
           colour= "#00BFC4")+
  geom_rug(data = data_vdj_log2[data_vdj_log2$steroid=="HSP",], 
           aes(x=NULL), 
           sides=  "l",
           colour="#F8766D")+
  theme(legend.position = "bottom", 
        legend.title = element_text(face="bold", size=10),
        legend.key=element_rect(size=10, color=NA), 
        legend.key.size=unit(8,"mm"),
        legend.text=element_text(size=10), 
        legend.direction = "horizontal",
        legend.box = "horizontal" )+
  stat_compare_means(method="wilcox.test", 
                     label.x = c(1.3), 
                     label.y = 16.1,
                     label="p.format",
                     size=3)

g6 <- ggplot(data_vdj_log2, aes(x=steroid, y=TRG, fill=steroid))+
  ggtitle("TRG")+
  geom_boxplot(width=c(0.4), lwd=0.2, outlier.size = 0.5)+
  scale_y_continuous("TRG counts")+
  scale_x_discrete("")+
  scale_fill_manual(values=c("#F8766D","#00BFC4"), 
                    labels=c("HSP (n=45)","LSP (n=30)"))+
  geom_rug(data = data_vdj_log2[data_vdj_log2$steroid=="LSP",], 
           aes(x=NULL), 
           sides=  "r",
           colour= "#00BFC4")+
  geom_rug(data = data_vdj_log2[data_vdj_log2$steroid=="HSP",], 
           aes(x=NULL), 
           sides=  "l",
           colour="#F8766D")+
  theme(legend.position = "bottom", 
        legend.title = element_text(face="bold", size=10),
        legend.key=element_rect(size=10, color=NA), 
        legend.key.size=unit(8,"mm"),
        legend.text=element_text(size=10), 
        legend.direction = "horizontal",
        legend.box = "horizontal" )+
  stat_compare_means(method="wilcox.test", 
                     label.x = c(1.3), 
                     label.y = 16.1,
                     label="p.format",
                     size=3)

g7 <- ggplot(data_vdj_log2, aes(x=steroid, y=TRD, fill=steroid))+
  ggtitle("TRD")+
  geom_boxplot(width=c(0.4), lwd=0.2, outlier.size = 0.5)+
  scale_y_continuous("TRD counts")+
  scale_x_discrete("")+
  scale_fill_manual(values=c("#F8766D","#00BFC4"), 
                    labels=c("HSP (n=45)","LSP (n=30)"))+
  geom_rug(data = data_vdj_log2[data_vdj_log2$steroid=="LSP",], 
           aes(x=NULL), 
           sides=  "r",
           colour= "#00BFC4")+
  geom_rug(data = data_vdj_log2[data_vdj_log2$steroid=="HSP",], 
           aes(x=NULL), 
           sides=  "l",
           colour="#F8766D")+
  theme(legend.position = "bottom", 
        legend.title = element_text(face="bold", size=10),
        legend.key=element_rect(size=10, color=NA), 
        legend.key.size=unit(8,"mm"),
        legend.text=element_text(size=10), 
        legend.direction = "horizontal",
        legend.box = "horizontal" )+
  stat_compare_means(method="wilcox.test", 
                     label.x = c(1.3), 
                     label.y = 16.1,
                     label="p.format",
                     size=3)

plot_grid(g1,g2,g3,g4,g5,g6,g7,ncol = 4, labels="AUTO")

