################################################################################
# Name: expression_TCR_BCR_boxplot                                             #
#                                                                              #
# Author: Jean Resende (jean.s.s.resende@gmail.com)                            #
# Creation date: 2023/06/19                                                    #
# Last update date: 2023/06/20                                                 #
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
data_vdj_BT <- data_vdj
data_vdj_BT$b <- data_vdj_BT$IGK + data_vdj_BT$IGL + data_vdj_BT$IGH
data_vdj_BT$t <- data_vdj_BT$TRB + data_vdj_BT$TRA + data_vdj_BT$TRG + data_vdj_BT$TRD

data_vdj_BT$IGK[1] + data_vdj_BT$IGL[1] + data_vdj_BT$IGH[1]
data_vdj_BT$b[1]
data_vdj_BT$TRB[1] + data_vdj_BT$TRA[1] + data_vdj_BT$TRG[1] + data_vdj_BT$TRD[1]
data_vdj_BT$t[1]

data_vdj_BT$bnorm <- data_vdj_BT$b/metadata$reads[idx]
data_vdj_BT$tnorm <- data_vdj_BT$t/metadata$reads[idx]

data_vdj_log2 <- data_vdj_BT
data_vdj_log2[,9:12] <- log2(data_vdj_log2[,9:12]+ 1)

data_vdj_log2 <- na.omit(data_vdj_log2)

g1 <- ggplot(data_vdj_log2, aes(x=steroid, y=b, fill=steroid))+
  ggtitle("Non-normalized BCR")+
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
                     label.y = 17,
                     label="p.format",
                     size=3)

g2 <- ggplot(data_vdj_log2, aes(x=steroid, y=t, fill=steroid))+
  ggtitle("Non-normalized TCR")+
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
                     label.y = 17,
                     label="p.format",
                     size=3)

g3 <- ggplot(data_vdj_log2, aes(x=steroid, y=bnorm, fill=steroid))+
  ggtitle("Normalized BCR")+
  geom_boxplot(width=c(0.4), lwd=0.2, outlier.size = 0.5)+
  scale_y_continuous("expression BCR")+
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
        legend.box = "horizontal" )+g2
  stat_compare_means(method="wilcox.test", 
                     label.x = c(1.3), 
#                     label.y = 17,
                     label="p.format",
                     size=3)

g4 <- ggplot(data_vdj_log2, aes(x=steroid, y=tnorm, fill=steroid))+
  ggtitle("Normalized TCR")+
  geom_boxplot(width=c(0.4), lwd=0.2, outlier.size = 0.5)+
  scale_y_continuous("expression TCR")+
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
#                    label.y = 17,
                     label="p.format",
                     size=3)


plot_grid(g1,g2,g3,g4,ncol = 2, labels="AUTO")
rm(list = ls())

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

data_vdj[,1:7] <- log2(data_vdj[,1:7]+1) 

data_vdj$IGK <- data_vdj$IGK/metadata$reads[idx]
data_vdj$IGL <- data_vdj$IGL/metadata$reads[idx]
data_vdj$IGH <- data_vdj$IGH/metadata$reads[idx]
data_vdj$TRB <- data_vdj$TRB/metadata$reads[idx]
data_vdj$TRA <- data_vdj$TRA/metadata$reads[idx]
data_vdj$TRG <- data_vdj$TRG/metadata$reads[idx]
data_vdj$TRD <- data_vdj$TRD/metadata$reads[idx]

data_vdj <- na.omit(data_vdj)

g1 <- ggplot(data_vdj, aes(x=steroid, y=IGK, fill=steroid))+
  ggtitle("Normalized IGK")+
  geom_boxplot(width=c(0.4), lwd=0.2, outlier.size = 0.5)+
  scale_y_continuous("Expression(log2)")+
  scale_x_discrete("")+
  scale_fill_manual(values=c("#F8766D","#00BFC4"), 
                    labels=c("HSP (n=45)","LSP (n=30)"))+
  geom_rug(data = data_vdj[data_vdj$steroid=="LSP",], 
           aes(x=NULL), 
           sides=  "r",
           colour= "#00BFC4")+
  geom_rug(data = data_vdj[data_vdj$steroid=="HSP",], 
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
                     #label.y = 0.0005,
                     label="p.format",
                     size=3)

g2 <- ggplot(data_vdj, aes(x=steroid, y=IGL, fill=steroid))+
  ggtitle("Normalized IGL")+
  geom_boxplot(width=c(0.4), lwd=0.2, outlier.size = 0.5)+
  scale_y_continuous("Expression(log2)")+
  scale_x_discrete("")+
  scale_fill_manual(values=c("#F8766D","#00BFC4"), 
                    labels=c("HSP (n=45)","LSP (n=30)"))+
  geom_rug(data = data_vdj[data_vdj$steroid=="LSP",], 
           aes(x=NULL), 
           sides=  "r",
           colour= "#00BFC4")+
  geom_rug(data = data_vdj[data_vdj$steroid=="HSP",], 
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
                     #label.y = 0.0005,
                     label="p.format",
                     size=3)

g3 <- ggplot(data_vdj, aes(x=steroid, y=IGH, fill=steroid))+
  ggtitle("Normalized IGH")+
  geom_boxplot(width=c(0.4), lwd=0.2, outlier.size = 0.5)+
  scale_y_continuous("Expression(log2)")+
  scale_x_discrete("")+
  scale_fill_manual(values=c("#F8766D","#00BFC4"), 
                    labels=c("HSP (n=45)","LSP (n=30)"))+
  geom_rug(data = data_vdj[data_vdj$steroid=="LSP",], 
           aes(x=NULL), 
           sides=  "r",
           colour= "#00BFC4")+
  geom_rug(data = data_vdj[data_vdj$steroid=="HSP",], 
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
                     #label.y = 0.0005,
                     label="p.format",
                     size=3)

g4 <- ggplot(data_vdj, aes(x=steroid, y=TRA, fill=steroid))+
  ggtitle("Normalized TRA")+
  geom_boxplot(width=c(0.4), lwd=0.2, outlier.size = 0.5)+
  scale_y_continuous("Expression(log2)")+
  scale_x_discrete("")+
  scale_fill_manual(values=c("#F8766D","#00BFC4"), 
                    labels=c("HSP (n=45)","LSP (n=30)"))+
  geom_rug(data = data_vdj[data_vdj$steroid=="LSP",], 
           aes(x=NULL), 
           sides=  "r",
           colour= "#00BFC4")+
  geom_rug(data = data_vdj[data_vdj$steroid=="HSP",], 
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
                     #label.y = 0.0005,
                     label="p.format",
                     size=3)

g5 <- ggplot(data_vdj, aes(x=steroid, y=TRB, fill=steroid))+
  ggtitle("Normalized TRB")+
  geom_boxplot(width=c(0.4), lwd=0.2, outlier.size = 0.5)+
  scale_y_continuous("Expression(log2)")+
  scale_x_discrete("")+
  scale_fill_manual(values=c("#F8766D","#00BFC4"), 
                    labels=c("HSP (n=45)","LSP (n=30)"))+
  geom_rug(data = data_vdj[data_vdj$steroid=="LSP",], 
           aes(x=NULL), 
           sides=  "r",
           colour= "#00BFC4")+
  geom_rug(data = data_vdj[data_vdj$steroid=="HSP",], 
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
                     #label.y = 0.0005,
                     label="p.format",
                     size=3)

g6 <- ggplot(data_vdj, aes(x=steroid, y=TRG, fill=steroid))+
  ggtitle("Normalized TRG")+
  geom_boxplot(width=c(0.4), lwd=0.2, outlier.size = 0.5)+
  scale_y_continuous("Expression(log2)")+
  scale_x_discrete("")+
  scale_fill_manual(values=c("#F8766D","#00BFC4"), 
                    labels=c("HSP (n=45)","LSP (n=30)"))+
  geom_rug(data = data_vdj[data_vdj$steroid=="LSP",], 
           aes(x=NULL), 
           sides=  "r",
           colour= "#00BFC4")+
  geom_rug(data = data_vdj[data_vdj$steroid=="HSP",], 
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
                     #label.y = 0.0005,
                     label="p.format",
                     size=3)

g7 <- ggplot(data_vdj, aes(x=steroid, y=TRD, fill=steroid))+
  ggtitle("Normalized TRD")+
  geom_boxplot(width=c(0.4), lwd=0.2, outlier.size = 0.5)+
  scale_y_continuous("Expression(log2)")+
  scale_x_discrete("")+
  scale_fill_manual(values=c("#F8766D","#00BFC4"), 
                    labels=c("HSP (n=45)","LSP (n=30)"))+
  geom_rug(data = data_vdj[data_vdj$steroid=="LSP",], 
           aes(x=NULL), 
           sides=  "r",
           colour= "#00BFC4")+
  geom_rug(data = data_vdj[data_vdj$steroid=="HSP",], 
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
                     #label.y = 0.0005,
                     label="p.format",
                     size=3)

plot_grid(g1,g2,g3,g4,g5,g6,g7,ncol = 4, labels="AUTO")
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
#data_vdj_BT <- data_vdj
#data_vdj_BT$b <- data_vdj_BT$IGK + data_vdj_BT$IGL + data_vdj_BT$IGH
#data_vdj_BT$t <- data_vdj_BT$TRB + data_vdj_BT$TRA + data_vdj_BT$TRG + data_vdj_BT$TRD

#data_vdj_BT$IGK[1] + data_vdj_BT$IGL[1] + data_vdj_BT$IGH[1]
#data_vdj_BT$b[1]
#data_vdj_BT$TRB[1] + data_vdj_BT$TRA[1] + data_vdj_BT$TRG[1] + data_vdj_BT$TRD[1]
#data_vdj_BT$t[1]

#data_vdj[,1:7] <- log2(data_vdj[,1:7]+1) 

data_vdj$IGK <- data_vdj$IGK/metadata$reads[idx]
data_vdj$IGL <- data_vdj$IGL/metadata$reads[idx]
data_vdj$IGH <- data_vdj$IGH/metadata$reads[idx]
data_vdj$TRB <- data_vdj$TRB/metadata$reads[idx]
data_vdj$TRA <- data_vdj$TRA/metadata$reads[idx]
data_vdj$TRG <- data_vdj$TRG/metadata$reads[idx]
data_vdj$TRD <- data_vdj$TRD/metadata$reads[idx]

data_vdj_log10 <- data_vdj
data_vdj_log10[,1:7] <- log10(data_vdj_log10[,1:7]+1)

data_vdj_log2 <- data_vdj
data_vdj_log2[,1:7] <- log2(data_vdj[,1:7]+1)

data_vdj_log2 <- na.omit(data_vdj_log2)
data_vdj_log10 <- na.omit(data_vdj_log10)

#data_vdj <- na.omit(data_vdj)

g2 <- ggplot(data_vdj_log10, aes(x=steroid, y=IGK, fill=steroid))+
  ggtitle("Normalized IGK")+
  geom_boxplot(width=c(0.4), lwd=0.2, outlier.size = 0.5)+
  scale_y_continuous("Expression(log10)")+
  scale_x_discrete("")+
  scale_fill_manual(values=c("#F8766D","#00BFC4"), 
                    labels=c("HSP (n=45)","LSP (n=30)"))+
  geom_rug(data = data_vdj_log10[data_vdj_log10$steroid=="LSP",], 
           aes(x=NULL), 
           sides=  "r",
           colour= "#00BFC4")+
  geom_rug(data = data_vdj_log10[data_vdj_log10$steroid=="HSP",], 
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
                     #label.y = 0.0005,
                     label="p.format",
                     size=3)
