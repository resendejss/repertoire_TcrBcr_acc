################################################################################
## ANALISE = fracao de linfocitos - plots
## PROJETO = repertorio de TCR e BCR na coorte TCGA-ACC
##
## CRIACAO DO SCRIPT (data) = 31/05/23
## ULTIMA ATUALIZACAO (data) = 20/06/2023
## 
## RESPONSAVEL = Jean Resende
## Big Data
################################################################################
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(magrittr)

# -- importar metodos de fracoes celulares
#cibersort <- read.csv("cibersort.csv")
#estimate <- read.csv("estimate.csv")
#ips <- read.csv("ips.csv")
#timer <- read.csv("timer.csv")

# -- metadata
load("metadata.RData")

col.seqGreys.reads <- colorRamp2(breaks = seq(21000000,83000000,length.out=9),
                                 colors = brewer.pal(9,"Greys"))

col.seqGreys <- colorRamp2(breaks = seq(0,10,length.out=9),
                           colors = brewer.pal(9, "Greys"))

col.bin.steroid <- c("HSP"="#F8766D", "LSP"="#00BFC4")
col.bin.vital <- c("Alive"="white", "Dead"="black")
col.bin.gender <- c("female"="white","male"="black")
col.bin.stage <- c("stage i"="#d9d9d9", "stage ii"="#969696",
                   "stage iii"="#525252", "stage iv"="#000000",
                   "not reported"="#ffffff")
col.bin.cortisol <- c("Yes"="black", "No"="white")
col.bin.imm_sub <- c("C1"="#8dd3c7","C2"="#ffffb3","C3"="#bebada",
                     "C4"="#b3de69","C5"="#80b1d3","C6"="#fdb462")
col.bin.oth_horm <- c("Mineralcorticoids"="#e41a1c","Sexual"="#377eb8",
                      "No"="#bdbdbd")

col.NR5A1 <- colorRamp2(breaks = seq(-2,2,length.out=9),
                        colors = rev(brewer.pal(9,"RdBu")))

col.ha <- HeatmapAnnotation(steroid=metadata$steroid,
                            NR5A1 = metadata$NR5A1,
                            immune_subtype=metadata$Immune.Subtype,
                            other.hormones=metadata$Other.Hormones,
                            stage=metadata$stage,
                            gender=metadata$gender,
                            vital=metadata$Vital.Status,
                            cortisol.excess=metadata$cortisol.excess,
                            count_TCR=metadata$tcr,
                            count_BCR=metadata$bcr,
                            qtd_reads=metadata$reads,
                            col=list(steroid=col.bin.steroid,
                                     NR5A1=col.NR5A1,
                                     stage=col.bin.stage,
                                     gender=col.bin.gender,
                                     vital=col.bin.vital,
                                     immune_subtype=col.bin.imm_sub,
                                     cortisol.excess=col.bin.cortisol,
                                     other.hormones=col.bin.oth_horm,
                                     count_TCR=col.seqGreys,
                                     count_BCR=col.seqGreys,
                                     qtd_reads=col.seqGreys.reads))

rm(col.bin.cortisol,col.bin.gender,col.bin.imm_sub,col.bin.oth_horm,
   col.bin.stage,col.bin.steroid,col.bin.vital,col.NR5A1,col.seqGreys,
   col.seqGreys.reads)

# -- heatmap epic
load("res_epic.RData")
rowSums(res.epic)
data <- res.epic
identical(rownames(data),gsub("-..R-A29S-07","",metadata$barcode))
data <- as.matrix(t(data))

range(data)

ht.epic <- Heatmap(data,
                   top_annotation = col.ha,
                   show_column_names = FALSE,
                   cluster_columns = FALSE,
                   cluster_rows = F,
                   col=colorRamp2(breaks = seq(0,1, length.out=9),
                                  colors = brewer.pal(9,"YlOrRd")),
                   row_title_gp = gpar(fontsize=10),
                   row_title_side = "left",
                   row_names_side = "right",
                   row_names_gp = gpar(fontsize=10))

head(colnames(data))
head(metadata$barcode)

# -- heatmap quantiseq
quantiseq <- read.csv("quantiseq.csv")
rowSums(quantiseq[,3:ncol(quantiseq)])
data <- quantiseq
rownames(data) <- data$ID
data <- data[,3:ncol(data)]
data <- as.matrix(t(data))
data <- data[-nrow(data),]
#data <- log2(data + 1)
range(data)

ht.quantiseq <- Heatmap(data,
                        #top_annotation = col.ha,
                        show_column_names = FALSE,
                        cluster_columns = FALSE,
                        cluster_rows = F,
                        col=colorRamp2(breaks = seq(0,0.3, length.out=9),
                                       colors = brewer.pal(9,"OrRd")),
                        row_title_gp = gpar(fontsize=10),
                        row_title_side = "left",
                        row_names_side = "right",
                        row_names_gp = gpar(fontsize=10))

ht_list = ht.epic %v% ht.quantiseq

# -- heatmap mcpcounter
mcp <- read.csv("mcp.csv")
rowSums(mcp[,3:ncol(mcp)])
data <- mcp
rownames(data) <- data$ID
data <- data[,c(3,4,6)]
data <- as.matrix(t(data))

range(data)
data <- log2(data)
data.z <- scale(data)
summary(data.z)
range(data)

ht.mcpcounter <- Heatmap(data.z,
                         #top_annotation = col.ha,
                         show_column_names = FALSE,
                         cluster_columns = F,
                         cluster_rows = F,
                         col=colorRamp2(breaks = seq(-2,2, length.out=9),
                                        colors = rev(brewer.pal(9,"RdYlBu"))),
                         row_title_gp = gpar(fontsize=10),
                         row_title_side = "left",
                         row_names_side = "right",
                         row_names_gp = gpar(fontsize=10))

head(colnames(data.z))
head(metadata$barcode)

ht_list = ht_list %v% ht.mcpcounter

# -- xcell
xcell <- read.csv("xcell.csv")
rowSums(xcell[,3:ncol(xcell)])
data <- xcell
rownames(data) <- data$ID
data <- data[,c(6,8:16,19,40,48)]
data <- as.matrix(t(data))

df.xcell <- as.data.frame(data)
df.xcell$cell_type <- c("B",rep("T",9),rep("B",3))

range(data)
data <- log2(data + 1)
data.z <- scale(data)
summary(data.z)
range(data.z)

data.z[data.z > 2] <- 2
data.z[data.z <(-2)] <- (-2)

ht.xcell <- Heatmap(data.z,
                         #top_annotation = col.ha,
                    split = df.xcell$cell_type,
                         show_column_names = FALSE,
                         cluster_columns = F,
                         cluster_rows = F,
                         col=colorRamp2(breaks = seq(-2,2, length.out=9),
                                        colors = rev(brewer.pal(9,"RdBu"))),
                         row_title_gp = gpar(fontsize=10),
                         row_title_side = "left",
                         row_names_side = "right",
                         row_names_gp = gpar(fontsize=10))

head(colnames(data.z))
head(metadata$barcode)

ht_list = ht_list %v% ht.xcell
################################################################################

# -- TRUST4
load("v.counts.RData")
load("d.counts.RData")
load("j.counts.RData")

### -- regiao v --

## -- construcao do data.v
data.v <- v.counts %>% as.matrix() %>% t()
all(colnames(data.v)%in% metadata$sample_id)

data.v <- log2(data.v + 1)

## -- transformacao z-score
data.v <- as.matrix(t(scale(t(data.v))))
summary(t(data.v))
#teste_z <- (data[1,]-mean(data[1,]))/sd(data[1,])
#all(data[1,] == teste_z)

## -- ordenacao de data - low e high steroid e total de TCR e BCR
data.v <- data.v[,metadata$sample_id]
data.v <- as.matrix(data.v)

## -- heatmap
### -- dados para argumentos
df.v <- as.data.frame(data.v)
df.v$cell_type <- c(rep("B",3),rep("T",4))

data.v[data.v > 2] <- 2
data.v[data.v <(-2)] <- (-2)

ht.v <- Heatmap(data.v,
                top_annotation = col.ha,
                split = df.v$cell_type,
                show_column_names = FALSE,
                cluster_columns = FALSE,
                cluster_rows = F,
                col=colorRamp2(breaks = seq(-2,2, length.out=9),
                               colors = rev(brewer.pal(9,"BrBG"))),
                row_title = "Region V",
                row_title_gp = gpar(fontsize=6),
                row_title_side = "left",
                row_names_side = "right",
                row_names_gp = gpar(fontsize=8))

rm(data.v,df.v,v.counts)

### -- region D --

## -- construcao do data.d
data.d <- d.counts %>% as.matrix() %>% t()
all(colnames(data.d)%in% metadata$sample_id)

data.d <- log2(data.d + 1)

## -- transformacao z-score
data.d <- data.d[rowSums(data.d)!=0,]
data.d <- as.matrix(t(scale(t(data.d))))
#teste_z <- (data[1,]-mean(data[1,]))/sd(data[1,])
#all(data[1,] == teste_z)

## -- ordenacao de data - low e high steroid e total de TCR e BCR
data.d <- data.d[,metadata$sample_id]
data.d <- as.matrix(data.d)

## -- heatmap
### -- dados para argumentos
df.d <- as.data.frame(data.d)
df.d$cell_type <- c("B","T","T")

data.d[data.d > 2] <- 2
data.d[data.d <(-2)] <- (-2)

ht.d <- Heatmap(data.d,
                split = df.d$cell_type,
                show_column_names = FALSE,
                cluster_columns = FALSE,
                cluster_rows = F,
                col=colorRamp2(breaks = seq(-2,2, length.out=9),
                               colors = rev(brewer.pal(9,"BrBG"))),
                row_title = "Region D",
                row_title_gp = gpar(fontsize=6),
                row_title_side = "left",
                row_names_side = "right",
                row_names_gp = gpar(fontsize=8))

ht_list = ht.v %v% ht.d
draw(ht_list, merge_legends=TRUE,annotation_legend_side = "top")

rm(d.counts,data.d,df.d,ht.d,ht.v)

### -- region J --

## -- construcao do data.j
data.j <- j.counts %>% as.matrix() %>% t()
all(colnames(data.j)%in% metadata$sample_id)

data.j <- log2(data.j + 1)

## -- transformacao z-score
data.j <- as.matrix(t(scale(t(data.j))))
#teste_z <- (data[1,]-mean(data[1,]))/sd(data[1,])
#all(data[1,] == teste_z)

## -- ordenacao de data - low e high steroid e total de TCR e BCR
data.j <- data.j[,metadata$sample_id]
data.j <- as.matrix(data.j)

## -- heatmap
### -- dados para argumentos
df.j <- as.data.frame(data.j)
df.j$cell_type <- c(rep("B",3),rep("T",4))

data.j[data.j > 2] <- 2
data.j[data.j <(-2)] <- (-2)

ht.j <- Heatmap(data.j,
                split = df.j$cell_type,
                show_column_names = FALSE,
                cluster_columns = FALSE,
                cluster_rows = F,
                col=colorRamp2(breaks = seq(-2,2, length.out=9),
                               colors = rev(brewer.pal(9,"BrBG"))),
                row_title = "Region J",
                row_title_gp = gpar(fontsize=6),
                row_title_side = "left",
                row_names_side = "right",
                row_names_gp = gpar(fontsize=8))

ht_list = ht_list %v% ht.j
draw(ht_list, merge_legends=TRUE,annotation_legend_side = "top")

rm(df.j, ht.j, j.counts)

ht_list2 = ht_list %v% ht.xcell
draw(ht_list2, merge_legends=TRUE,annotation_legend_side = "top")
