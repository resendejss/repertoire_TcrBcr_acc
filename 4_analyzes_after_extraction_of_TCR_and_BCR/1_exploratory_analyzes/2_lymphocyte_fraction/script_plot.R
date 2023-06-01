################################################################################
## ANALISE = fracao de linfocitos - plots
## PROJETO = repertorio de TCR e BCR na coorte TCGA-ACC
##
## CRIACAO DO SCRIPT (data) = 31/05/23
## ULTIMA ATUALIZACAO (data) = 31/05/23
## 
## RESPONSAVEL = Jean Resende
## Big Data
################################################################################
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

# -- importar metodos de fracoes celulares
cibersort <- read.csv("cibersort.csv")
epic <- read.csv("epic.csv")
estimate <- read.csv("estimate.csv")
ips <- read.csv("ips.csv")
mcp <- read.csv("mcp.csv")
quantiseq <- read.csv("quantiseq.csv")
timer <- read.csv("timer.csv")
xcell <- read.csv("xcell.csv")

# -- metadata
load("metada.RData")

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

# -- heatmap cibersort
data <- cibersort
rownames(data) <- data$ID
data <- data[,3:24]
data <- as.matrix(t(data))
#data <- log2(data + 1)

ht.cibersort <- Heatmap(data,
                        top_annotation = col.ha,
                        show_column_names = FALSE,
                        cluster_columns = FALSE,
                        cluster_rows = F,
                        col=colorRamp2(breaks = seq(-0,1, length.out=9),
                                       colors = brewer.pal(9,"YlOrRd")),
                        row_title_gp = gpar(fontsize=10),
                        row_title_side = "left",
                        row_names_side = "right",
                        row_names_gp = gpar(fontsize=10))

# -- heatmap epic
data <- epic
rownames(data) <- data$ID
data <- data[,3:ncol(data)]
data <- as.matrix(t(data))
#data <- log2(data + 1)

ht.epic <- Heatmap(data,
                   top_annotation = col.ha,
                   show_column_names = FALSE,
                   cluster_columns = FALSE,
                   cluster_rows = F,
                   col=colorRamp2(breaks = seq(-0,1, length.out=9),
                                  colors = brewer.pal(9,"YlOrRd")),
                   row_title_gp = gpar(fontsize=10),
                   row_title_side = "left",
                   row_names_side = "right",
                   row_names_gp = gpar(fontsize=10))

# -- estimate
data <- estimate
rownames(data) <- data$ID
data <- data[,3:ncol(data)]
data <- as.matrix(t(data))
# data <- log2(data + 1) --- parei aqui

ht.estimate <- Heatmap(data,
                       top_annotation = col.ha,
                       show_column_names = FALSE,
                       cluster_columns = FALSE,
                       cluster_rows = F,
                       col=colorRamp2(breaks = seq(-0,1, length.out=9),
                                      colors = rev(brewer.pal(9,"YlOrRd"))),
                       row_title_gp = gpar(fontsize=10),
                       row_title_side = "left",
                       row_names_side = "right",
                       row_names_gp = gpar(fontsize=10))

