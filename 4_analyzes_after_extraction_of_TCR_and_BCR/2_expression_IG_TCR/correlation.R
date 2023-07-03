################################################################################
## ANALISE = correlation
## PROJETO = repertorio de TCR e BCR na coorte TCGA-ACC
##
## CRIACAO DO SCRIPT (data) = 31/05/03
## ULTIMA ATUALIZACAO (data) = 20/06/2023
## 
## RESPONSAVEL = Jean Resende
## Big Data
################################################################################
library(corrplot)

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

lsp_counts <- data_vdj[data_vdj$steroid=="LSP",1:7]
lsp_counts <- lsp_counts[!is.na(lsp_counts$IGK),]

hsp_counts <- data_vdj[data_vdj$steroid=="HSP",1:7]
hsp_counts <- hsp_counts[!is.na(hsp_counts$IGK),]



