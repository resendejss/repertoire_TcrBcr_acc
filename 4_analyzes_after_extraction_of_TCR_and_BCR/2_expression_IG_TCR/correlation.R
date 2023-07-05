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

# -- normalizacao
vdj_norm <- data_vdj
#vdj_norm[,1:7] <- log2(vdj_norm[,1:7]+1)
identical(rownames(vdj_norm),metadata$sample_id)
idx <- match(metadata$sample_id, rownames(vdj_norm))
vdj_norm[,1:7] <- vdj_norm[,1:7]/metadata$reads[idx]

lsp_expression <- vdj_norm[vdj_norm$steroid=="LSP",1:7]
lsp_expression <- lsp_expression[!is.na(lsp_expression$IGK),]

hsp_expression <- vdj_norm[vdj_norm$steroid=="HSP",1:7]
hsp_expression <- hsp_expression[!is.na(hsp_expression$IGK),]

# -- entropy
load("../afterExtraction_trust4/data.entropy.RData")

idx <- match(gsub("_report","",rownames(data.entropy)), metadata$sample_id)
data.entropy$steroid <- metadata$steroid[idx]

data.entropy$IGH[is.na(data.entropy$IGH)] <- 0
data.entropy$IGL[is.na(data.entropy$IGL)] <- 0
data.entropy$IGK[is.na(data.entropy$IGK)] <- 0
data.entropy$TRA[is.na(data.entropy$TRA)] <- 0
data.entropy$TRB[is.na(data.entropy$TRB)] <- 0
data.entropy$TRD[is.na(data.entropy$TRD)] <- 0
data.entropy$TRG[is.na(data.entropy$TRG)] <- 0

lsp_entropy <- data.entropy[data.entropy$steroid=="LSP",c(1,11:16)]
lsp_entropy <- lsp_entropy[!is.na(lsp_entropy$IGK),]

hsp_entropy <- data.entropy[data.entropy$steroid=="HSP",c(1,11:16)]
hsp_entropy <- hsp_entropy[!is.na(hsp_entropy$IGK),]


colnames(lsp_counts) <- paste(colnames(lsp_counts),"counts", sep = "_")
colnames(hsp_counts) <- paste(colnames(hsp_counts),"counts", sep = "_")
colnames(lsp_expression) <- paste(colnames(lsp_expression),"expression", sep = "_")
colnames(hsp_expression) <- paste(colnames(hsp_expression),"expression", sep = "_")
colnames(lsp_entropy) <- paste(colnames(lsp_entropy),"entropy", sep = "_")
colnames(hsp_entropy) <- paste(colnames(hsp_entropy),"entropy", sep = "_")


rownames(lsp_entropy) <- gsub("_report","",rownames(lsp_entropy))
rownames(hsp_entropy) <- gsub("_report","",rownames(hsp_entropy))

table(rownames(lsp_counts) == rownames(lsp_expression))
table(rownames(lsp_counts) == rownames(lsp_entropy))
idx <- match(rownames(lsp_counts), rownames(lsp_entropy))
lsp_entropy <- lsp_entropy[idx,]
table(rownames(lsp_counts) == rownames(lsp_entropy))

table(rownames(hsp_counts) == rownames(hsp_expression))
table(rownames(hsp_counts) == rownames(hsp_entropy))
idx <- match(rownames(hsp_counts), rownames(hsp_entropy))
hsp_entropy <- hsp_entropy[idx,]
table(rownames(hsp_counts) == rownames(hsp_entropy))

gsub("_counts","",colnames(lsp_counts)) == gsub("_expression","",colnames(lsp_expression))
gsub("_counts","",colnames(lsp_counts)) == gsub("_entropy","",colnames(lsp_entropy))
lsp_entropy <- lsp_entropy[,c(2,3,1,5,4,6,7)]
gsub("_counts","",colnames(lsp_counts)) == gsub("_entropy","",colnames(lsp_entropy))

gsub("_counts","",colnames(hsp_counts)) == gsub("_expression","",colnames(hsp_expression))
gsub("_counts","",colnames(hsp_counts)) == gsub("_entropy","",colnames(hsp_entropy))
hsp_entropy <- hsp_entropy[,c(2,3,1,5,4,6,7)]
gsub("_counts","",colnames(hsp_counts)) == gsub("_entropy","",colnames(hsp_entropy))

lsp <- cbind(lsp_counts, lsp_expression, lsp_entropy)
hsp <- cbind(hsp_counts, hsp_expression, hsp_entropy)


lsp <- lsp[,colSums(lsp) != 0]
hsp <- hsp[,colSums(hsp) != 0]

res_lsp <- cor(lsp)
res_hsp <- cor(hsp)

teste <- cor.mtest(lsp, conf.level=0.95)
corrplot(res_lsp, p.mat =  teste$p, method = "color", order = "AOE", type = "upper", 
         diag = F, tl.col = "black", insig = "label_sig", pch.cex = 0.8)

teste <- cor.mtest(hsp, conf.level=0.95)
corrplot(res_hsp, p.mat =  teste$p, method = "color", order = "AOE", type = "upper", 
         diag = F, tl.col = "black", insig = "label_sig", pch.cex = 0.8)
 