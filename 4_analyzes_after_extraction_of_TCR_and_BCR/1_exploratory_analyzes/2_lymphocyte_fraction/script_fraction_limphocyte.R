################################################################################
## ANALISE = fracao de linfocitos
## PROJETO = repertorio de TCR e BCR na coorte TCGA-ACC
##
## CRIACAO DO SCRIPT (data) = 09/05/23
## ULTIMA ATUALIZACAO (data) = 13/06/23
## 
## RESPONSAVEL = Jean Resende
## Big Data
################################################################################
# -- instalacao do IOBR -- #####################################################

## -- instalacao de dependencias
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

depens<-c('tibble', 'survival', 'survminer', 'sva', 'limma', "DESeq2","devtools",
          'limSolve', 'GSVA', 'e1071', 'preprocessCore', 'ggplot2', "biomaRt",
          'ggpubr', "devtools", "tidyHeatmap", "caret", "glmnet", "ppcor", "timeROC","pracma")

for(i in 1:length(depens)){
  depen<-depens[i]
  if (!requireNamespace(depen, quietly = TRUE))
    BiocManager::install(depen,update = FALSE)
}


## -- instalacao do IOBR
if (!requireNamespace("remotes", quietly = TRUE))
  install("remotes")

if (!requireNamespace("IOBR", quietly = TRUE))
  remotes::install_github("IOBR/IOBR",ref="master")

## -- instalacao do UCSCXenaTools
if (!requireNamespace("UCSCXenaTools", quietly = TRUE))
  BiocManager::install("UCSCXenaTools")

# -- 
library(UCSCXenaTools)
library(IOBR)
library(tidyr)

## -- download da matriz de expressao genica atraves do xenabrowser
eset_acc<-XenaGenerate(subset = XenaCohorts =="GDC TCGA Adrenocortical Cancer (ACC)") %>% 
  XenaFilter(filterDatasets    = "TCGA-ACC.htseq_counts.tsv") %>% 
  XenaQuery() %>%
  XenaDownload() %>% 
  XenaPrepare()

eset_acc[1:5,1:5]

## -- anotacao e snormalizacao em TPM
eset_acc$Ensembl_ID<-substring(eset_acc$Ensembl_ID, 1, 15) # remocao da versao
eset_acc<-column_to_rownames(eset_acc, var = "Ensembl_ID")
eset_acc <- (2^eset_acc)+1 # revertendo para a contagem original (pois estao em log2)
eset_acc<-count2tpm(countMat = eset_acc, idType = "Ensembl") # normalizacao em TPM
eset_acc[1:5,1:5]

## -- ordenacao conforme low e high steroid
load("metadata.RData")

idx <- match(gsub("-..R-A29S-07","",metadata$barcode),colnames(eset_acc))
idx <- idx[!is.na(idx)]
eset_acc <- eset_acc[,idx]
#idx2 <- match(colnames(eset_acc),substring(metadata$barcode, 1, 16))

save(eset_acc, file = "eset_acc_20230607.RData")

## -- assinaturas associadas ao TME - metodo PCA
load("eset_acc_20230607.RData")
sig_tme <- calculate_sig_score(pdata = NULL,
                               eset = eset_acc,
                               signature = signature_tme,
                               method = "pca",
                               mini_gene_count = 2)
colnames(sig_tme)

sig_hallmark <- calculate_sig_score(pdata = NULL,
                                    eset = eset_acc,
                                    signature = hallmark,
                                    method = "ssgsea",
                                    mini_gene_count = 2)
colnames(sig_hallmark)

## -- fracao celular 
### -- cibersort
load("eset_acc.RData")
eset_acc_tpm_log2 <- log2(eset_acc+1) # transformacao em log2

cibersort <- deconvo_tme(eset = eset_acc_tpm_log2, method = "cibersort", arrays = F, perm = 200)

#res<-cell_bar_plot(input = cibersort, title = "CIBERSORT Cell Fraction")

#pdf(file = "cellBarPlot_cibersort.pdf", width = 8.27, height = 11.69)
#cell_bar_plot(input = cibersort, title = "CIBERSORT Cell Fraction")
#dev.off()

write.csv(cibersort, file = "cibersort.csv")

### -- quantiseq
load("eset_acc_20230607.RData")
eset_acc_tpm_log2 <- log2(eset_acc+1) # transformacao em log2

quantiseq <- deconvo_tme(eset = eset_acc_tpm_log2,
                         tumor = TRUE, arrays = FALSE, scale_mrna = TRUE,
                         method = "quantiseq")
#res <- cell_bar_plot(input = quantiseq, title = "quanTIseq Cell Fraction")

#pdf(file = "cellBarPlot_quantiseq.pdf", width = 8.27, height = 11.69)
#cell_bar_plot(input = quantiseq, title = "quanTIseq Cell Fraction")
#dev.off()

write.csv(quantiseq, file = "quantiseq.csv")

################################################################################
### -- epic
load("eset_acc_20230607.RData")
eset_acc_tpm_log2 <- log2(eset_acc+1) # transformacao em log2

epic <- deconvo_tme(eset = eset_acc_tpm_log2, method = "epic", arrays = FALSE)

rowSums(epic[,-1])

#pdf(file = "cellBarPlot_epic.pdf", width = 8.27, height = 11.69)
#cell_bar_plot(input = epic, title = "EPIC Cell Fraction")
#dev.off()

write.csv(epic, file = "epic.csv")

### -- mcpcounter
load("eset_acc_20230607.RData")
eset_acc_tpm_log2 <- log2(eset_acc+1) # transformacao em log2

mcp <- deconvo_tme(eset = eset_acc_tpm_log2, method = "mcpcounter")
rowSums(mcp[,-1])
summary(mcp)

write.csv(mcp, file = "mcp.csv")

### -- xcell
load("eset_acc_20230607.RData")
eset_acc_tpm_log2 <- log2(eset_acc+1) # transformacao em log2
xcell <- deconvo_tme(eset_acc_tpm_log2, method = "xcell", arrays = FALSE)

write.csv(xcell, file = "xcell.csv")

### -- estimate
load("eset_acc.RData")
eset_acc_tpm_log2 <- log2(eset_acc+1) # transformacao em log2
estimate <- deconvo_tme(eset = eset_acc_tpm_log2, method = "estimate")

write.csv(estimate, file = "estimate.csv")

### -- timer
load("eset_acc.RData")
eset_acc_tpm_log2 <- log2(eset_acc+1) # transformacao em log2
timer <- deconvo_tme(eset = eset_acc_tpm_log2, method = "timer", group_list = rep("acc",dim(eset_acc_tpm_log2)[2]))

write.csv(timer, file = "timer.csv")

### -- ips
load("eset_acc.RData")
eset_acc_tpm_log2 <- log2(eset_acc+1) # transformacao em log2
ips <- deconvo_tme(eset = eset_acc_tpm_log2, method = "ips", plot= FALSE)

write.csv(ips, file = "ips.csv")



