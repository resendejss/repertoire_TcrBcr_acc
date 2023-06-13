################################################################################
## ANALISE = fracao de linfocitos - EPIC
## PROJETO = repertorio de TCR e BCR na coorte TCGA-ACC
##
## CRIACAO DO SCRIPT (data) = 13/06/23
## ULTIMA ATUALIZACAO (data) = 13/06/23
## 
## RESPONSAVEL = Jean Resende
## Big Data
################################################################################

library(EPIC)
library(UCSCXenaTools)
library(IOBR)
library(tidyr)

#install.packages("devtools")
#devtools::install_github("GfellerLab/EPIC", build_vignettes=TRUE)

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
################################################################################

load("eset_acc_20230607.RData")
eset_acc_tpm_log2 <- log2(eset_acc+1) # transformacao em log2

res.bref <- EPIC(eset_acc, reference = BRef, withOtherCells = F)
res_bref <- res.bref$cellFractions
save(res_bref, file = "res_bref.RData")

res.tref <- EPIC(eset_acc, reference = TRef, withOtherCells = F)
res_tref <- res.tref$cellFractions
save(res_tref, file = "res_tref.RData")




