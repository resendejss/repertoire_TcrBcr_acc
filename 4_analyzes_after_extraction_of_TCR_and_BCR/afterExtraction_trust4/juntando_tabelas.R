################################################################################
## post-processing TRUST4
## jean resende
## big data
################################################################################
arquivos <- list.files("output_trust-stats", pattern = ".tsv$")

for (arquivo in arquivos) {
  nome_objeto <- gsub(".tsv", "", arquivo)
  assign(nome_objeto, read.table(paste("output_trust-stats/", arquivo, sep = ""), sep = "\t", header = FALSE))
}

rm(arquivo,arquivos,nome_objeto)


nomes_tabelas <- ls()

# Definir nomes das colunas
nomes_colunas <- c("chain","Abundance","Richness","CPK","Entropy","Clonality")

# Usar o lapply e a função setNames para nomear as colunas de todas as tabelas
lista_tabelas <- lapply(nomes_tabelas, function(nome) setNames(get(nome), nomes_colunas))

# Substituir as tabelas originais pelos novos data frames com nomes de colunas
for (i in 1:length(nomes_tabelas)) {
  assign(nomes_tabelas[i], lista_tabelas[[i]])
}

rm(i, nomes_colunas,nomes_tabelas,lista_tabelas)

lista_tabelas <- mget(ls(pattern = "^1"))

# Junte as tabelas usando rbind
tabela_total <- do.call(rbind, lista_tabelas)

nomes_tabelas <- names(lista_tabelas)
tabela_total$nome_tabela <- rep(nomes_tabelas, sapply(lista_tabelas, nrow))

# Visualize a tabela final
tabela_total
rownames(tabela_total) <- NULL
colnames(tabela_total) <- c("chain","Abundance","Richness","CPK","Entropy",
                            "Clonality","sample_id")

save(tabela_total, file = "tabela_total.RData")

## -- abundance -- 
abundance <- tabela_total[,c("chain","Abundance","sample_id")]
head(abundance)

data <- data.frame(IGH = rep(NA,76),
                   IGHM = rep(NA,76),
                   IGHD = rep(NA,76),
                   IGHG3 = rep(NA,76),
                   IGHG1 = rep(NA,76),
                   IGHA1 = rep(NA,76),
                   IGHG2 = rep(NA,76),
                   IGHG4 = rep(NA,76),
                   IGHE = rep(NA,76),
                   IGHA2 = rep(NA,76),
                   IGK = rep(NA,76),
                   IGL = rep(NA,76),
                   TRA = rep(NA,76),
                   TRB = rep(NA,76),
                   TRG = rep(NA,76),
                   TRD = rep(NA,76))

rownames(data) <- nomes_tabelas

abundance$Abundance[abundance$chain=="IGH"&abundance$sample_id=="130723_UNC9-SN296_0386_BC2E4WACXX_ATCACG_L003_report"]

for (i in nomes_tabelas) {
  data$IGH[rownames(data)==i] <- abundance$Abundance[abundance$chain=="IGH"&abundance$sample_id==i]
  data$IGHM[rownames(data)==i] <- abundance$Abundance[abundance$chain=="IGHM"&abundance$sample_id==i]
  data$IGHD[rownames(data)==i] <- abundance$Abundance[abundance$chain=="IGHD"&abundance$sample_id==i]
  data$IGHG3[rownames(data)==i] <- abundance$Abundance[abundance$chain=="IGHG3"&abundance$sample_id==i]
  data$IGHG1[rownames(data)==i] <- abundance$Abundance[abundance$chain=="IGHG1"&abundance$sample_id==i]
  data$IGHA1[rownames(data)==i] <- abundance$Abundance[abundance$chain=="IGHA1"&abundance$sample_id==i]
  data$IGHG2[rownames(data)==i] <- abundance$Abundance[abundance$chain=="IGHG2"&abundance$sample_id==i]
  data$IGHG4[rownames(data)==i] <- abundance$Abundance[abundance$chain=="IGHG4"&abundance$sample_id==i]
  data$IGHE[rownames(data)==i] <- abundance$Abundance[abundance$chain=="IGHE"&abundance$sample_id==i]
  data$IGHA2[rownames(data)==i] <- abundance$Abundance[abundance$chain=="IGHA2"&abundance$sample_id==i]
  data$IGK[rownames(data)==i] <- abundance$Abundance[abundance$chain=="IGK"&abundance$sample_id==i]
  data$IGL[rownames(data)==i] <- abundance$Abundance[abundance$chain=="IGL"&abundance$sample_id==i]
  data$TRA[rownames(data)==i] <- abundance$Abundance[abundance$chain=="TRA"&abundance$sample_id==i]
  data$TRB[rownames(data)==i] <- abundance$Abundance[abundance$chain=="TRB"&abundance$sample_id==i]
  data$TRG[rownames(data)==i] <- abundance$Abundance[abundance$chain=="TRG"&abundance$sample_id==i]
  data$TRD[rownames(data)==i] <- abundance$Abundance[abundance$chain=="TRD"&abundance$sample_id==i]
}

data.abundance <- data
save(data.abundance, file = "data.abundance.RData")

## -- richness -- 
richness <- tabela_total[,c("chain","Richness","sample_id")]
head(richness)

data <- data.frame(IGH = rep(NA,76),
                   IGHM = rep(NA,76),
                   IGHD = rep(NA,76),
                   IGHG3 = rep(NA,76),
                   IGHG1 = rep(NA,76),
                   IGHA1 = rep(NA,76),
                   IGHG2 = rep(NA,76),
                   IGHG4 = rep(NA,76),
                   IGHE = rep(NA,76),
                   IGHA2 = rep(NA,76),
                   IGK = rep(NA,76),
                   IGL = rep(NA,76),
                   TRA = rep(NA,76),
                   TRB = rep(NA,76),
                   TRG = rep(NA,76),
                   TRD = rep(NA,76))

rownames(data) <- nomes_tabelas

richness$Richness[richness$chain=="IGH"&richness$sample_id=="130723_UNC9-SN296_0386_BC2E4WACXX_ATCACG_L003_report"]

for (i in nomes_tabelas) {
  data$IGH[rownames(data)==i] <- richness$Richness[richness$chain=="IGH"&richness$sample_id==i]
  data$IGHM[rownames(data)==i] <- richness$Richness[richness$chain=="IGHM"&richness$sample_id==i]
  data$IGHD[rownames(data)==i] <- richness$Richness[richness$chain=="IGHD"&richness$sample_id==i]
  data$IGHG3[rownames(data)==i] <- richness$Richness[richness$chain=="IGHG3"&richness$sample_id==i]
  data$IGHG1[rownames(data)==i] <- richness$Richness[richness$chain=="IGHG1"&richness$sample_id==i]
  data$IGHA1[rownames(data)==i] <- richness$Richness[richness$chain=="IGHA1"&richness$sample_id==i]
  data$IGHG2[rownames(data)==i] <- richness$Richness[richness$chain=="IGHG2"&richness$sample_id==i]
  data$IGHG4[rownames(data)==i] <- richness$Richness[richness$chain=="IGHG4"&richness$sample_id==i]
  data$IGHE[rownames(data)==i] <- richness$Richness[richness$chain=="IGHE"&richness$sample_id==i]
  data$IGHA2[rownames(data)==i] <- richness$Richness[richness$chain=="IGHA2"&richness$sample_id==i]
  data$IGK[rownames(data)==i] <- richness$Richness[richness$chain=="IGK"&richness$sample_id==i]
  data$IGL[rownames(data)==i] <- richness$Richness[richness$chain=="IGL"&richness$sample_id==i]
  data$TRA[rownames(data)==i] <- richness$Richness[richness$chain=="TRA"&richness$sample_id==i]
  data$TRB[rownames(data)==i] <- richness$Richness[richness$chain=="TRB"&richness$sample_id==i]
  data$TRG[rownames(data)==i] <- richness$Richness[richness$chain=="TRG"&richness$sample_id==i]
  data$TRD[rownames(data)==i] <- richness$Richness[richness$chain=="TRD"&richness$sample_id==i]
}

data.richness <- data
save(data.richness, file = "data.richness.RData")

## -- CPK -- 
cpk <- tabela_total[,c("chain","CPK","sample_id")]
head(cpk)
rownames(cpk) <- NULL

data <- data.frame(IGH = rep(NA,76),
                   IGHM = rep(NA,76),
                   IGHD = rep(NA,76),
                   IGHG3 = rep(NA,76),
                   IGHG1 = rep(NA,76),
                   IGHA1 = rep(NA,76),
                   IGHG2 = rep(NA,76),
                   IGHG4 = rep(NA,76),
                   IGHE = rep(NA,76),
                   IGHA2 = rep(NA,76),
                   IGK = rep(NA,76),
                   IGL = rep(NA,76),
                   TRA = rep(NA,76),
                   TRB = rep(NA,76),
                   TRG = rep(NA,76),
                   TRD = rep(NA,76))

rownames(data) <- nomes_tabelas

cpk$CPK[cpk$chain=="IGH"&cpk$sample_id=="130723_UNC9-SN296_0386_BC2E4WACXX_ATCACG_L003_report"]

for (i in nomes_tabelas) {
  data$IGH[rownames(data)==i] <- cpk$CPK[cpk$chain=="IGH"&cpk$sample_id==i]
  data$IGHM[rownames(data)==i] <- cpk$CPK[cpk$chain=="IGHM"&cpk$sample_id==i]
  data$IGHD[rownames(data)==i] <- cpk$CPK[cpk$chain=="IGHD"&cpk$sample_id==i]
  data$IGHG3[rownames(data)==i] <- cpk$CPK[cpk$chain=="IGHG3"&cpk$sample_id==i]
  data$IGHG1[rownames(data)==i] <- cpk$CPK[cpk$chain=="IGHG1"&cpk$sample_id==i]
  data$IGHA1[rownames(data)==i] <- cpk$CPK[cpk$chain=="IGHA1"&cpk$sample_id==i]
  data$IGHG2[rownames(data)==i] <- cpk$CPK[cpk$chain=="IGHG2"&cpk$sample_id==i]
  data$IGHG4[rownames(data)==i] <- cpk$CPK[cpk$chain=="IGHG4"&cpk$sample_id==i]
  data$IGHE[rownames(data)==i] <- cpk$CPK[cpk$chain=="IGHE"&cpk$sample_id==i]
  data$IGHA2[rownames(data)==i] <- cpk$CPK[cpk$chain=="IGHA2"&cpk$sample_id==i]
  data$IGK[rownames(data)==i] <- cpk$CPK[cpk$chain=="IGK"&cpk$sample_id==i]
  data$IGL[rownames(data)==i] <- cpk$CPK[cpk$chain=="IGL"&cpk$sample_id==i]
  data$TRA[rownames(data)==i] <- cpk$CPK[cpk$chain=="TRA"&cpk$sample_id==i]
  data$TRB[rownames(data)==i] <- cpk$CPK[cpk$chain=="TRB"&cpk$sample_id==i]
  data$TRG[rownames(data)==i] <- cpk$CPK[cpk$chain=="TRG"&cpk$sample_id==i]
  data$TRD[rownames(data)==i] <- cpk$CPK[cpk$chain=="TRD"&cpk$sample_id==i]
}

data.cpk <- data
save(data.cpk, file = "data.cpk.RData")

## -- entropy -- 
entropy <- tabela_total[,c("chain","Entropy","sample_id")]
head(entropy)
rownames(entropy) <- NULL

data <- data.frame(IGH = rep(NA,76),
                   IGHM = rep(NA,76),
                   IGHD = rep(NA,76),
                   IGHG3 = rep(NA,76),
                   IGHG1 = rep(NA,76),
                   IGHA1 = rep(NA,76),
                   IGHG2 = rep(NA,76),
                   IGHG4 = rep(NA,76),
                   IGHE = rep(NA,76),
                   IGHA2 = rep(NA,76),
                   IGK = rep(NA,76),
                   IGL = rep(NA,76),
                   TRA = rep(NA,76),
                   TRB = rep(NA,76),
                   TRG = rep(NA,76),
                   TRD = rep(NA,76))

rownames(data) <- nomes_tabelas

entropy$Entropy[entropy$chain=="IGH"&entropy$sample_id=="130723_UNC9-SN296_0386_BC2E4WACXX_ATCACG_L003_report"]

for (i in nomes_tabelas) {
  data$IGH[rownames(data)==i] <- entropy$Entropy[entropy$chain=="IGH"&entropy$sample_id==i]
  data$IGHM[rownames(data)==i] <- entropy$Entropy[entropy$chain=="IGHM"&entropy$sample_id==i]
  data$IGHD[rownames(data)==i] <- entropy$Entropy[entropy$chain=="IGHD"&entropy$sample_id==i]
  data$IGHG3[rownames(data)==i] <- entropy$Entropy[entropy$chain=="IGHG3"&entropy$sample_id==i]
  data$IGHG1[rownames(data)==i] <- entropy$Entropy[entropy$chain=="IGHG1"&entropy$sample_id==i]
  data$IGHA1[rownames(data)==i] <- entropy$Entropy[entropy$chain=="IGHA1"&entropy$sample_id==i]
  data$IGHG2[rownames(data)==i] <- entropy$Entropy[entropy$chain=="IGHG2"&entropy$sample_id==i]
  data$IGHG4[rownames(data)==i] <- entropy$Entropy[entropy$chain=="IGHG4"&entropy$sample_id==i]
  data$IGHE[rownames(data)==i] <- entropy$Entropy[entropy$chain=="IGHE"&entropy$sample_id==i]
  data$IGHA2[rownames(data)==i] <- entropy$Entropy[entropy$chain=="IGHA2"&entropy$sample_id==i]
  data$IGK[rownames(data)==i] <- entropy$Entropy[entropy$chain=="IGK"&entropy$sample_id==i]
  data$IGL[rownames(data)==i] <- entropy$Entropy[entropy$chain=="IGL"&entropy$sample_id==i]
  data$TRA[rownames(data)==i] <- entropy$Entropy[entropy$chain=="TRA"&entropy$sample_id==i]
  data$TRB[rownames(data)==i] <- entropy$Entropy[entropy$chain=="TRB"&entropy$sample_id==i]
  data$TRG[rownames(data)==i] <- entropy$Entropy[entropy$chain=="TRG"&entropy$sample_id==i]
  data$TRD[rownames(data)==i] <- entropy$Entropy[entropy$chain=="TRD"&entropy$sample_id==i]
}

data.entropy <- data
save(data.entropy, file = "data.entropy.RData")

## -- clonality -- 
clonality <- tabela_total[,c("chain","Clonality","sample_id")]
head(clonality)
rownames(Clonality) <- NULL

data <- data.frame(IGH = rep(NA,76),
                   IGHM = rep(NA,76),
                   IGHD = rep(NA,76),
                   IGHG3 = rep(NA,76),
                   IGHG1 = rep(NA,76),
                   IGHA1 = rep(NA,76),
                   IGHG2 = rep(NA,76),
                   IGHG4 = rep(NA,76),
                   IGHE = rep(NA,76),
                   IGHA2 = rep(NA,76),
                   IGK = rep(NA,76),
                   IGL = rep(NA,76),
                   TRA = rep(NA,76),
                   TRB = rep(NA,76),
                   TRG = rep(NA,76),
                   TRD = rep(NA,76))

rownames(data) <- nomes_tabelas

clonality$Clonality[clonality$chain=="IGH"&clonality$sample_id=="130723_UNC9-SN296_0386_BC2E4WACXX_ATCACG_L003_report"]

for (i in nomes_tabelas) {
  data$IGH[rownames(data)==i] <- clonality$Clonality[clonality$chain=="IGH"&clonality$sample_id==i]
  data$IGHM[rownames(data)==i] <- clonality$Clonality[clonality$chain=="IGHM"&clonality$sample_id==i]
  data$IGHD[rownames(data)==i] <- clonality$Clonality[clonality$chain=="IGHD"&clonality$sample_id==i]
  data$IGHG3[rownames(data)==i] <- clonality$Clonality[clonality$chain=="IGHG3"&clonality$sample_id==i]
  data$IGHG1[rownames(data)==i] <- clonality$Clonality[clonality$chain=="IGHG1"&clonality$sample_id==i]
  data$IGHA1[rownames(data)==i] <- clonality$Clonality[clonality$chain=="IGHA1"&clonality$sample_id==i]
  data$IGHG2[rownames(data)==i] <- clonality$Clonality[clonality$chain=="IGHG2"&clonality$sample_id==i]
  data$IGHG4[rownames(data)==i] <- clonality$Clonality[clonality$chain=="IGHG4"&clonality$sample_id==i]
  data$IGHE[rownames(data)==i] <- clonality$Clonality[clonality$chain=="IGHE"&clonality$sample_id==i]
  data$IGHA2[rownames(data)==i] <- clonality$Clonality[clonality$chain=="IGHA2"&clonality$sample_id==i]
  data$IGK[rownames(data)==i] <- clonality$Clonality[clonality$chain=="IGK"&clonality$sample_id==i]
  data$IGL[rownames(data)==i] <- clonality$Clonality[clonality$chain=="IGL"&clonality$sample_id==i]
  data$TRA[rownames(data)==i] <- clonality$Clonality[clonality$chain=="TRA"&clonality$sample_id==i]
  data$TRB[rownames(data)==i] <- clonality$Clonality[clonality$chain=="TRB"&clonality$sample_id==i]
  data$TRG[rownames(data)==i] <- clonality$Clonality[clonality$chain=="TRG"&clonality$sample_id==i]
  data$TRD[rownames(data)==i] <- clonality$Clonality[clonality$chain=="TRD"&clonality$sample_id==i]
}

data.clonality <- data
save(data.clonality, file = "data.clonality.RData")
