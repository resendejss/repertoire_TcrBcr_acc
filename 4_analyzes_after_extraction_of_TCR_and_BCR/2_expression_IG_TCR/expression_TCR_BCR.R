################################################################################
# Name: expression_TCR_BCR                                                     #
#                                                                              #
# Author: Jean Resende (jean.s.s.resende@gmail.com)                            #
# Creation date: 2023/06/15                                                    #
# Last update date: 2023/06/17                                                 #
################################################################################

load("v.counts.RData")
load("metadata.RData")

data <- v.counts
idx <- match(metadata$sample_id, rownames(data))
data <- data[idx,] 

identical(metadata$sample_id,rownames(data))

boxplot(log2(data+1))

length(rowSums(data))

data$IGK <-  data$IGK / metadata$reads
data$IGL <-  data$IGL / metadata$reads
data$IGH <-  data$IGH / metadata$reads
data$TRB <-  data$TRB / metadata$reads
data$TRA <-  data$TRA / metadata$reads
data$TRG <-  data$TRG / metadata$reads
data$TRD <-  data$TRD / metadata$reads

save(data, file = "v.normalizations.RData")
rm(list = ls())
################################################################################

load("d.counts.RData")
load("metadata.RData")

data <- d.counts
idx <- match(metadata$sample_id, rownames(data))
data <- data[idx,] 

identical(metadata$sample_id,rownames(data))

boxplot(log2(data+1))

data$IGK <-  data$IGK / metadata$reads
data$IGL <-  data$IGL / metadata$reads
data$IGH <-  data$IGH / metadata$reads
data$TRB <-  data$TRB / metadata$reads
data$TRA <-  data$TRA / metadata$reads
data$TRG <-  data$TRG / metadata$reads
data$TRD <-  data$TRD / metadata$reads

save(data, file = "d.normalizations.RData")
rm(list = ls())
################################################################################

load("j.counts.RData")
load("metadata.RData")

data <- j.counts
idx <- match(metadata$sample_id, rownames(data))
data <- data[idx,] 

identical(metadata$sample_id,rownames(data))

boxplot(log2(data+1))

data$IGK <-  data$IGK / metadata$reads
data$IGL <-  data$IGL / metadata$reads
data$IGH <-  data$IGH / metadata$reads
data$TRB <-  data$TRB / metadata$reads
data$TRA <-  data$TRA / metadata$reads
data$TRG <-  data$TRG / metadata$reads
data$TRD <-  data$TRD / metadata$reads

save(data, file = "j.normalizations.RData")
rm(list = ls())
################################################################################
