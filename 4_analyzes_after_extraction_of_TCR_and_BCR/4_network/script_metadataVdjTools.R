
load("metadata.RData")

sampleNames <- dir("~/Documentos/input_vdjtools_vdjformat/")
head(paste("~/Documentos/input_vdjtools_vdjformat/",sampleNames, sep = ""))

metadataVdjTools <- data.frame(
  file_name = paste("~/Documentos/input_vdjtools_vdjformat/",sampleNames, sep = ""),
  sample_id = gsub("_report_vdjformat.tsv","",sampleNames)
)

idx <- match(metadataVdjTools$sample_id, metadata$sample_id)
metadataVdjTools$steroid <- metadata$steroid[idx]

head(metadataVdjTools)

write.table(metadataVdjTools, file = "metadataVdjTools.txt", row.names = F)
?write.table
