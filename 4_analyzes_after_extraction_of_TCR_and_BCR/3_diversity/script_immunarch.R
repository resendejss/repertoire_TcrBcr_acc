################################################################################

#install.packages("immunarch")

library(immunarch)

immdata <- repLoad("~/Documentos/outputTrust4_report/")
load("metadata.RData")

metadata$Sample <- paste(metadata$sample_id, "_report", sep = "")
metadata <- metadata[,c(16, 1:15)]
metadata$steroid
metadata$Sample[is.na(metadata$steroid)]
metadata <- metadata[!is.na(metadata$steroid),]
metadata$steroid

#metadata <- na.omit(metadata)
immdata[[2]] <- metadata

# estatistica basica
exp_vol <- repExplore(immdata$data, .method = "volume")
exp_len <- repExplore(immdata$data, .method = "len", .col = "aa")
exp_cnt <- repExplore(immdata$data, .method = "count")

p1 <- vis(exp_vol, .by = c("steroid"), .meta = immdata$meta)
p2 <- vis(exp_len, .by = c("steroid"), .meta = immdata$meta, .points = F)
p3 <- vis(exp_cnt, .by = c("steroid"), .meta = immdata$meta)

# sobreposicao de repertorios
ov <- repOverlap(immdata$data)
vis(ov, "circos", annotationTrack="grid")
vis_heatmap2(ov, show_rownames=F, show_colnames=F)

imm_ov1 <- repOverlap(immdata$data, .method = "public", .verbose = F)
vis_heatmap2(imm_ov1, show_rownames=F, show_colnames=F)
vis(imm_ov1)

imm_ov2 <- repOverlap(immdata$data, .method = "jaccard", .verbose = F)
vis_heatmap2(imm_ov2, show_rownames=F, show_colnames=F)

pr.aav <- pubRep(immdata$data, "aa+v", .verbose = F)
write.csv(pr.aav, file="")
