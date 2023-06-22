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

vis(exp_vol, .by = c("steroid"), .meta = immdata$meta)
vis(exp_vol, .by = c("gender"), .meta = immdata$meta)
vis(exp_vol, .by = c("stage"), .meta = immdata$meta)
vis(exp_vol, .by = c("Immune.Subtype"), .meta = immdata$meta)
vis(exp_vol, .by = c("cortisol.excess"), .meta = immdata$meta)

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
write.csv(pr.aav, file="pr.avv.csv")


clon_homeo <- repClonality(immdata$data, "clonal.prop")
imm_rare <- repClonality(immdata$data, .method = "rare")
vis(imm_rare)
vis(imm_rare, .by = "steroid", .meta = immdata$meta)

imm_top <- repClonality(immdata$data, .method = "top", .head = c(10, 100, 1000, 3000))
vis(imm_top, .by = "steroid", .meta = immdata$meta)

imm_hom <- repClonality(immdata$data, .method = "homeo")
vis(imm_hom)
vis(imm_hom, .by = "steroid", .meta = immdata$meta)

