data.z <- mcp[,-1]
data.z <- scale(data.z)
summary(data.z)

#### -- z-score
data <- t(data.z)
colnames(data) <- mcp$ID

pdf(file = "cellBarPlot_mcp_zscore.pdf", width = 8.27, height = 11.69)
barplot(data, 
        col=colors()[c(23,89,12,30,5,18,75,90,120,150)] , 
        border="white", 
        space=0.04, 
        font.axis=2, 
        xlab="group",
        legend=rownames(data))
dev.off()

#### -- score
data <- t(mcp[,-1])
colnames(data) <- mcp$ID

pdf(file = "cellBarPlot_mcp.pdf", width = 20, height = 16)
barplot(data, 
        horiz = F,
        las = 2,
        cex.axis = 0.2,
        col=colors()[c(23,89,12,30,5,18,75,90,120,150)] , 
        border="white", 
        space=0.04, 
        font.axis=2, 
        xlab="group",
        legend=rownames(data))
dev.off()

### -- xcell
load("eset_acc.RData")
eset_acc_tpm_log2 <- log2(eset_acc+1) # transformacao em log2
xcell <- deconvo_tme(eset_acc_tpm_log2, method = "xcell", arrays = FALSE)

data <- t(xcell[,-1])
colnames(data) <- xcell$ID

pdf(file = "cellBarPlot_xcell.pdf", width = 8.27, height = 11.69)
barplot(data, 
        horiz = F,
        las = 2,
        cex.axis = 0.2,
        #col=colors()[c(23,89,12,30,5,18,75,90,120,150)] , 
        border="white", 
        space=0.04, 
        font.axis=2, 
        xlab="group",
        legend=rownames(data))
dev.off()

### -- estimate
load("eset_acc.RData")
eset_acc_tpm_log2 <- log2(eset_acc+1) # transformacao em log2
estimate<-deconvo_tme(eset = eset_acc_tpm_log2, method = "estimate")

### -- timer
load("eset_acc.RData")
eset_acc_tpm_log2 <- log2(eset_acc+1) # transformacao em log2
timer<-deconvo_tme(eset = eset_acc_tpm_log2, method = "timer", group_list = rep("acc",dim(eset_acc_tpm_log2)[2]))

### -- ips
load("eset_acc.RData")
eset_acc_tpm_log2 <- log2(eset_acc+1) # transformacao em log2
ips<-deconvo_tme(eset = eset_acc_tpm_log2, method = "ips", plot= FALSE)




