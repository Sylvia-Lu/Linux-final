data = read.csv("1.raw_data/raw_counts_training.csv",header = T,row.names = 1)
mirna = as.data.frame(data)
group = factor(c(rep("HCC",5),rep("Control",5)))
# Meta data
pdata = read.csv("1.raw_data/pdata_total.csv",header = T,row.names = 1)
pdata = pdata[c(21:25,31:35),]

# Differential gene expression
BiocManager::install("edgeR")
library(edgeR)
dgelist <- DGEList(counts=mirna, group = group)
dgelist <- calcNormFactors(dgelist)
dgelist <- estimateCommonDisp(dgelist)
dgelist <- estimateTagwiseDisp(dgelist)
norm <- cpm(dgelist,log = F)
dge <- exactTest(dgelist)
data1 <- dge$table
keep <- (abs(data1$logFC)>2) & (data1$PValue<0.05)
DEG <- data1[keep,]
write.csv(data1,"2.DGE/1.total.csv")
write.csv(DEG,"2.DGE/2.DEG.csv")
write.csv(norm,"2.DGE/3.norm.csv")