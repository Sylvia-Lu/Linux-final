data = read.csv("2.DGE/3.norm.csv",header = T,row.names = 1)
mirna = as.data.frame(data)
group = factor(c(rep("HCC",5),rep("Control",5)))
# Meta data
pdata = read.csv("1.raw_data/pdata_total.csv",header = T,row.names = 1)
pdata = pdata[c(21:25,31:35),]

# WGCNA
library(WGCNA)
exp = mirna
exp = apply(exp,2,as.numeric)
exp <- t(exp)
gsg <- goodSamplesGenes(exp)
gsg$allOK

if(!gsg$allOK){
  if(sum(!gsg$goodGenes)>0){
    printFlush(paste("Removing genes:",paste(names(exp)[!gsg$goodGenes],collapse=",")))
  }
  if (sum(!gsg$goodSamples)>0){
    printFlush(paste("Removing samples:",paste(rownames(exp)[!gsg$goodSamples],collapse=","))) }
  exp <- exp[gsg$goodSamples,gsg$goodGenes]
}

#soft threshold
powers  =  c(c(1:10), seq(from = 12, to = 30, by = 2))
sft = pickSoftThreshold(exp,networkType = "signed",corFnc = "cor",verbose = 5,powerVector = powers)
pdf(file="3.WGCNA/01.soft_threshold.pdf",width=10,height=6) 
par(mfrow=c(1,2))
cex1 <- 0.9
plot(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft threshold (power)",
     ylab="Scale free topology model fit, signed R^2", type="n",
     main=paste("Scale independence")) 
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red")
abline(h=0.8,col="red")
plot(sft$fitIndices[,1],sft$fitIndices[,5],
     xlab="Soft threshold (power)",ylab="Mean connectivity", type="n",main=paste("Mean connectivity")) 
text(sft$fitIndices[,1],sft$fitIndices[,5],labels=powers,cex=cex1,col="red") 
dev.off()

#adjacency matrix + TOM
wgcna_parameters = list(powers = 18)
wgcna_parameters$minModSize = 50
wgcna_parameters$minHeight = 0.15
wgcna_parameters$bsize = 10000 
wgcna_parameters$ds = 4 
wgcna_parameters$networkType = "signed"  
wgcna_parameters$corFnc = "pearson"
wgcna_parameters$pamStage = FALSE
net <- blockwiseModules(exp, maxBlockSize=wgcna_parameters$bsize, 
                        networkType=wgcna_parameters$networkType, corType = wgcna_parameters$corFnc ,  
                        power = wgcna_parameters$powers, mergeCutHeight= wgcna_parameters$minHeight, 
                        nThreads=8, saveTOMFileBase="3.WGCNA/02.WGCNA_TOM", 
                        saveTOMs=TRUE, minModuleSize= wgcna_parameters$minModSize, pamStage=wgcna_parameters$pamStage, 
                        reassignThreshold=1e-6, verbose = 3, deepSplit=wgcna_parameters$ds)
load(file = "3.WGCNA/02.WGCNA_TOM-block.1.Rdata")
TOM = as.matrix(TOM)
geneTree = hclust(1-as.dist(TOM), method="average") 
tree = cutreeHybrid(dendro = geneTree, minClusterSize= wgcna_parameters$minModSize, pamStage=wgcna_parameters$pamStage, cutHeight = 0.999, 
                    deepSplit=wgcna_parameters$ds, distM=as.matrix(1-as.dist(TOM)))
merged = mergeCloseModules(exprData= exp, colors = tree$labels, cutHeight=wgcna_parameters$minHeight)
colors = labels2colors(merged$colors)
table(colors)
length(table(colors))
pdf("3.WGCNA/03.final_modules.pdf",width = 20, height = 12)
plotDendroAndColors(geneTree,colors,groupLabels = "modules",cex.colorLabels = 0.5,addGuide=T,dendroLabels=F)
dev.off()
save(file = "3.WGCNA/04.module_results.Rdata",wgcna_parameters,tree,merged,colors,sft)

#Module eigengenes
MEs0 = moduleEigengenes(expr = exp, colors, softPower = wgcna_parameters$powers)$eigengenes
MEs <- orderMEs(MEs0)
text <- cbind(rownames(MEs),MEs) 
colnames(text)[1] <- "samples"
write.csv(text,file="3.WGCNA/05.module_eigengenes.csv",quote=F,row.names=F)

names(MEs) <- substring(names(MEs),3)
MEDiss <- 1-cor(MEs)
METree <- hclust(as.dist(MEDiss),method="average") 
pdf(file="3.WGCNA/06.modules_cluster_tree.pdf",width=7,height=5) 
plot(METree,main="Clustering of module eigengenes",xlab="",sub="") 
dev.off()
save(file = "3.WGCNA/07.module_eigengenes.Rdata",MEs,MEDiss,METree)

#pre-processing of datMeta
all_colors = unique(colors)
all_colors = all_colors[-6]
allTraits<- as.data.frame(pdata$source_name_ch1)
names(allTraits) = 'group'
rownames(allTraits) = rownames(pdata)
allTraits$group=as.character(allTraits$group)
allTraits$group[allTraits$group=="Normal plasma"] <- 0
allTraits$group[!allTraits$group=="0"] <- 1
sapply(1:dim(allTraits)[2],function(x){
  allTraits[,x]<<-  as.numeric(allTraits[,x])
})

#module_traits_relationships
moduleTraitCor <- cor(MEs,allTraits,use = "p",method = "pearson")
nSamples <- nrow(exp[,net$goodGenes])
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,nSamples)
moduleTraitPvalue_fdr <- p.adjust(moduleTraitPvalue, method = 'fdr')
textMatrix <- paste(signif(moduleTraitCor,2),"\n(",signif(moduleTraitPvalue_fdr,1),")",sep="") 
dim(textMatrix) <- dim(moduleTraitCor)
rowLabels <- paste("ME",names(MEs),sep="") 
pdf(file="3.WGCNA/08.modules_traits_relationships.pdf",10,8) 
labeledHeatmap(Matrix=moduleTraitCor,textMatrix=textMatrix, xLabels=colnames(allTraits), yLabels=rowLabels, ySymbols=names(MEs),setStdMargins=FALSE, xLabelsAngle=90,zlim=c(-1,1))
dev.off()

text <- paste("cor=",round(moduleTraitCor,4), ";p-value=",round(moduleTraitPvalue_fdr,8),sep="")
dim(text) <- dim(moduleTraitCor)
rownames(text) <- rownames(moduleTraitCor)
colnames(text) <- colnames(moduleTraitCor)
text <- cbind(rownames(text),text)
colnames(text)[1] <- "modules" 
write.csv(text,file="3.WGCNA/08.modules_traits_relationships.csv",
          quote=F,row.names=F)
save(file = "./3.WGCNA/09.module_traits_relationships.Rdata",exp,allTraits,moduleTraitCor,moduleTraitPvalue)

keep <- intersect(names(moduleTraitCor[which(abs(rowSums(moduleTraitCor)) >= 0.8),]),names(moduleTraitPvalue[which(rowSums(moduleTraitPvalue) <= 0.05),]))
#Gene module membership
modNames <- names(MEs)
geneModuleMembership <- cor(exp[,net$goodGenes],MEs,use="p",method = "pearson") 
nSamples <- nrow(exp[,net$goodGenes])
MMPvalue <- corPvalueStudent(geneModuleMembership,nSamples) 
colnames(geneModuleMembership) <- paste("MM",modNames,sep="") 
colnames(MMPvalue) <- paste("p.MM",modNames,sep="")
text <- paste("cor=",round(geneModuleMembership,4), ";p-value=",round(MMPvalue,4),sep="")
dim(text) <- dim(geneModuleMembership)
row.names(geneModuleMembership) = rownames(norm)
row.names(text) = rownames(norm)
colnames(text) = modNames
write.csv(text,file="3.WGCNA/10.genes_module_membership.csv")

#Gene trait significance
modNames <- names(MEs)
status <- as.data.frame(allTraits$group)
names(status) <- "status"
geneTraitSignificance <- as.data.frame(cor(exp[,net$goodGenes],status,use="p")) 
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance),nSamples)) 
names(geneTraitSignificance) <- paste("GS.",names(status),sep="")
names(GSPvalue) <- paste("p.GS.",names(status),sep="")
gene_colors <- as.data.frame(colors)
geneTraitSignificance_colors <- cbind(geneTraitSignificance,gene_colors)
row.names(geneTraitSignificance) <- rownames(norm)
row.names(geneTraitSignificance_colors) <- rownames(norm)

#GS vs MM
dir.create("3.WGCNA/11.MM_vs_status")
for(i in 1:(ncol(MEs)-1)) {
  which.module <- all_colors[i]
  column <- match(which.module,modNames)
  moduleGenes <- colors[net$blockGenes[[1]]]==which.module
  dir <- "3.WGCNA/11.MM_vs_status/" 
  pdf(file=paste(dir,"11.",which.module,"_MM_vs_status.pdf",sep=""))
  par(mfrow = c(1,1))
  verboseScatterplot(abs(geneModuleMembership[moduleGenes,column]),   
                     abs(geneTraitSignificance[moduleGenes,1]),
                     xlab=paste("Module membership (MM) in",which.module,"module"), ylab="Gene significance for status",
                     main=paste("Module membership vs. gene significance"), col=which.module)
  dev.off() 
} 
text <- cbind(geneTraitSignificance,GSPvalue)
text <- cbind(rownames(text),text)
colnames(text)[1] <- "genes"
write.csv(text,file="3.WGCNA/12.genes_trait_significance.csv",
          quote=F, row.names=F)
save(file = "3.WGCNA/13.GS_vs_MM.Rdata",geneModuleMembership,MMPvalue,geneTraitSignificance,GSPvalue,geneTraitSignificance_colors)

dir.create("3.WGCNA/14.module_hub_genes")
setwd("3.WGCNA/14.module_hub_genes")
#Output hub genes with high GS and MM
for(i in 1:(ncol(MEs)-1)) {
  module <- all_colors[i]
  inModule <- colors[net$blockGenes[[1]]]==module 
  genename <- rownames(norm[net$goodGenes,]) 
  modGenes <- genename[inModule]
  gene_trait_significance <- geneTraitSignificance_colors[modGenes,]
  high_GS <- rownames(gene_trait_significance[abs(gene_trait_significance$GS.status)>0.7,])
  rownames_MM <- paste("MM",all_colors[i],sep="")
  high_MM <- rownames(geneModuleMembership[abs(geneModuleMembership[,rownames_MM])>0.7,])                             
  module_hub_genes <- intersect(high_GS,high_MM)
  deg <- rownames(DEG)
  common <- intersect(deg,module_hub_genes)
  write.table(module_hub_genes,paste(module,"_high_GS&MM_hubgene.txt",sep=""),sep="\t",quote=F,row.names=F)
  write.table(common,paste(module,"hubgene.txt",sep=""),sep="\t",quote=F,row.names=F)
}
