########################################################################################################
#Get all genes into WGCNA but remove the genes with low basemean values
#navigate into the correct directory
setwd("/usr4/bi594/ebrenner/ondemand/GitHub-final")
#source("http://bioconductor.org/biocLite.R") #To download DESeq package (you can comment these lines out, they only need to be run once ever)
#biocLite("WGCNA")
#biocLite("flashClust")
#If using R version or greater you need to use BiocManager instead:
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("WGCNA")
#^This should download both WGCNA and flashClust
library(WGCNA) #Version 1.70-3
library(impute) #Version 1.64.0
library(flashClust) #Version 1.01-2
library(DESeq2) #Version  1.30.1

countData <- read.table("B_psygmophilum_counts.txt")
head(countData)
length(countData[,1])
#28265
#changes names as we did with deseq for ease
names(countData)=c( "Control_1.1", "Control_1", "Control_2.1", "Control_2", "Control_3.1", "Control_3", "Control_4.1", "Control_4", "Cool_1", "Cool_2", "Cool_3", "Cool_4","Heat_1.1", "Heat_1", "Heat_2.1", "Heat_2", "Heat_3.1", "Heat_3", "Heat_4.1", "Heat_4")
##row.names(countData)=sub("", "isogroup", rownames(countData))
head(countData)
#create call data function, not really needed since blind is true
treat=c( "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Cool", "Cool", "Cool", "Cool","Heat", "Heat", "Heat", "Heat", "Heat", "Heat", "Heat", "Heat")
g=data.frame(treat)
g
colData<- g

dds<-DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~ treat) 
dds<-DESeq(dds)

res<- results(dds)
#filter for contigs with average(baseMean) >3. gets rid of trash counts that would just be annoying
res3<-res[res$baseMean>3, ]
dim(res) #28265
dim(res3) #15962

#could use vst here instead for large dataset. get rlog data (better transformation when size factors vary across samples)
rld <- rlogTransformation(dds, blind=TRUE, fitType="local")
head(rld)
rld_wg=(assay(rld))
head(rld_wg)
nrow(rld_wg)
#28265
#filter like before
rldFiltered=(assay(rld))[(rownames((assay(rld))) %in% rownames(res3)),]
nrow(rldFiltered)
#15962
write.csv( rldFiltered,file="Bpsy_wgcna_allgenes.csv",quote=F,row.names=T)
#now we have our filtered data to take into WGCNA



####First part of tutorial:Data input and cleaning

options(stringsAsFactors=FALSE)
allowWGCNAThreads()

dat=read.csv("Bpsy_wgcna_allgenes.csv")
head(dat) 
rownames(dat)<-dat$X
head(dat)
dat$X=NULL
head(dat)
names(dat)
nrow(dat)
#15962
datExpr0 = as.data.frame(t(dat))

gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK #if TRUE, no outlier genes, if false run the script below



### Outlier detection incorporated into trait measures. check names are correct 
traitData= read.csv("Bpsy_Traits_Num.csv", row.names=1)
dim(traitData)
head(traitData)
names(traitData)

# Form a data frame analogous to expression data that will hold the clinical traits.
dim(datExpr0)
rownames(datExpr0)
rownames(traitData)=rownames(datExpr0)
traitData$Sample= NULL 
# datTraits=allTraits
datTraits=traitData

table(rownames(datTraits)==rownames(datExpr0)) #should return TRUE if datasets align correctly, otherwise your names are out of order
head(datTraits)
head(datExpr0)

#sample dendrogram and trait heat map showing outliers. how linked genes are
A=adjacency(t(datExpr0),type="signed") #default is unsigned, we want signed
# this calculates the whole network connectivity we choose signed because we care about direction of gene expression
k=as.numeric(apply(A,2,sum))-1
# standardized connectivity
Z.k=scale(k)
thresholdZ.k=-2.5 # often -2.5
outlierColor=ifelse(Z.k<thresholdZ.k,"red","black")
sampleTree = flashClust(as.dist(1-A), method = "average")
# Convert traits to a color representation where red indicates high values
traitColors=data.frame(numbers2colors(datTraits,signed=FALSE))
dimnames(traitColors)[[2]]=paste(names(datTraits))
datColors=data.frame(outlierC=outlierColor,traitColors)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree,groupLabels=names(datColors), colors=datColors,main="Sample dendrogram and trait heatmap")


#everything we did above to make plot was pretty standard
# Remove outlying samples from expression and trait data
# remove.samples= Z.k<thresholdZ.k | is.na(Z.k)
# datExpr=datExpr0[!remove.samples,]
# datTraits=datTraits[!remove.samples,]

#we save this because it can take a lot to actually make this with large dataset
save(datExpr0, datTraits, file="Bpsy_Samples_Traits_ALL.RData")

################Moving on!  Network construction and module detection - this section can take a lot of time you might consider running it on a cluster for a larger dataset
library(WGCNA)
library(flashClust)
options(stringsAsFactors = FALSE)
#enableWGCNAThreads use this in base R
allowWGCNAThreads() 
lnames = load(file="Bpsy_Samples_Traits_ALL.RData")

#Figure out proper SFT
# Choose a set of soft-thresholding powers
powers = c(seq(1,14,by=2), seq(15,30, by=0.5)); #may need to adjust these power values to hone in on proper sft value
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr0, powerVector = powers, networkType="signed", verbose = 2) #want smallest value, closest to 0.9 (but still under)

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
#run lines 153 to 164 together
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",ylim=c(0,1),
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

softPower=9 #smallest value to plateau at ~0.85
adjacency=adjacency(datExpr0, power=softPower,type="signed") #must change method type here too!!
#translate the adjacency into topological overlap matrix and calculate the corresponding dissimilarity:
TOM= TOMsimilarity(adjacency,TOMType = "signed")
dissTOM= 1-TOM

library(flashClust)
geneTree= flashClust(as.dist(dissTOM), method="average")
sizeGrWindow(10,6)
# pdf(file="dendrogram_thresh16.5_signed_1868.pdf", width=20, height=20)
plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE,hang=0.04)
# dev.off()
#each leaf corresponds to a gene, branches grouping together densely are interconnected, highly co-expressed genes

minModuleSize=90 #we only want large modules
dynamicMods= cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=2, pamRespectsDendro= FALSE, minClusterSize= minModuleSize)
table(dynamicMods)

#dynamicMods
#1    2    3    4    5    6    7    8    9   10   11   12   13   14 
#3940 1974 1781 1604 1139 1047  955  920  784  744  458  218  207  191 

dynamicColors= labels2colors(dynamicMods)
#plot dendrogram and colors underneath, pretty sweet
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang= 0.05, main= "Gene dendrogram and module colors")

#Merg modules whose expression profiles are very similar
#calculate eigengenes
MEList= moduleEigengenes(datExpr0, colors= dynamicColors,softPower = 9)
MEs= MEList$eigengenes
#Calculate dissimilarity of module eigenegenes
MEDiss= 1-cor(MEs)
#Cluster module eigengenes
METree= flashClust(as.dist(MEDiss), method= "average")

save(dynamicMods, MEList, MEs, MEDiss, METree, file= "Network_bpsy_nomerge.RData")

lnames = load(file = "Network_bpsy_nomerge.RData")
#plot
sizeGrWindow(7,6)
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")

MEDissThres= 0.4
abline(h=MEDissThres, col="red")

merge= mergeCloseModules(datExpr0, dynamicColors, cutHeight= MEDissThres, verbose =3)

mergedColors= merge$colors
mergedMEs= merge$newMEs

pdf(file="MergeNetwork.pdf", width=20, height=20)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)
dev.off()

moduleColors= mergedColors
colorOrder= c("turquoise", standardColors(50))
moduleLabels= match(moduleColors, colorOrder)-1
MEs=mergedMEs

#save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file= "Network_signed_0.6.RData")

###############Relating modules to traits and finding important genes
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
lnames = load(file = "Bpsy_Samples_Traits_ALL.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "Network_signed_0.6.RData");
lnames = load(file = "Network_bpsy_nomerge.RData");
lnames

nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)
table(moduleColors)

#moduleColors
#brown        cyan greenyellow     magenta   turquoise 
#4731         191        3176        2659        5205 

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

#represent module trait correlations as a heatmap
quartz()
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
#plot looks like we need to cut higher to merge more groups.
###based on a lot of red occuring tight together. create better modules
#####changed to MEDissThres= 0.4
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

#Gene relationship to trait and important modules:
#############################################################################
#Magenta module
# Define variable weight containing the weight column of datTrait - leave weight as variable, but change names in first 2 commands
weight = as.data.frame(datTraits$Cool); #change Lipidrobust to your trait name
names(weight) = "Cool" #change based on the module you are running below
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr0, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="")


#Gene-trait significance correlation plots. change module to look at different ones
# par(mfrow=c(2,3))
module = "magenta"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("ModMem in", module, "module"),
                   ylab = "Gene Sig for Cool",
                   main = paste("MM vs. GS\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

#Making VSD files by module for GO plot functions
vs=t(datExpr0)
cands=names(datExpr0[moduleColors=="magenta"]) #black  blue brown green  grey  pink   red 
#looking at genes in this module and subsetting them
c.vsd=vs[rownames(vs) %in% cands,]
head(c.vsd)
nrow(c.vsd) #should correspond to module size
table(moduleColors)
#moduleColors
#brown        cyan greenyellow     magenta   turquoise 
#4731         191        3176        2659        5205 
head(c.vsd)
#creates csv file with subsetted data
write.csv(c.vsd,"rlog_MMmagenta.csv",quote=F)

##############################heatmap of module expression with bar plot of eigengene, no resorting of samples...
#names(dis)
sizeGrWindow(8,7);
which.module="magenta" #pick module of interest
ME=MEs[, paste("ME",which.module, sep="")]
genes=datExpr0[,moduleColors==which.module ] #replace where says subgene below to plot all rather than just subset

#quartz()
# par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
par(mfrow=c(2,1), mar=c(0.3, 5.5, 5, 2))
plotMat(t(scale(genes) ),nrgcols=30,rlabels=F, clabels=rownames(genes), rcols=which.module)
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="sample")
#eigengene is overall expression
#this is a cool plot where you can see that genes in this module are upregulated in the cool treatment


#Gene relationship to trait and important modules: Gene Significance and Module membership
allkME =as.data.frame(signedKME(t(dat), MEs))
head(allkME)
vsd=read.csv(file="rlog_MMmagenta.csv", row.names=1)
head(vsd)

library(pheatmap)

############################################
#top 100 genes
whichModule="magenta"
top=100

datME=MEs
vsd <- read.csv("Bpsy_wgcna_allgenes.csv", row.names=1)
head(vsd)
datExpr=t(vsd)
modcol=paste("kME",whichModule,sep="")
head(vsd)
sorted=vsd[order(allkME[,modcol],decreasing=T),]
hubs=sorted[1:top,]
# attaching gene names
summary(hubs)



contrasting = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
#quartz()
pheatmap(hubs,scale="row",col=contrasting,border_color=NA, main=paste(whichModule,"top",top,"kME",sep=""))

###fisher for GO
# two ways to do GO: Fischers- every you in the module or not or KME: how well gene belongs, more continuous
##########fisher of module vs whole dataset
library(WGCNA)
vsd <- read.csv("Bpsy_wgcna_allgenes.csv", row.names=1)
head(vsd)
options(stringsAsFactors=FALSE)
data=t(vsd)
allkME =as.data.frame(signedKME(data, MEs))

whichModule="magenta" # name your color and execute to the end

length(moduleColors)
inModule=data.frame("module"=rep(0,nrow(vsd)))
row.names(inModule)=row.names(vsd)
genes=row.names(vsd)[moduleColors == whichModule]
inModule[genes,1]=1
sum(inModule[,1]) #should be same number as in darkgrey module
head(inModule)
row.names(inModule)=sub("", "isogroup", rownames(inModule))
write.csv(inModule,file=paste(whichModule,"magenta_fisher.csv",sep=""),quote=F)

##KME method. input for delta ranks
modColName=paste("kME",whichModule,sep="")
modkME=as.data.frame(allkME[,modColName])
row.names(modkME)=row.names(allkME)
names(modkME)=modColName
row.names(modkME)=sub("", "isogroup", rownames(modkME))
write.csv(modkME,file=paste(whichModule,"magenta_kME.csv",sep=""),quote=F)

#############################################################################
#greenyellow module

#Gene relationship to trait and important modules:
# Define variable weight containing the weight column of datTrait - leave weight as variable, but change names in first 2 commands
weight = as.data.frame(datTraits$Heat); #change Lipidrobust to your trait name
names(weight) = "Heat" #change based on the module you are running below
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr0, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="")
#Gene-trait significance correlation plots. change module to look at different ones
# par(mfrow=c(2,3))
module = "greenyellow"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("ModMem in", module, "module"),
                   ylab = "Gene Sig for Cool",
                   main = paste("MM vs. GS\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

#Making VSD files by module for GO plot functions
vs=t(datExpr0)
cands=names(datExpr0[moduleColors=="greenyellow"])  
#looking at genes in this module and subsetting them
c.vsd=vs[rownames(vs) %in% cands,]
head(c.vsd)
nrow(c.vsd) #should correspond to module size
table(moduleColors)
#moduleColors
#brown        cyan greenyellow     magenta   turquoise 
#4731         191        3176        2659        5205 
head(c.vsd)
#creates csv file with subsetted data
write.csv(c.vsd,"rlog_MMgreenyellow.csv",quote=F)

##############################heatmap of module expression with bar plot of eigengene, no resorting of samples...
#names(dis)
sizeGrWindow(8,7);
which.module="greenyellow" #pick module of interest
ME=MEs[, paste("ME",which.module, sep="")]
genes=datExpr0[,moduleColors==which.module ] #replace where says subgene below to plot all rather than just subset

#quartz()
# par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
par(mfrow=c(2,1), mar=c(0.3, 5.5, 5, 2))
plotMat(t(scale(genes) ),nrgcols=30,rlabels=F, clabels=rownames(genes), rcols=which.module)
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="sample")
#eigengene is overall expression
#this is a cool plot where you can see that genes in this module are upregulated in the pH7.6 treatment


#Gene relationship to trait and important modules: Gene Significance and Module membership
allkME =as.data.frame(signedKME(t(dat), MEs))
head(allkME)
vsd=read.csv(file="rlog_MMgreenyellow.csv", row.names=1)
head(vsd)

library(pheatmap)

############################################
#top 100 genes
whichModule="greenyellow"
top=100

datME=MEs
vsd <- read.csv("Bpsy_wgcna_allgenes.csv", row.names=1)
head(vsd)
datExpr=t(vsd)
modcol=paste("kME",whichModule,sep="")
head(vsd)
sorted=vsd[order(allkME[,modcol],decreasing=T),]
hubs=sorted[1:top,]
# attaching gene names
summary(hubs)



contrasting = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
#quartz()
pheatmap(hubs,scale="row",col=contrasting,border_color=NA, main=paste(whichModule,"top",top,"kME",sep=""))

###fisher for GO
# two ways to do GO: Fischers- every you in the module or not or KME: how well gene belongs, more continuous
##########fisher of module vs whole dataset
library(WGCNA)
vsd <- read.csv("Bpsy_wgcna_allgenes.csv", row.names=1)
head(vsd)
options(stringsAsFactors=FALSE)
data=t(vsd)
allkME =as.data.frame(signedKME(data, MEs))

whichModule="greenyellow" # name your color and execute to the end

length(moduleColors)
inModule=data.frame("module"=rep(0,nrow(vsd)))
row.names(inModule)=row.names(vsd)
genes=row.names(vsd)[moduleColors == whichModule]
inModule[genes,1]=1
sum(inModule[,1]) #should be same number as in darkgrey module
head(inModule)
row.names(inModule)=sub("", "isogroup", rownames(inModule))
write.csv(inModule,file=paste(whichModule,"greenyellow_fisher.csv",sep=""),quote=F)

##KME method. input for delta ranks
modColName=paste("kME",whichModule,sep="")
modkME=as.data.frame(allkME[,modColName])
row.names(modkME)=row.names(allkME)
names(modkME)=modColName
row.names(modkME)=sub("", "isogroup", rownames(modkME))
write.csv(modkME,file=paste(whichModule,"greenyellow_kME.csv",sep=""),quote=F)
######--------------------end--------------------#######
