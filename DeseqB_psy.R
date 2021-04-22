# source("http://bioconductor.org/biocLite.R")
# biocLite("DESeq2")

#set your working directory
setwd("C:/Users/corey/Downloads/Rscriptanddata/New_Oc9/Gene_Expression") #you will need to change to your own directory

###conduct array quality metrics to detect and remove outliers
library(DESeq2) #Version 1.30.1
library(affycoretools) #Version 1.62.0
library(arrayQualityMetrics) #Version 3.46.0
library(genefilter) #Version 1.72.1
library(Biobase) #Version 2.50.0
library(ggplot2) #Version 3.3.3
library(dplyr) #Version 1.0.5
library(tidyverse) #Version 1.3.0
library(pheatmap) #Version 1.0.12
library(vegan) #Version 2.5.7
library(ggrepel) #Version 0.9.1
library(RColorBrewer) #Version 1.1.2
library(gplots) #Version 3.1.1
library(VennDiagram) #Version 1.6.20

#########################################################################################
#Outliers
#read in counts 
countData <- read.table("B_psygmophilum_counts.txt")
head(countData)
length(countData[,1])
#28265

###rename the annoying name. Add isogroup in the row name
names(countData)=c( "Control_1.1", "Control_1", "Control_2.1", "Control_2", "Control_3.1", "Control_3", "Control_4.1", "Control_4", "Cool_1", "Cool_2", "Cool_3", "Cool_4","Heat_1.1", "Heat_1", "Heat_2.1", "Heat_2", "Heat_3.1", "Heat_3", "Heat_4.1", "Heat_4")
row.names(countData)=sub("", "isogroup", rownames(countData))
head(countData)

#Make sure you can succesfully read in the iso2gene tab separated file
gg=read.table("B_psygmophilum_sequenceID_to_isogroup.tab",sep="\t",quote="", row.names=1) 

#Check for outliers

####First set the working directory and create a conditions table

setwd("C:/Users/Corey/Downloads/Rscriptanddata/New_Oc9/Gene_Expression/outlier")
v=setwd("C:/Users/Corey/Downloads/Rscriptanddata/New_Oc9/Gene_Expression/outlier")
treat=c( "Control_1.1", "Control_1", "Control_2.1", "Control_2", "Control_3.1", "Control_3", "Control_4.1", "Control_4", "Cool_1", "Cool_2", "Cool_3", "Cool_4","Heat_1.1", "Heat_1", "Heat_2.1", "Heat_2", "Heat_3.1", "Heat_3", "Heat_4.1", "Heat_4")
g=data.frame(treat)
g
colData= g

####Create our model for us. How does gene expression vary by treatment
dds=DESeqDataSetFromMatrix(countData=countData,
                           colData = g,
                           design = ~treat)
####vsd normalizing data
vsd.ge=assay(vst(dds))
rl=vst(dds)
e=ExpressionSet(assay(rl), AnnotatedDataFrame(as.data.frame(colData(rl))))
arrayQualityMetrics(e,outdir=v,intgroup=c("treat"),force=T)
warnings()
    #dev.off() #for windows if there is a bug
    #double-click index.html to look at array quality to see if there are outliers
### outlier found in Cool1

 #####you only ever need to run the above code once. Outliers are decided at the beginning. 

##########################################################################################
#Identifying differentially expressed genes
###rewriting some of the above steps reading in data because outliers is only run one time

setwd("C:/Users/corey/Downloads/Rscriptanddata/New_Oc9/Gene_Expression")

###read in counts Same stuff again after outliers
countData <- read.table("B_psygmophilum_counts.txt")
head(countData)
length(countData[,1])
#28265

names(countData)=c( "Control_1.1", "Control_1", "Control_2.1", "Control_2", "Control_3.1", "Control_3", "Control_4.1", "Control_4", "Cool_1", "Cool_2", "Cool_3", "Cool_4","Heat_1.1", "Heat_1", "Heat_2.1", "Heat_2", "Heat_3.1", "Heat_3", "Heat_4.1", "Heat_4")
#row.names(countData)=sub("", "isogroup", rownames(countData))
##ignoring this above step for now because it screwed up the GO analysis
head(countData)

##total number of raw counts. Useful to report
totalCounts=colSums(countData)
totalCounts
barplot(totalCounts, col=c("green", "green", "green", "green", "green", "green", "green", "green", "blue", "blue", "blue", "blue", "red", "red", "red", "red", "red", "red", "red", "red"), ylab="raw counts")

#Control_1.1   Control_1 Control_2.1   Control_2 Control_3.1 
#1285858     1259276     1229172      806901      931442 
#Control_3 Control_4.1   Control_4      Cool_1      Cool_2 
#571218      728048      603137      484972      713239 
#Cool_3      Cool_4    Heat_1.1      Heat_1    Heat_2.1 
#1113555     1184420     1060861     1144023      983978 
#Heat_2    Heat_3.1      Heat_3    Heat_4.1      Heat_4 
#1154335     1154315      606975     1186170      684762  
min(totalCounts) #484972
max(totalCounts)  # 1285858

###As we did when checking for outliers, we create the same conditions table and deseq object
treat=c( "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Cool", "Cool", "Cool", "Cool","Heat", "Heat", "Heat", "Heat", "Heat", "Heat", "Heat", "Heat")
g=data.frame(treat)
g
colData<- g
dds<-DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~treat) #can only test for the main effects of site, pco2, temp

#one step DESeq. Can run ones that are way more complicated 
dds<-DESeq(dds)
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing

head(dds)
res<- results(dds)

###Look at dispersions plot. Want it to look like a hockey stick
plotDispEsts(dds, main="Dispersion plot Symbiont Stress")


#################################################################################################
# Cool vs control pairwise comparisons. Starting to look for DE genes
colData$cool<-factor(colData$treat, levels=c("Cool","Control")) #treatment first control second
##second term is the "control"
rescool <- results(dds, contrast=c("treat","Cool","Control"))
###how many FDR < 10%?
table(rescool$padj<0.05)
# 0.1= 7163
# 0.05= 5834
# 0.01= 4057

###p adjusted value accounts for comparsions across different treatments.Account for all the times you are comparing things

#####This is a summary of the number of genes up and down-regulated relative to the control at a p adjusted value of 0.1
summary(rescool)

##another way to look at it 
nrow(rescool[rescool$padj<0.05 & !is.na(rescool$padj),])  # Num significantly differentially expressed genes excluding the no/low count genes   #5834

### mean of normalized counts. far right is more highly expressed genes, 
###more likely to find DE genes a=among those highly expressed. hard to find when low coverage
plotMA(rescool, main="Cool vs Control")
plotMA(rescool, main="Cool vs Control", ylim=c(-5,5))

##put results in dataframe
results <- as.data.frame(rescool)
head(results)
##another way to look at up and down regulated genes
nrow(rescool[rescool$padj<0.1 & rescool$log2FoldChange > 0 & !is.na(rescool$padj),])
nrow(rescool[rescool$padj<0.1 & rescool$log2FoldChange < 0 & !is.na(rescool$padj),])
#UP in Cool 3963
#DOWN in Cool 3200

###write data in a table
write.table(rescool, file="Cool_2021.txt", quote=F, sep="\t")
###read the table back in
cd <- read.table("Cool_2021.txt")
head(cd)

##read it back in to make a GO table

###this is taking each isogroup and making ranked p value with directionality
cd
go_input_Cool = cd %>%
  tibble::rownames_to_column(var = "iso") %>%
  mutate(mutated_p = -log(pvalue)) %>%
  mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
  na.omit() %>%
  select(iso, mutated_p_updown)
head(go_input_Cool)
colnames(go_input_Cool) <- c("gene", "pval")
head(go_input_Cool)
write.csv(go_input_Cool, file="Cool_GO.csv", quote=F, row.names=FALSE)

#########################################################################################################
# Heat vs Control. Same thing as above but with the other treatment
colData$heat<-factor(colData$treat, levels=c("Heat","Control")) #treatment first control second
##second term is the "control"
resheat <- results(dds, contrast=c("treat","Heat","Control"))
#how many FDR < 10%?
table(resheat$padj<0.1)
# 0.1= 5415
# 0.05= 4330
# 0.01= 2860
								 
plotMA(resheat, main="Heat vs Control")
plotMA(resheat, main="Heat vs Control", ylim=c(-5,5))

results2 <- as.data.frame(resheat)
head(results2)

nrow(resheat[resheat$padj<0.1 & resheat$log2FoldChange > 0 & !is.na(resheat$padj),])
nrow(resheat[resheat$padj<0.1 & resheat$log2FoldChange < 0 & !is.na(resheat$padj),])
#UP in heat 2852
#DOWN in heat 2563

###write data in a table
write.table(resheat, file="Heat_2021.txt", quote=F, sep="\t")

cd2 <- read.table("Heat_2021.txt")
head(cd2)

##make the GO table for MWU for heat
head(cd2)

go_input_Heat = cd2 %>%
  tibble::rownames_to_column(var = "iso") %>%
  mutate(mutated_p = -log(pvalue)) %>%
  mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
  na.omit() %>%
  select(iso, mutated_p_updown)
head(go_input_Heat)
colnames(go_input_Heat) <- c("gene", "pval")
head(go_input_Heat)
write.csv(go_input_Heat, file="Heat_GO.csv", quote=F, row.names=FALSE)


###############################################################################################
##############################################################################
#--------------get pvals
# binding two columns and pulling out pvalues
valcool=cbind(rescool$pvalue, rescool$padj)
head(valcool)
colnames(valcool)=c("pval.cool", "padj.cool")
length(valcool[,1])
table(complete.cases(valcool))

valheat=cbind(resheat$pvalue, resheat$padj)
head(valheat)
colnames(valheat)=c("pval.heat", "padj.heat")
length(valheat[,1])
table(complete.cases(valheat))

##############################################################################################
#make rlogdata and pvals table
# normalization method, important for heatmap and PCA
rlog=rlogTransformation(dds, blind=TRUE) 
rld=assay(rlog)
head(rld)
colnames(rld)=paste(colData$treat)
head(rld)
length(rld[,1])

#bind together pvalues and rlog data 
rldpvals=cbind(rld,valcool, valheat)
head(rldpvals)
dim(rldpvals)
# [1] 28265    24
table(complete.cases(rldpvals))
#FALSE  TRUE 
#8768 19497

write.csv(rldpvals, "Bpsy2021_RLDandPVALS.csv", quote=F)

colnames(rld)=paste(colData$treat)
head(rld)

#############################################################################################
# Sample distance heatmap. Sample to sample. How similar are each sample
sampleDists <- as.matrix(dist(t(rld)))
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          margin=c(10, 10), main="Sample Distance Matrix")


#heat map of sample distances
#only use columns 1 through 20 to cut off pvalues becase we don't want to visulise them
rldpvals <- read.csv(file="Bpsy2021_RLDandPVALS.csv", row.names=1)
head(rldpvals)
rld=rldpvals[,1:20]
head(rld)

sampleDists <- dist(t(rld))
sampleDistMatrix <- as.matrix( sampleDists )
treat=c( "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Cool", "Cool", "Cool", "Cool","Heat", "Heat", "Heat", "Heat", "Heat", "Heat", "Heat", "Heat")
colnames(sampleDistMatrix)=paste(treat)
rownames(sampleDistMatrix)=paste(treat)

#Same heatmap as before but with color scales
heat.colors = colorRampPalette(rev(c("blue","yellow","red")),bias=0.3)(100)
pheatmap(sampleDistMatrix,color = heat.colors,cex=0.9,border_color=NA,cluster_rows=T,cluster_cols=T)



#################################################################################
# VENN Diagram to include both up and down regulated genes in common for PC02
### First we create a series of genes that are up and down-regulated from each of our earlier pairwise comparisons.
cool_up=row.names(rescool[rescool$padj<0.1 & !is.na(rescool$padj) & rescool$log2FoldChange>0,])
length(cool_up) 
#3963
cool_down=row.names(rescool[rescool$padj<0.1 & !is.na(rescool$padj) & rescool$log2FoldChange<0,])
length(cool_down) 
#3200
heat_up=row.names(resheat[resheat$padj<0.1 & !is.na(resheat$padj) & resheat$log2FoldChange>0,])
length(heat_up) 
#2852
heat_down=row.names(resheat[resheat$padj<0.1 & !is.na(resheat$padj) & resheat$log2FoldChange<0,])
length(heat_down) 
#2563

###overall number of genes with correct p value and no NA
cool=row.names(rescool[rescool$padj<0.1 & !is.na(rescool$padj),])
heat=row.names(resheat[resheat$padj<0.1 & !is.na(resheat$padj),])

###We can look at the overall number of up or down genes across the two samples without repeated isogroups. 
####UP. Total number of up genes without repeated isogroups
pdegs1_up=union(heat_up,cool_up)
length(pdegs1_up)
#6163

####DOWN. Total number of down genes without repeated isogroups
pdegs1_down=union(heat_down,cool_down)
length(pdegs1_down)
#5060

#ALL
pdegs1=union(cool,heat)
length(pdegs1)
#10112

###Venn diagram code
####do UP, DOWN, ALL. Can change list to plot different things
candidates=list("Cool"=cool_up, "Heat"=heat_up)

prettyvenn=venn.diagram(
  x = candidates,
  filename=NULL,
  col = "transparent",
  fill = c("coral2", "forestgreen"),
  alpha = 0.5,
  # label.col = c("darkred", "white", "darkgreen", "white", "white", "white", "blue4"),
  cex = 2.5,
  fontfamily = "sans",
  fontface = "bold",
  cat.default.pos = "text",
  cat.col = c("darkred", "darkgreen"),
  cat.cex = 2.5,
  cat.fontfamily = "sans",
  cat.dist = c(0.08, 0.08),
  cat.pos = 1
);
grid.draw(prettyvenn)


###########################################################################################
#PCA looks at distance between samples. Focusing on PCA 1 and 2
rld_t=t(rld)
pca <- prcomp(rld_t,center = TRUE, scale. = TRUE)
head(pca)
li <- pca$sdev^2 / sum(pca$sdev^2)
pc1v <- round(li[1] * 100, 1)
pc2v <- round(li[2] * 100, 1)
pca_s <- as.data.frame(pca$x)
head(pca_s)
pca_s <- pca_s[,c(1,2)]
pca_s$Samples = row.names(pca_s)
pca_s$treat=colData$treat
head(pca_s)

cbPalette <- c("darkgoldenrod2",  "darkolivegreen3", "dodgerblue3")
ggplot(pca_s, aes(PC1, PC2, color = treat, pch = treat)) +
  geom_point(size=3) +
  #  geom_text_repel(aes(label=Samples)) +
  scale_colour_manual(values=cbPalette)+
  theme_bw() +
  # geom_density2d(alpha=.5)+
  geom_polygon(alpha=.2)+
  xlab(paste0("PC1: ",pc1v,"% variance")) +
  ylab(paste0("PC2: ",pc2v,"% variance")) 
  head(pca)
  
# tests if distances is significantly different. Treatment is having signficant effect  
adonis(pca$x ~ treat, data = pca_s, method='eu', na.rm = TRUE)
          # Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
# treat      2     168902  84451  3.8999 0.31451  0.001 ***
# Residuals  17    368133  21655         0.68549          
# Total      19    537035                1.00000 


###########################################################################################
# Heatmaps for genes based on expression in each treatment

rldpvals <- read.csv(file="Bpsy2021_RLDandPVALS.csv", row.names=1)
head(rldpvals)
rld_site= rldpvals[,1:20]
head(rld_site)

#Heat map for top 100 DE in heat
topnum= 100
top100=head(rldpvals[order(rldpvals$padj.heat), ],topnum)

p.val=0.1
conds=top100[top100$padj.heat<=p.val & !is.na(top100$padj.heat),]
length(conds[,1])

exp=conds[,1:20]
means=apply(exp,1,mean) 
explc=exp-means 
colnames(explc)=c( "Control_1.1", "Control_1", "Control_2.1", "Control_2", "Control_3.1", "Control_3", "Control_4.1", "Control_4", "Cool_1", "Cool_2", "Cool_3", "Cool_4","Heat_1.1", "Heat_1", "Heat_2.1", "Heat_2", "Heat_3.1", "Heat_3", "Heat_4.1", "Heat_4")

ccol=colorRampPalette(rev(c("red","chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
col0=colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)

pheatmap(explc,cluster_cols=T,scale="row",color=col0, show_rownames = F)

# Heat map for top 100 DE in cool
topnum= 100
top100=head(rldpvals[order(rldpvals$padj.cool), ],topnum)

p.val=0.1
conds=top100[top100$padj.cool<=p.val & !is.na(top100$padj.cool),]
length(conds[,1])

exp=conds[,1:20]
means=apply(exp,1,mean) 
explc=exp-means 
colnames(explc)=c( "Control_1.1", "Control_1", "Control_2.1", "Control_2", "Control_3.1", "Control_3", "Control_4.1", "Control_4", "Cool_1", "Cool_2", "Cool_3", "Cool_4","Heat_1.1", "Heat_1", "Heat_2.1", "Heat_2", "Heat_3.1", "Heat_3", "Heat_4.1", "Heat_4")

ccol=colorRampPalette(rev(c("red","chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
col0=colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)

pheatmap(explc,cluster_cols=T,scale="row",color=col0, show_rownames = F)

##########################################################################################
##########################################################################################
#Heatmap for the genes in common

rldpvals <- read.csv(file="Bpsy2021_RLDandPVALS.csv", row.names=1)

p.val=0.1 
conds=rldpvals[rldpvals$padj.cool<=p.val & !is.na(rldpvals$padj.cool) & rldpvals$padj.heat<=p.val & !is.na(rldpvals$padj.heat),]
rld_data= conds[,c(1:20)]


means=apply(rld_data,1,mean) # means of rows
explc=rld_data-means # subtracting them
colnames(explc)=c( "Control_1.1", "Control_1", "Control_2.1", "Control_2", "Control_3.1", "Control_3", "Control_4.1", "Control_4", "Cool_1", "Cool_2", "Cool_3", "Cool_4","Heat_1.1", "Heat_1", "Heat_2.1", "Heat_2", "Heat_3.1", "Heat_3", "Heat_4.1", "Heat_4")

ccol=colorRampPalette(rev(c("red","chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
col0=colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)

pheatmap(explc,cluster_cols=T,scale="row",color=col0, show_rownames = F, main = "Common genes")
