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

setwd("C:/Users/corey/Downloads/Rscriptanddata/New_Oc9/Gene_Expression/Symbiont-Thermal-Stress")


GO=read.table(file="BP_brownbrown_kME.csv", sep="\t", header=T)

ribo=filter(GO, term=="GO:0000302")
nrow(ribo[abs(ribo$value)<1,])
sigribo=ribo[abs(ribo$value)<1,]
colnames(rld)=c( "Control_1.1", "Control_1", "Control_2.1", "Control_2", "Control_3.1", "Control_3", "Control_4.1", "Control_4", "Cool_1", "Cool_2", "Cool_3", "Cool_4","Heat_1.1", "Heat_1", "Heat_2.1", "Heat_2", "Heat_3.1", "Heat_3", "Heat_4.1", "Heat_4")
hm=subset(rld, rownames(rld) %in% sigribo$seq)

annotations <- read.table("B_psygmophilum_isogroup_to_genename.tab",sep="\t",quote="", row.names=1)
hm1<-transform(merge(hm,annotations,by=0))
hm.z = data.matrix(hm)
hm.z = sweep(hm.z, 1L, rowMeans(hm.z), check.margin = FALSE)
hm.z.sx = apply(hm.z, 1L, sd)
hm.z = sweep(hm.z, 1L, hm.z.sx, "/", check.margin = FALSE)
hm.z = data.matrix(hm.z)
colour = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
heatmap.2(hm.z, col = colour, Rowv = TRUE, Colv = TRUE, scale = "row", 
          dendrogram = "both",
          trace = "none", labRow = hm1$gene,cexRow=.7,
)
dev.off()
