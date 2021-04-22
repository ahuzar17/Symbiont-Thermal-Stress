GO=read.table(file="CC_BvsD_GO_T.csv", header=T)

ribo=filter(GO, term=="GO:0005840")
nrow(ribo[abs(ribo$value)>1,])
sigribo=ribo[abs(ribo$value)>1,]
hm=subset(HFNHvsLFNHrld, rownames(HFNHvsLFNHrld) %in% sigribo$seq)

annotations <- read.table("BvsD_T.txt")
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
          margin = c(8,38))
