##Three Way Data Comparison by visualizing log2fold changes against a control.
##SET VARIABLES##
controlSample <- "PtcG4_D6"
liveSample1 <- "7AGFP_D5"
liveSample2 <- "7AGFP_D8"

#Run DESeq2 using 3 way comparison test (LRT)#
dds <- DESeqDataSetFromMatrix(rawdata, info, ~GenotypeDay)
ddsTimeSeries <- DESeq(dds, test = "LRT" , reduced = ~1)

#Write results of Control V Group 1 and Control V Group 2 to variables res1Table & res2Table
res1 <- results(ddsTimeSeries, contrast= c("GenotypeDay", controlSample, liveSample1))
write.csv(res1, "tempFile1.csv")
res1Table <- read.csv("tempFile1.csv")
file.remove("tempFile1.csv")
res2 <- results(ddsTimeSeries, contrast = c("GenotypeDay", controlSample, liveSample2))
write.csv(res2, "tempFile2.csv")
res2Table <- read.csv("tempFile2.csv")
file.remove("tempFile2.csv")

#Merge Data sets# 
resMerged <- merge(x = res1Table, y= res2Table, by = 0)

#Clean Merged Table#
resMergedClean <- subset(resMerged, select = c(X.x, log2FoldChange.x, log2FoldChange.y))
row.names(resMergedClean) <- resMergedClean$X.x
resMergedClean <- subset(resMergedClean, select = -c(X.x))
colnames(resMergedClean) <- c(paste0(liveSample1,"_log2FoldChange"), paste0(liveSample2,"_log2FoldChange"))
resMergedFinal <- na.omit(resMergedClean)
write.csv(resMergedFinal, paste0(liveSample1,"_&_",liveSample2,"_V_",controlSample,".csv"))


##SELECT GENE/S FOR ANALYSIS##
genesForAnalysis <- list("FBgn0261529", "FBgn0261581", "FBgn0261797") ## <- Paste Gene List Here

#Pull row.names data for all genes specified above#
wantedGenesTable <- resMergedFinal[paste0(genesForAnalysis), ]
#Melt Matrix#
wantedGenesTable <- melt(as.matrix(wantedGenesTable))
#Rename Columns#
colnames(wantedGenesTable) <- c("Gene_ID", "Genotype", "Log2FoldChange") ## <- May need to edit wanted column names.
#Rename X.axis labels#
wantedGenesTable$Genotype <- ifelse(grepl("D5", wantedGenesTable$Genotype), "7AD5","7AD8") ##May Need to edit wanted X axis label names.
#GG plot##
log2FoldChangePlot <- ggplot(wantedGenesTable, aes(x =Genotype, y =Log2FoldChange, color= Genotype)) +
  geom_point() +
  facet_grid(~Gene_ID) +
  theme(axis.text.x = element_text(angle = 90))
print(log2FoldChangePlot)