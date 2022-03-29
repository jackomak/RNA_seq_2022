##Three Way Data Comparison by visualizing log2fold changes against a control.

###LIBRARIES###
rawdata <- read.csv("<RAWCOUNTS_FILE>", header = T, row.names = 1 ) # <- Rawcounts matrix file name.
info <- read.csv("<INFO_FILE>", header = T, row.names = 1) # <- Info matrix file name.
experimentalGroup <- "<EXPERIMENTAL_NAME>" # <- What does the rawcounts file contain? eg. "All_WD_SAMPLES".
experimentalCondition <- "<CONDITION_TO_TEST>" # <- Which experimental condition to test present on <INFO_FILE> eg. GenotypeDay.
controlSample <- "<CONTROL_SAMPLE>" # <- Control genotype from rawcounts matrix file for two way comparison.
liveSample1 <- "<LIVE_SAMPLE>" # <- Genotype #1 for analysis.
liveSample2 <- "<LIVE SAMPLE2>" # Genotype #2 for analysis.
#Select genes for analysis#
genesForAnalysis <- list("FBgn0261529", "FBgn0261581", "FBgn0261797") ## <- Paste Gene List Here

###RUN DESEQ2### using 3 way comparison test (LRT)###
dds <- DESeqDataSetFromMatrix(rawdata, info, ~GenotypeDay) # <- may need to change ~ to experimental condition here.
ddsTimeSeries <- DESeq(dds, test = "LRT" , reduced = ~1)

###WRITE RESULTS TO VARIBLES###
res1 <- results(ddsTimeSeries, contrast= c("GenotypeDay", controlSample, liveSample1))
write.csv(res1, "tempFile1.csv")
res1Table <- read.csv("tempFile1.csv")
file.remove("tempFile1.csv")
res2 <- results(ddsTimeSeries, contrast = c("GenotypeDay", controlSample, liveSample2))
write.csv(res2, "tempFile2.csv")
res2Table <- read.csv("tempFile2.csv")
file.remove("tempFile2.csv")

###MERGE AND CLEAN DATASETS### 
resMerged <- merge(x = res1Table, y= res2Table, by = 0)
resMergedClean <- subset(resMerged, select = c(X.x, log2FoldChange.x, log2FoldChange.y))
row.names(resMergedClean) <- resMergedClean$X.x
resMergedClean <- subset(resMergedClean, select = -c(X.x))
colnames(resMergedClean) <- c(paste0(liveSample1,"_log2FoldChange"), paste0(liveSample2,"_log2FoldChange"))
resMergedFinal <- na.omit(resMergedClean)
write.csv(resMergedFinal, paste0(liveSample1,"_&_",liveSample2,"_V_",controlSample,".csv"))

###PULL EXPRESSION DATA FOR WANTED GENES###
wantedGenesTable <- resMergedFinal[paste0(genesForAnalysis), ]
#Melt Matrix#
wantedGenesTable <- melt(as.matrix(wantedGenesTable))
#Rename Columns#
colnames(wantedGenesTable) <- c("Gene_ID", "Genotype", "Log2FoldChange") ## <- May need to edit wanted column names.
#Rename X.axis labels#
wantedGenesTable$Genotype <- ifelse(grepl("D5", wantedGenesTable$Genotype), "7AD5","7AD8") ##May Need to edit wanted X axis label names.

###PLOT LOG2FOLDCHANGE CHANGE AGAINST COMMON CONTROL###
log2FoldChangePlot <- ggplot(wantedGenesTable, aes(x =Genotype, y =Log2FoldChange, color= Genotype)) +
  geom_point() +
  facet_grid(~Gene_ID) +
  theme(axis.text.x = element_text(angle = 90))
print(log2FoldChangePlot)
