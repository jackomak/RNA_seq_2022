
#String Sets for heatmap.


#Path to rawdata file.
rawdataPath <- "./Data/ALL_Tissues_LFC_Database.xlsx"

#List of names to describe the availible genotypes for the heatmap in the GUI.
genotypesForAnalysisNames <- list("RasYki (D5)","RasYki (D8)","Feritin (D6)",
                                  "Feritin (D8)", "Feritin WT looking (D6)",
                                  "ImpL2 RNAi (D6)", "ImpL2 RNAi (D8)")

#List of genotypes stored in the raw data excel file that describe column headings for each genotype.
genotypesForAnalysisIDs <- list("RasYki_D5","RasYki_D8","Fer12OG_D6",
                                "Fer12OG_D8","Fer12WT_D6","ImpL2i_D6",
                                "ImpL2i_D8") 

#Select which genotypes to initially mark as checked on the GUI.
initallySelectedGenotypes <- list("RasYki_D5", "RasYki_D8", "Fer12OG_D6", "Fer12OG_D8")

#List of names to describe which heatmaps can be built and visualized on the GUI.
tissueNames <- list("FH 2022 Wing disc", "FH 2022 Salivary Gland", 
                    "FH 2022 Brain")

#List of sheet names in the raw data excel file availible for heatmap creation.
tissueValues <- list("Wingdisc", "Salivarygland", "Brain")

#Select which datasets to intially mark as checked by the GUI.
initallySelectedTissues <- list("Wingdisc", "Salivarygland", "Brain")

#Function hardcode for individual gene count viewer.
formatIndividualGenecells <- function(df){
  df$Genotype <- ifelse(grepl("PtcG4_D6", df$Genotype), "PtcG4_D6",
                 ifelse(grepl("RasYki_D5", df$Genotype), "RasYki_D5", 
                 ifelse(grepl("RasYki_D8", df$Genotype), "RasYki_D8",
                 ifelse(grepl("Fer12OG_D6", df$Genotype), "Fer12OG_D6",
                 ifelse(grepl("Fer12OG_D8", df$Genotype), "Fer12OG_D8",
                 ifelse(grepl("Fer12WT_D6", df$Genotype), "Fer12WT_D6",
                 ifelse(grepl("ImpL2i_D6", df$Genotype), "ImpL2i_D6", 
                 ifelse(grepl("ImpL2i_D8", df$Genotype), "ImpL2i_D8",
                 ifelse(grepl("Yw_D5", df$Genotype), "Yw_D5", "Yw_D8")))))))))
}

#factor levels for individual gene count viewer
brainLevels <- c("Yw_D5", "RasYki_D5", "RasYki_D8", "ImpL2i_D6", "ImpL2i_D8")
WDSGLevels <- c("PtcG4_D6", "RasYki_D5", "RasYki_D8", "Fer12OG_D6", "Fer12OG_D8", "Fer12WT_D6", "ImpL2i_D6", "ImpL2i_D8")