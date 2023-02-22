
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

