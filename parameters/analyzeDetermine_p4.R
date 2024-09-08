setwd("c:/Users/mimib/Documents/Uni/Masterarbeit/Stub_NGS_simulation/Simulation_true_positives/src/")
library(ggplot2)
library(cowplot)
library(ldsep)
rm(list = ls()[ls() != "rGenome0"])

source("functions/testingFunctions.R")
pathDat <- "../../../Data/"
outPath <- paste0(pathDat, "simulatedData/vcfSimulation/Parameters/p4/")

# pUncorrVec <- seq(0.05,1,0.05)
pUncorrVec <- seq(0,0.02,0.001)


datComp <- tibble(p4=numeric(), mean=numeric(), median=numeric(), nHaplo=numeric())


corrList <- list()
mafList <- list()
for(index in seq_along(pUncorrVec)){
  print(index)
  load(paste0(outPath, "data4VCF_p4_", pUncorrVec[index], ".RData"))
  tabSNPs <- test_dat4VCF_dosages
  tabSNPs <- as.data.frame(tabSNPs[, names(tabSNPs)[grep("Indv", names(tabSNPs))]])
  
  nSNPS_inHaplotype <- nrow(tabSNPs)
  
  #linkage
  corrSeedSNP <- linkage(tabSNPs)$ldmat
  #upper diagonal of correlation matrix, r squared
  corrSeedSNP <- corrSeedSNP[upper.tri(corrSeedSNP)]
  corrSeedSNP <- as.vector(corrSeedSNP)**2
  
  # store linkage date in list for plotting
  corrSeedSNP <- tibble(corr=corrSeedSNP)
  corrList[[length(corrList)+1]] <- corrSeedSNP
  
  #MAF
  mafs <- tibble(maf = minorAlleleFrequencies(tabSNPs))
  mafList[[length(mafList)+1]] <- mafs
  
  #Haplotypes
  tmp <- dat4VCF0[,grep("Indv", names(dat4VCF0))]
  nHaplotypes <-  haplotypes(phasings_toNumeric(tmp))
  datComp <- datComp %>% add_case(p4=pUncorrVec[index], mean=round(mean(corrSeedSNP$corr),3), median=round(median(corrSeedSNP$corr),3), nHaplo=nHaplotypes)
  
  
}
save(corrList, datComp, mafList, file = paste0(outPath, "Analysis.RData"))
###############################################
# outPath <- paste0(pathDat, "simulatedData/vcfSimulation/Parameters/p4/")
# plot(datComp$nHaplo~datComp$p4)
# corrPlot <- plot_grid(plotlist = corrList[c(rep(TRUE,12), rep(FALSE,length(corrList)-12))], ncol=4)
# print(corrPlot)
# 
# corrPlot <- plot_grid(plotlist = corrList[c(rep(FALSE,12), rep(TRUE,length(corrList)-12))], ncol=4)
# print(corrPlot)
# 
# corrPlot <- plot_grid(plotlist = mafList[c(rep(TRUE,12), rep(FALSE,length(mafList)-12))], ncol=4)
# print(corrPlot)
# 
# corrPlot <- plot_grid(plotlist = mafList[c(rep(FALSE,12), rep(TRUE,length(mafList)-12))], ncol=4)
# print(corrPlot)
