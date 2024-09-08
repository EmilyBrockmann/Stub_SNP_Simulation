setwd("c:/Users/mimib/Documents/Uni/Masterarbeit/Stub_NGS_simulation/Simulation_true_positives/src/")
library(ggplot2)
library(cowplot)
library(ldsep)
rm(list = ls()[ls() != "rGenome0"])

source("functions/testingFunctions.R")
pathDat <- "../../../Data/"
outPath <- paste0(pathDat, "simulatedData/vcfSimulation/Parameters/p1/")

# changeIndvProbBinomVec <- seq(0.05,1,0.05)
changeIndvProbBinomVec <- seq(0.0125,0.1,0.0125)
par(mfrow=c(ceiling(length(changeIndvProbBinomVec)/2),2))


datComp <- tibble(p1=numeric(), mean=numeric(), median=numeric(), nHaplo=numeric())


corrList <- list()
mafList <- list()
for(index in seq_along(changeIndvProbBinomVec)){
  load(paste0(outPath, "data4VCF_p1_", changeIndvProbBinomVec[index], ".RData"))
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
  datComp <- datComp %>% add_case(p1=changeIndvProbBinomVec[index], mean=round(mean(corrSeedSNP$corr),3), median=round(median(corrSeedSNP$corr),3), nHaplo=nHaplotypes)
  
  
  }
save(corrList, datComp, mafList, file = paste0(outPath, "Analysis.RData"))
###############################################
# outPath <- paste0(pathDat, "simulatedData/vcfSimulation/Parameters/p1/")
# 
# #pCorr <- ggplot(data=corrSeedSNP,aes(corr))+geom_histogram(alpha=0.3,bins=20, fill="cadetblue")+labs(title = paste0("r for p1 = ", changeIndvProbBinomVec[index]))
# #ggplot(data=mafs, aes(maf))+
# geom_histogram(alpha=0.3,bins=20, fill="cadetblue")+labs(title = paste0("maf for p1 = ", changeIndvProbBinomVec[index]))+
#   geom_vline(xintercept = 0.15)
# corrPlot <- plot_grid(plotlist = corrList, ncol=4)
# print(corrPlot)

