#when re-running program reference genome doesn't have to be imported again
rm(list = ls()[ls() != "rGenome0"])

#MAIN
library(updog)
library(dplyr)
library(seqinr)
library(tibble)
library(ldsep)
library(purrr)
library(truncnorm)

# setwd("/media/rna/WHAKAHAO/Emily/Simulation/Simulation_true_positives/src/")
# setwd("C:/Users/Emily/EmilyGit/Uni/Masterarbeit/Stub_NGS_simulation/Simulation_true_positives/src/")
setwd("C:/Users/mimib/Documents/Uni/Masterarbeit/Stub_NGS_simulation/Simulation_true_positives/src/")



pathDat <- "../../../Data/"

#label for all output of simulation
labelRun <-  "Try2"
# labelRun <-  "TEST2"
VCFnam <- paste0("data4VCF_", labelRun)
outPath <- paste0(pathDat, "simulatedData/vcfSimulation/")

###################################################################################
# SIMULATON OF SNPS ON A CHROMOSOME OF A REFERENCE GENOME FOR SEVERAL INDIVIDUALS
# THE ALTERNATE ALLELE IS SAMPLED RANDOMLY (NO MULTIALLELIC SNPS)
# current runtime: Time difference of 17.7753 mins
#New plan: Haplotype blocks which are on average 100kbp, one "seed-SNP" is sampled via rgeno, the other SNPs within the haploblock are defined by it.
#Initial SNPs minor allele frequency is generated by gene density.
###################################################################################
t0 <- Sys.time()

# set.seed(40624)
set.seed(1234567)
#seed for fixing gene annotations
# set.seed(24445)


#IMPORTS PARAMETERS
source("parameters/parameters.R")
# source("parameters/testParameters.R")

#IMPORTS REFERENCE GENOME and it's annotation file and creates subset of them 
#according to simulation parameters
source("programs/importReference.R")


#ANNOTATION FILES: HELP FILES
#extracts data frame of genes and nested list of annotations, that are all in the 
#same gene from annotation file (output: annotGenes, annotList)
source("functions/annotation_functions.R")
source("programs/processAnnotationFile.R")

#-----------------------------------------------------------------------------------
#SOURCE SIMULATION

#returns haploLim (vector containing boundaries of haplotypes in format which can be used in default version of cut() function)
#haplotype boundaries are not within genes
source("programs/simulateHaplotypeBoundaries.R")


#simulates SNP positions for each haplotype block
#returns dat4VCF (final data frame) with its first two variables containing SNP positions and the corresponding haplotype Block
source("programs/simulateSnpPositions.R")

#adds reference allele and alternative allele to dat4VCF
source("functions/simulateAlleles_functions.R")
source("programs/simulateAlternativeAlleles.R")

#simulates dosage and phasing for all SNPs for all individuals
source("functions/simulateDosages_functions.R")
source("functions/simulatePhasing_functions.R")
source("programs/simulateGenotypes.R")


#EXPORT
#-----------------------------------------------------------------------------------
# backup of Haploblock number and snp type (correlated/uncorrelated), then removal from data set in prepatation of VCF-file generation
dat4VCF0 <- dat4VCF
dat4VCF <- dat4VCF[,!names(dat4VCF) %in% c("snpType", "haplotype")]
dat4VCF$POS <- as.integer(dat4VCF$POS)

# export data which will later be used in bash script to generate genotype
write.table(dat4VCF, file = paste0(outPath, VCFnam,  ".txt"), sep = "\t", quote = F, row.names = F,)

# save meta information for VCF
meta <- c(rGenomeName, chromosome, length(rGenome),c(paste("Indv", 1:nIndv, sep="_", collapse = "\t")))
names(meta) <- c("referenceGenome", "chromosomeName", "chromosomeLength","sampleNames")
write.table(meta, file = paste0(outPath, VCFnam,  "_meta.txt"), sep = "\t", quote = F, row.names = T, col.names = F)

#calculate time difference
timeDiff <- Sys.time()-t0
print(timeDiff)
parameters$runTime <- timeDiff
rm(t0)

#save Timing
timesave <-  tibble(label="dat4VCF in R", time=timeDiff)
write.table(timesave, file = paste0(pathDat, "simulatedData/timingSimulation_", labelRun,".txt"), row.names = F)

#save backup data for testing
save(dat4VCF0, test_dat4VCF_dosages, parameters, haploLim,file = paste0(outPath, VCFnam, ".RData"))


# save corresponding reference for contig/chromosome 
sequInd <- (min(dat4VCF$POS)):(max(dat4VCF$POS))
length(rGenome0[[chromosome]])
rGenome0_subset0 <- rGenome0[[chromosome]][sequInd]#[as.numeric(min(datVCF@fix[,"POS"])): as.numeric(max(datVCF@fix[,"POS"]))]
rGenome0_subset <- as.SeqFastadna(rGenome0_subset0, Annot = paste0(">", chromosome), name = chromosome)

setwd(paste0(outPath))
write.fasta(toupper(rGenome0_subset), names = chromosome, file.out = paste0(pathDat,"ReferenceGenomes/Genome_for_sim_",VCFnam,".fa"))