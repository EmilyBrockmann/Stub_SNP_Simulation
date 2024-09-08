#Simulates data for 1000 SNPs with different values for p4
#!!!!!!!!!!!!!!!!!!!!!
VCFnam0 <- "data4VCF_p4"
pathDat <- "../../../Data/"
outPath <- paste0(pathDat, "simulatedData/vcfSimulation/Parameters/p4/")

#p4<0.02
#set p4=0.005
#PARAM p4
#-----------------------------------------------------------------------------------
set.seed(40624)
nIndv <- 100
ploidy <- 4

#TODO Choose Chromosome
#format as in vcf
# chromosome <- "chr05"
chromosome <- "Chr5"

#(Uitdewilligen 2013)
transitTransvRatio <- 1.55  
densityExon <- 1/24
densityIntron <- 1/15

nhaplo <- 1

maxPOS <- 15000*nhaplo

# The average minor allele frequency (MAF) of a variant was 0.14 (Uitdewilligen 2013)
meanMAF <- 0.14

#probability of if an individuals dosage of a SNP is changed as compared to the seed SNP (?)
changeIndvProbBinom <- 0.05

#probability of SNP in haplotype block being inverted or not correlated to seed SNP
pInversion <- 0.0

# (0.05 already too large if changeIndvProbBinom <- 0.05)???
pUncorrelatedVec <- seq(0,0.02,0.001)

probInversion <- 0

invert.phasing <- TRUE
is.testrun <- TRUE

for(p4 in seq_along(pUncorrelatedVec)){
  pUncorrelated <- pUncorrelatedVec[p4]
  VCFnam <- paste0(VCFnam0, "_", pUncorrelated)
  
  
  
  #SAVE PARAMETERS
  # Get all object names in the current environment
  allObjects <- ls()
  
  # Initialize an empty list to store single value variables
  parameters <- list()
  
  # Iterate over each object name
  for (objNames in allObjects) {
    # Get the actual object
    obj <- get(objNames)
    
    # Check if the object is a single value (atomic vector of length 1)
    if (is.atomic(obj) && length(obj) == 1) {
      # Add the object to the list with its name as the key
      parameters[[objNames]] <- obj
    }
  }
  rm(allObjects, objNames)
  
  rm(list = ls()[!ls() %in% c("rGenome0", "parameters", names(parameters), "pUncorrelatedVec")])
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
  #save backup data for testing
  dat4VCF0 <- dat4VCF
  save(dat4VCF0, test_dat4VCF_dosages, parameters, haploLim,file = paste0(outPath, VCFnam, ".RData"))
}




#-----------------------------------------------------------------------------------------------------------
