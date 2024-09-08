#PARAM
#-----------------------------------------------------------------------------------
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

nhaplo <- 12

# maxPOS <- 100000*nhaplo
maxPOS <- 10000*nhaplo

#TODO: move/check model
# alleleFrequency0 <- seq(0.05,0.5,length.out =10)
# The average minor allele frequency (MAF) of a variant was 0.14 (Uitdewilligen 2013)
meanMAF <- 0.14

#probability of if an individuals dosage of a SNP is changed as compared to the seed SNP (?)
changeIndvProbBinom <- 0.05

#probability of SNP in haplotype block being inverted or not correlated to seed SNP
pUncorrelated <- 0.005

# #removed Inversion of singele SNPs
# pInversion <- 0.0
# probInversion <- 0.00


invert.phasing <- FALSE
is.testrun <- TRUE
#-----------------------------------------------------------------------------------------------------------
# Get all object names in the current environment
all_objects <- ls()

# Initialize an empty list to store single value variables
parameters <- list()

# Iterate over each object name
for (obj_name in all_objects) {
  # Get the actual object
  obj <- get(obj_name)
  
  # Check if the object is a single value (atomic vector of length 1)
  if (is.atomic(obj) && length(obj) == 1) {
    # Add the object to the list with its name as the key
    parameters[[obj_name]] <- obj
  }
}
rm(all_objects, obj_name)