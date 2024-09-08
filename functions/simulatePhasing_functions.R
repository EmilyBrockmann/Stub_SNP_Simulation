sample_phasing <- function(nAltAllele, ploidy){
  # FUNCTION SAMPLES RANDOM PHASING FOR GIVEN NUMBER OF ALTERNATIVE ALLELES
  # :param nAltAllele [int]: genotype (number of alternative alleles) of individual that was sampled using rgeno()
  # :param ploidy [int]: ploidy of organism
  # :returns [string]: phased genotype of an individual, which are randomly sampled from all possible configurations
  
  #list of all posible genotypes (0|0|0|0, 1|0|0|0 etc.)
  phasedGT <- expand.grid(rep(list(0:1),ploidy))
  #selects genotypes that have requestet numbers of ones
  phasedGT <- phasedGT[rowSums(phasedGT)==nAltAllele,]
  #randomly selects one of the phasing options
  randomOption <- sample(1:nrow(phasedGT), 1)
  return(paste0(phasedGT[randomOption,], collapse = "|"))
}


invert_phasing <- function(phasing, split = "\\|"){
  #FUNCTION CHANGES ALL ONES INTO ZEROS AND ALL ZEROS TO ONES IN PHASING STRING
  # :param phasing: [str] phasing of SNP e.g. "0|1|1|0"
  # :param split: symbol seperating zeros and ones in phasing
  # :returns: [str] inverted phasing of SNP e.g. "1|0|0|1"
  
  #split up string into vector
  phasingVector <- as.numeric(unlist(strsplit(phasing, split = split)))
  #inverts vector
  phasingVectorInv <- sapply(phasingVector, function(x) ifelse(x==1,0,1))
  #transforms back into string
  return(paste0(phasingVectorInv, collapse = "|"))
}

similar_phasing <- function(seedSNPphasing, nNew, split = "\\|"){
  # FUNCTION GENERATES PHASING THAT IS SIMILAR TO THAT OF SEED SNP
  # :param seedSNPphasing: [str] phasing of seed SNP e.g. "0|1|1|0"
  # :param nNew: [int] marker dosage of new SNP (same individual)
  # :returns: [str] phased "similar" genotype of correlating SNP for a individual
  seedVec <- as.numeric(unlist(strsplit(seedSNPphasing, split = split)))
  if(sum(seedVec) == nNew){
    return(seedSNPphasing)
  }
  else if(sum(seedVec) > nNew){
    while(sum(seedVec) > nNew){
      # shifts random 1 to a 0 until genotype is the new desired one
      i <- sample(which(seedVec == 1), 1)
      seedVec[i] <- 0
    }
  }
  else if(sum(seedVec) < nNew){
    while(sum(seedVec) < nNew){
      # shifts random 0 to a 1 until genotype is the new desired one
      i <- sample(which(seedVec == 0), 1)
      seedVec[i] <- 1
    }
  }
  return(paste0(seedVec, collapse = "|"))
}