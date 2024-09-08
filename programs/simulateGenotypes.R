# 04 simulates genotypes-----------------------------------------
save(dat4VCF, annotGene, haploLim, file = paste0(outPath, "log/savePOS.RData"))

#space for genotypes in data set (empty matrix)
individuals <- paste0("Indv", 1:nIndv)

#storage of all individual's phasings for all SNPs
dat4VCF_phasings <- tibble(POS=numeric(), haplotype = character(), snpType=character(), !!!set_names(rep(list(character()), length(individuals)), individuals))
#temporary storage of dosage for testing purposes
test_dat4VCF_dosages <- tibble(POS=numeric(), haplotype = character(), snpType=character(), !!!set_names(rep(list(numeric()), length(individuals)), individuals))


#haplotype blocks
for(hapType in unique(dat4VCF$haplotype)){
  
  #start: simulates dosages of all individuals for all SNPs in this haploblock-------------------------------------------
  
  #get position of SNPs within haplotype
  snpHaplo <- dat4VCF$POS[dat4VCF$haplotype == hapType]
  
  #empty data frame for dosages in this haplotype
  datHaploDosage <- as.data.frame(matrix(,nrow = length(snpHaplo), ncol = nIndv))
  names(datHaploDosage) <- individuals
  datHaploDosage <- datHaploDosage %>% mutate(POS=snpHaplo, .before =1)
  rm(snpHaplo)
  
  #3 Types of SNPs: similar to seed SNP, similar to seed SNP but "inverted", uncorrelated to seed SNP
  # datHaploDosage <- datHaploDosage %>% mutate(snpType=sample(c("corr", "inv", "uncorr"), size=nrow(datHaploDosage), 
  #                                                            replace=TRUE, prob=c((1-pInversion-pUncorrelated), pInversion, pUncorrelated)),.before =2)
  # 
  datHaploDosage <- datHaploDosage %>% mutate(snpType=sample(c("corr", "uncorr"), size=nrow(datHaploDosage), 
                                                             replace=TRUE, prob=c((1-pUncorrelated), pUncorrelated)),.before =2)
  

  #chooses if allele in reference is minor allele
  #TODO: alleleFrequencies
  # freq <- alleleFrequency0[as.numeric(hapType)]
  if(is.testrun==TRUE){
    freq <- meanMAF
  }else{
    freq <- rtruncnorm(1, a=0.001,b=0.499, mean=meanMAF)
  }
  alleleFrequency <- sample(c(freq, 1-freq), 1)
  
  #seed SNP from which the others are generated
  seedSNP <- rgeno(nIndv, ploidy, model="hw", allele_freq = alleleFrequency)
  
  # inversionIndexes <- list()
  # inversionIndexesNames <- c()
  
  #GENERATES DOSAGE 
  for(pos in datHaploDosage$POS){
    
    #01 simulates correlating vector of marker dosages for "correlating" and "correlating and inverted snps"
    # if(datHaploDosage$snpType[datHaploDosage$POS == pos] %in% c("corr", "inv")){
    if(datHaploDosage$snpType[datHaploDosage$POS == pos] == "corr"){
      #selects amount of individuals in which dosage is varied for this SNP (using seed SNP as reference)
      nShift <- rbinom(1, size = nIndv, prob = changeIndvProbBinom)
      datHaploDosage[datHaploDosage$POS == pos ,individuals] <- shift_random_individuals(seedSNP, nShift)
    }
    
    # #02 inverts SNP dosages for "inverted" SNP type -> not used
    # if(datHaploDosage$snpType[datHaploDosage$POS == pos] == "inv"){
    #   inversionList <- invert_dosages(datHaploDosage[datHaploDosage$POS == pos ,individuals], prob = probInversion, returnInverted = T)
    #   # print(str(inversionList[[1]]))
    #   datHaploDosage[datHaploDosage$POS == pos ,individuals] <- inversionList[[1]]
    #   #indexes of inverted individuals for this SNP
    #   names(inversionList) <- c("dosage", as.character(pos))
    #   inversionIndexes <- inversionIndexes %>% append(inversionList[[2]])
    #   inversionIndexesNames <- c(inversionIndexesNames, as.character(pos))
    #   
    #   #03 simulates random dosages for uncorrelated SNPs
    # }
    else if(datHaploDosage$snpType[datHaploDosage$POS == pos] == "uncorr"){
      alleleFrequencyUncorr <- runif(1, min=0.001, max = 0.999)
      datHaploDosage[datHaploDosage$POS == pos ,individuals] <- rgeno(nIndv, ploidy, model="hw", allele_freq = alleleFrequencyUncorr)
    }
  }
  # names(inversionIndexes) <- inversionIndexesNames
  rm(nShift, pos, alleleFrequency)
  
  #end:   simulates dosages of all individuals for all SNPs in this haploblock-------------------------------------------
  
  
  
  #start: generates phasings for all SNPs and all individuals in this haplotype--------------------------------------------------------
  seedSNP_phasings <- sapply(rep(2, nIndv), function(x) return(sample_phasing(x, ploidy = ploidy)))
  
  #storage of phasings for all individuals and all SNPs in this Haplotype 
  datHaploPhasing <- datHaploDosage
  datHaploPhasing <- datHaploPhasing %>% mutate(haplotype = hapType, .before = 2)
  
  for(indv in individuals){
    
    #set of phasing options(one for every marker dosage for this individual)
    tmp <- sapply(0:ploidy, function(x) return(similar_phasing(seedSNP_phasings[which(individuals == indv)], nNew=x)), simplify = "array")
    # tmpInv <- rev(as.vector(sapply(tmp, function(x) invert_phasing(x))))
    
    #transforms all snps of datHaploDosage into phasings for this individual
    datHaploPhasing[[indv]] <- sapply(1:nrow(datHaploPhasing), function(snp){
      dosage <- datHaploDosage[snp,indv]
      # print(dosage)
      #01 same phasing per dosage for correlating markers
      if(datHaploPhasing$snpType[snp] == "corr"){
        return(tmp[dosage+1])
      }
      # #02 same phasing per dosage for correlating markers
      # #TODO
      # if(datHaploPhasing$snpType[snp] == "inv"){
      #   if(invert.phasing == TRUE){
      #     #checks if snp was inverted
      #     if(as.character(snp) %in% names(inversionIndexes) && indv %in% inversionIndexes[[as.character(snp)]]){
      #       return(tmpInv[dosage+1])
      #     }
      #     return(tmp[dosage+1])
      #   }
      #   else{return(tmp[dosage+1])}
      #   
      #   
      #   #03 simulates random phasing for uncorrelated SNPs
      # }
      else if(datHaploPhasing$snpType[snp] == "uncorr"){
        return(sample_phasing(dosage, ploidy))   }
    }, simplify="array")
    
    rm(tmp)
  }
  #end:  generates phasings for all SNPs and all individuals in this haplotype
  
dat4VCF_phasings <- dat4VCF_phasings %>% add_row(datHaploPhasing)
datHaploDosage <- datHaploDosage %>% mutate(haplotype=hapType)
test_dat4VCF_dosages <- test_dat4VCF_dosages %>% add_row(datHaploDosage)
rm(seedSNP, seedSNP_phasings, indv)  

if(!is.testrun && (as.numeric(hapType) %% 10) == 0){
  print(paste0("Haploblock ", hapType, " of ", length(unique(dat4VCF$haplotype)), " complete"))
}
if(((as.numeric(hapType) %% 1000) == 0)){
  save(dat4VCF_phasings, test_dat4VCF_dosages, file = paste0(outPath, "log/saveUntilHapblock_", hapType, ".RData"))
}
if(is.testrun){
  print(paste0("Haploblock ", hapType, " of ", length(unique(dat4VCF$haplotype)), " complete"))
  if(((as.numeric(hapType) %% 10) == 0)){
    save(dat4VCF_phasings, test_dat4VCF_dosages, file = paste0(outPath, "log/test_saveUntilHapblock_", hapType, ".RData"))
  }
}
}
save(dat4VCF_phasings, dat4VCF, test_dat4VCF_dosages, parameters, file = paste0(outPath, "log/saveAfterHapblockGeneration.RData"))
rm(hapType, datHaploDosage, datHaploPhasing, individuals, alleleFrequencyUncorr, freq)

dat4VCF <- dat4VCF %>% left_join(dat4VCF_phasings)
dat4VCF <- dat4VCF[order(dat4VCF$POS),]
