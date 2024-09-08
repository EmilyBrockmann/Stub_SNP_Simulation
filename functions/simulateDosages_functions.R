shift_dosage <- function(dosage, ploidy=4){
  # FUNCTION SHIFTS MARKER DOSAGE OF A SINGLE SNP OF A SINGLE INDIVIDUAL BY ONE
  # :param dosage: number in interval [0,1,...,4]
  # :param ploidy [int]: ploidy of organism
  # :returns: [int] number smaller or larger than dosage if ploidy allows it
  
  if(dosage == 0){
    return(1)
  }
  else if(dosage == ploidy){
    return(ploidy-1)
  }
  else{
    return(sample(c(dosage-1, dosage+1), 1))
  }
}

shift_random_individuals <- function(vec, x){
  # FUNCTION SHIFTS MARKER DOSAGE OF A SINGLE SNP OF A CERTAIN NUMBER OF RANDOMLY SELECTED INDIVIDUALS
  # :param vec: marker dosage for a SNP for all individuals ("seed SNP")
  # :param x: number of individuals of which the genotype is changed
  # :returns: marker dosage for a SNP for all individuals whith shifted dosage compared to reference vector "vec"
  if(x > length(vec)){
    x <- length(vec)
  }
  index <- sample(1:length(vec), replace = F, size = x)
  for(i in index){
    vec[i] <- shift_dosage(vec[i])
  }
  return(vec)
}

# 
# invert_dosages0 <- function(vec, ploidy=4, prob=1){
#   # :param vec: [nIndv x 1] marker dosage 
#   # :param ploidy: ploidy of organism
#   # :param prob: [double] probability of an individuals dosage being inverted
#   # :returns: vector of dosages for a single SNP for all individuals in which some dosages are inverted
#   dosages0 <- 0:ploidy
#   dosages <- cbind(dosages0, rev(dosages0))
#   vecInverted <- unlist(sapply(vec, function(x){
#     inv <- dosages[,2][dosages[,1] == x]
#     return(sample(c(x, inv), size=1, prob = c(1-prob, prob)))
#   }, simplify = "array"))
#   return(vecInverted)
# }


invert_dosages <- function(vec, ploidy=4, prob=1, returnInverted=F){
  # :param vec: [nIndv x 1] marker dosage 
  # :param ploidy: ploidy of organism
  # :param prob: [double] probability of an individuals dosage being inverted
  # :returns: vector of dosages for a single SNP for all individuals in which some dosages are inverted
  dosages0 <- 0:ploidy
  dosages <- cbind(dosages0, rev(dosages0))
  # print(dosages)
  invertedIndexes <- c()
  vecInverted <- vec
  for(i in 1:length(vec)){
    inv <- dosages[,2][dosages[,1] == as.numeric(vec[i])]
    # print(dosages[,1] == vec[i])
    #chooses if the SNP is inverted for this individual
    invertBool <- sample(c(FALSE, TRUE), size=1, prob = c(1-prob, prob))
    
    if(invertBool){
      vecInverted[i] <- inv
      invertedIndexes <- c(invertedIndexes, paste0("Indv",i))
    }
  }
  if(returnInverted){
    return(list(vecInverted, list(invertedIndexes)))
  }
  else{return(vecInverted)}
}


inversionCorrection <- function(tab, minDist=300){
  #FUNCTION SHIFTS WHICH POSITIONS ARE ASSIGNED AS "inv" so they have a minumum distance to other snps
  # :param tab: data frame which contains column POS (position of snps) and a column "snpType" ("corr", "inv", "uncorr")
  
  tab <- tab %>% as.data.frame(arrange(tab, POS))
  
  #distance to previous SNP
  tab0 <- tab %>% mutate(dist1 = c(NA, sapply(2:nrow(tab), function(x){
    return(tab$POS[x]-tab$POS[x-1])
  }, simplify = "array")))
  #distance to next SNP
  tab0 <- tab0 %>% mutate(dist2 = c(sapply(1:nrow(tab)-1, function(x){
    return(tab$POS[x+1]-tab$POS[x])
  }, simplify = "array"), NA))
  
}



