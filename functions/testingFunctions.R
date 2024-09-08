#DUMP----------------------------------------------------------------------------------------------------
testPrint <- function(x, var=c()){
  print(paste0("Test Nr. :",x))
  if(length(var) != 0){
    print(paste0("Variable value: ", var))
  }
}


get_noSNP <- function(data, string="0|0|0|0"){
  # help function that checks if any POS in  dat4VCF have the genotpye "0|0|0|0" for all individuals
  data <- data[, grep("Indv", names(data))]
  cols <- sapply(1:nrow(data), function(x) all(data[x,] == string))
  print(sum(cols))
}

#TEST DOSAGE--------------------------------------------------------------------------------------------------------
minorAlleleFrequencies <- function(tab, ploidy=4){
  # FUNCTION RETURNS MINOR ALLELE FREQUENCY FOR DATA SET OF DOSAGES WHICH HAS SNPS AS ROWS AND INDIVIDUALS AS COLUMNS
  mafs <- sapply(1:nrow(tab), function(x){
    freq <- sum(tab[x,])/(ncol(tab)*ploidy)
    return(ifelse(freq < 0.5, freq, 1-freq))})
  return(mafs)
}


linkage <- function(tab, ploidy=4){
  # :param tab: [nSNPs x nIndv]
  
  # nSNPs x nIndv x ploidy+1
  gpTmp <- array(0,dim = c(nrow(tab), ncol(tab), ploidy+1))
  for(i in 1:nrow(tab)){
    for(j in 1:ncol(tab)){
      gtVec <- rep(0.000000001, ploidy+1)
      # print(tab[i,j])
      gtVec[tab[i,j]+1] <- 1
      # print(gtVec)
      gpTmp[i,j,] <- gtVec
    }
  }
  return(ldfast(gpTmp, type = "r", se=F))
}
linkage0 <- function(tab){
  # :param tab: [nIndv x nSNPs]
  
  corr <- c()
  for(i in 1:ncol(tab)){
    for(j in 1:ncol(tab)){
      corr <- c(corr, cor(tab[,i], tab[,j])**2)
    }
  }
  return(corr)
}





#TEST PHASINGS---------------------------------------------------------------------------------------
phasings_toNumeric <- function(tab, split = "\\|", ploidy=4){
  # :param tab: nSNPxnIndv data frame of phasing information
  # :returns: nSNPx(ploidy*nIndv) data frame of numeric phasing information (1 column per homologue chromosome)
  numericTab <- c()
  nam <- c()
  for(j in 1:ncol(tab)){
    
    tmp <- t(sapply(tab[,j], function(x){
      return(as.numeric(unlist(strsplit(x, split = split))))
    }, simplify = "array"))
    nam <- c(nam,(paste0("Indv",j,"_Homol", 1:ploidy)))
    numericTab <- numericTab %>% cbind(tmp)
  }
  numericTab <- data.frame(numericTab, row.names = NULL)
  names(numericTab) <- nam
  return(numericTab)
}

haplotypes <- function(tab){
  # :param tab: nSNPx(ploidy*nIndv) data frame of numeric phasing information (1 column per homologue chromosome)
  # :returns: tabele of unique haplotypes
  duplicatedColumns <- duplicated(t(tab))
  uniqueColomns <- tab[,!duplicatedColumns]
  names(uniqueColomns) <- paste0("Haplotype", 1:ncol(uniqueColomns))
  print(paste0(ncol(uniqueColomns), " of ", 2**nrow(tab), " possible configurations in ", ncol(tab), " chromosomes"))
  # return(uniqueColomns)
  return(ncol(uniqueColomns))
}
