transition_Transversion_probabilities <- function(nucleotides0, ratio = 1, prob = NA){
  # FUNCTION GENERATES NUCLEOTIDE SUBSTITUTION MATRIX (for sampling alternative alleles)
  # :param nucleotides0: nucleotides as in reference genome
  # :param ratio: transition/transversion ratio
  # :param prob: as an alternative to ratio the probability of a transition can be given directly
  # 
  # :returns: nucleotide substiution matrix
  
  #converts ratio into probability for transition, if necessary
  prob <- ifelse(is.na(prob), ratio/(1+ratio), prob)
  
  #initialize matrix for substitution probabilities (all initially transversion probabilities)
  tab <- matrix(1-prob,4,4)
  
  #converts alphabet to upper case in order for comparison to work
  nucleotides <- toupper(nucleotides0)
  rownames(tab) <- nucleotides
  colnames(tab) <- nucleotides
  
  #adds transition probabilities
  tab[c("A", "G"), c("A", "G")] <- prob
  tab[c("C", "T"), c("C", "T")] <- prob
  
  #sets diagonal to zero
  tab[cbind(c(1:4), c(1:4))] <- 0 
  
  #names set back into format of reference genome if necessary
  rownames(tab) <- nucleotides0
  colnames(tab) <- nucleotides0
  
  return(tab)
}


add_Alleles <- function(positions, referenceGenome, nucleotides, subMat){
  # FUNCTION GENERATES ALT AND REF ALLELES FOR GIVEN POSITIONS
  # gets allele from reference genome for each position and samples alternative allele
  # :param positions: positions of simulated SNPs on a chromosome of the reference genome
  # :param referenceGenome: nucleotide sequence of a certain chromosome or contig
  # :param nucleotides: alphabet of reference genome
  # :param subMat: 4x4 matrix containing substition probabilities for nucleotides
  # 
  # :returns: 2D tibble with column of reference alles and randomly sampled alternative alleles
  
  refAllele <- referenceGenome[positions]
  # randomly samples alternative Allele
  altAllele <- sapply(refAllele, function(x){
    if(x == "n"){return("n")}
    probs <- subMat[x,]
    return(sample(nucleotides, 1, prob = probs))
  })
  names(altAllele) <- c()
  return(tibble(refAllele=toupper(refAllele), altAllele=toupper(altAllele)))
}