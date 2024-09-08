# 02 simulates positions-----------------------------------------
simulate_pos_within_haploblock <- function(haploBoundaries, haploIndex, density){
  # FUNCTION SAMPLES SNP POSITIONS
  # :param haploBoundaries: Vector of upper limit of positions within each haploblock, interval closed on the right as in cut()
  # :param haploIndex: Index of haploblock within interval 1:length(haploBoundaries)-1
  # :param density: density of SNPs in haploblock
  # 
  # :returns: positions of SNPs within the haploblock
  
  # positions included in haploblock
  haploblockPOS <- seq(haploBoundaries[haploIndex]+1, haploBoundaries[haploIndex+1], by = 1)
  
  # amount of SNPs
  nSNP <- round(length(haploblockPOS)*density)
  
  # positions of SNPs within haploblock
  snpPos <- sample(haploblockPOS, nSNP, replace = F)
  return(snpPos)
}

#initializes data frame from which vcf will be generated and adds SNP positions as the first column
dat4VCF <- tibble(POS=as.vector(unlist(sapply(1:nhaplo, function(x){
  haploStart <- haploLim[x]+1
  haploEnd <- haploLim[x+1]
  
  #calculates exon percentage within this haploblock
  exonPerc <- ifelse(nrow(rGenome_annot) != 0, exon_Percentage(rGenome_annot,haploStart, haploEnd), 0)
  # exonPerc <- exon_Percentage(rGenome_annot,haploStart, haploEnd)
  
  #generates SNP positions with densities according to density of exons in haploblock
  return(simulate_pos_within_haploblock(haploLim, x, density = exonPerc*densityExon+(1 - exonPerc)*densityIntron))
}, simplify = "array")))
)

# store haploblock as factor in data set
dat4VCF <- dat4VCF %>% mutate(haplotype=cut(POS, haploLim, labels = 1:nhaplo))

rm(simulate_pos_within_haploblock)