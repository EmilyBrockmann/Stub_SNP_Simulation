# 03 simulates alternate alleles-----------------------------------------

#alphabet for simulation of alternate alleles
alphabet <- unique(rGenome)
alphabet <- alphabet[alphabet != "n"]

#add transition vs. transversion probabilities
substitutionMatrix <- transition_Transversion_probabilities(alphabet, ratio = transitTransvRatio)

#adds reference allele and alternative allele for each SNP position
#using transition vs. transversion probabilities
dat4VCF <- dat4VCF %>% cbind(add_Alleles(dat4VCF$POS, rGenome, alphabet, substitutionMatrix))


rm(substitutionMatrix,alphabet)

print("Simulation of Alternative Alleles complete")