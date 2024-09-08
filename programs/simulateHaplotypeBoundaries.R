# 01 simulate limits of haploblocks-----------------------------------------
if(nrow(annotGene) == 0){
  #if no genes/annotations exist, boundaries are sampled from all positions (used for testing)
  possibleSNPPos0 <- tibble(startPOS=1, endPOS=maxPOS)
}else if(annotGene$start[1] == 1 | annotGene$end[nrow(annotGene)] == maxPOS){
  stop("Invalid gene boundaries in annotGene")
} else{
  #generates (closed) interval boundaries in between genes (gene end included as intervals that are closed on the right are generated later using cut())
  possibleSNPPos0 <- tibble(startPOS=c(1,annotGene$end), endPOS=c(annotGene$start-1, maxPOS))
}

#generates vector from which haploblock boundaries can be sampled (from space in between genes)
possibleSNPPos <- unlist(sapply(1:nrow(possibleSNPPos0), function(x){
  return(possibleSNPPos0$startPOS[x]:possibleSNPPos0$endPOS[x])
}))


#samples upper limit of interval, interval closed on the right as in cut()
haploLim <- sort(sample(possibleSNPPos, nhaplo-1, replace = F))

rm(possibleSNPPos0, possibleSNPPos)

#reformating for cut() function (adding labels/haploblock number to each SNP-position later)
haploLim <-c(0,haploLim, maxPOS)

print("Simulate haploblock boundaries complete")
