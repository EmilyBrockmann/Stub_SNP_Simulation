#DATA
#-----------------------------------------------------------------------------------
rGenomeName <- "He1_haplotype_genome.fasta.gz"

#Sources reference genome only if it isn't in environment yet
if(!"rGenome0" %in% ls()){
  rGenome0 <- read.fasta(paste0(pathDat,"ReferenceGenomes/", rGenomeName))
}

#import annotation file and add column names
rGenome_annot0 <- read.delim(paste0(pathDat,"ReferenceGenomes/He1.protein-coding.gene.gff3"), sep = "\t", 
                             comment.char = "#", header = F)
names(rGenome_annot0) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", 
                           "attributes")

#filter for chromosome
rGenome_annot <- rGenome_annot0 %>% filter(seqid==chromosome)
rGenome_annot <- rGenome_annot[,c("type", "start", "end")]
rGenome <- rGenome0[[chromosome]]

#subset for testing if it is a testrun, else:calculation of number of bp and amount of haplotypes to be generated
if(is.testrun){
  rGenome <- rGenome0[[chromosome]][1:maxPOS]
}else{
  #when not given by testing parameters: calculation of number of bp and amount of haplotypes to be generated
  maxPOS <- length(rGenome)
  nhaplo <- ifelse(maxPOS > averageHaplotypeLength, floor(maxPOS/averageHaplotypeLength),2)
}
rGenome_annot <- rGenome_annot %>% filter(start <= maxPOS & end <= maxPOS)

#-----------------------------------------------------------------------------------