

#adds length of segment for each row in annotation file
rGenome_annot <- rGenome_annot %>% mutate(lengthSegm = abs(start-end)+1)

#sorts annotation file rowwise in order to make sure that a segment that lies within another is also listed afterwards (eg. gene is followed by exons)
rGenome_annot <- rGenome_annot %>% arrange(start, desc(end))

###########################################################
#transforms annotation into list of data frames (each of the latter contain sequences nested in one of the others, e.g. gene and corresponding CDS and exons)

#storage for list of data frames
annotList <- list()

#storage of individual data frames, initializing of first entry
group <- rGenome_annot[1,]
tmpList <- tibble(group)
i <- 2

#adds subsequent rows to each data frame if they are within the higher order row "group" (e.g. gene) and 
#starts a new data frame if this is not the case
#result: list of data frames (in each of which the first row contains all of the following)
while(i <= nrow(rGenome_annot)){
  vec <- rGenome_annot[i,]
  
  if(within_annotationRow(vec, group)){
    tmpList <- tmpList %>% add_row(vec)
  }
  else{
    #save list of entries that are within current group
    annotList <- annotList %>% append(list(tmpList))
    
    #new group
    group <- vec
    tmpList <- tibble(group)
  }
  i <- i+1
}
#add final data frame to list
annotList <- annotList %>% append(list(tmpList))

rm(i, vec, group, tmpList)


#extracts genes from annotList in order to select haplotype boundaries that are not within genes later
annotGene <- do.call(rbind,lapply(annotList, function(x){
  x[which.max(x$lengthSegm),]
}))

#tests if all annotations in annotGene are labeled "gene" and sorts anntGene
if(length(annotGene) != 0){
  annotGene <- annotGene[order(annotGene$start, decreasing = F),]
  if(any(annotGene$type != "gene")){
    warning("Not all entries in gene file are genes")
  }
}

#SHOULD NOT BE THE CASE DUE TO SORTING OF FILES
# #checks if one gene is within another completely
# if(length(annotList) != sum(rGenome_annot$type == "gene")){
#   warning("extracted genes are not the same as total amount of genes\n (see Following gene: ... Was summarized in ... message for details")
#   allGenes <- rGenome_annot %>% subset(type == "gene")
#   indMissing <- which(!allGenes$start %in% annotGene$start)
#   for(i in 1:nrow(annotGene)){
#     for(j in indMissing){
#       if(within_annotationRow(allGenes[j,], annotGene[i,])){
#         print("Following gene:")
#         print(allGenes[j,])
#         print("Was summarised in")
#         print(annotGene[i,])
#       }
#     }
#   }
# }
print("Annotation File Processing complete")

