within_annotationRow <- function(probe, higherOrderRow = NA, start = NA, end = NA){
  # FUNCTION CHECKS IF X IS IN INTERVAL OR ROW OF ANNOTATION FILE
  # :param probe: single number index or row(vector) of annotation file data frame, which contain columns named "start" and "end"
  # :param higherOrderRow: row(vector) of annotation file data frame, which contain columns named "start" and "end"
  # :param start: alternative to higherOrderRow, gives "start" parameter
  # :param end: alternative to higherOrderRow, gives "end" parameter
  # 
  # :returns: TRUE if row is within higherOrderRow according to limits of sequence given in "start" and "end", FALSE otherwise
  
  #extracts start and end parameters from row in question if row is given
  if(!all(is.na(higherOrderRow)) && is.na(start) && is.na(end)){
    # print("case row")
    start <- higherOrderRow[,"start"]
    end <- higherOrderRow[,"end"]
  }
  #checks for invalid input
  else if((all(is.na(higherOrderRow)) && is.na(start) && is.na(end)) | !all(is.na(higherOrderRow)) && (!is.na(start) | !is.na(end))){
    warning("Invalid input of within_annotationRow")
    return(NULL)
  }
  
  #tests for single numbers if they're in interval/higherOrderRow
  if(length(probe) == 1 && is.numeric(probe)){
    return(probe >= start & probe <= end)
  }
  #tests for rows of annotation file if they're in interval/higherOrderRow
  else{return(probe[,"start"] >= start & probe[,"end"] <= end)}
}

get_annotationRows <- function(annotData, index){
  # FUNCTION RETURNS ALL ANNOTATION ROWS WHICH CONTAIN INDEX
  # :param annotData: annotation file data frame, which contain columns named "start", "end" and "type"
  # :param index: sequence index in question
  # :returns: subset of annotData
  
  mask <- sapply(1:nrow(annotData), function(x) within_annotationRow(index, annotData[x,]), simplify = "array")
  return(annotData[mask,])
} 

exon_Percentage <- function(annotData, iStart, iEnd){
  # FUNCTION CALCULATES PERCENTAGE OF EXONS IN INTERVAL
  # ASSUMPTION: WHOLE EXON MUST BE IN INTERVAL
  # :param annotData: annotation file data frame, which contain columns named "start", "end" and "type"
  # :param iStart: starting index of interval for which exon ratio is to be calculated
  # :param iEnd: last index of interval for which exon ratio is to be calculated
  # 
  # :returns: percentage of exons in interval
  
  #selects rows that describe exons
  exons <- annotData %>% filter(type=="exon")
  
  #shouldn't be the case for haplotypes sampled in simulateGrnotypes.R
  if(nrow(get_annotationRows(annotData, iEnd)) != 0){
    tmp <- get_annotationRows(annotData, iEnd)
    if(any(iEnd != tmp$end)){
      print(paste0(iStart, " ",iEnd))
      warning("End within gene, exon not considered in calculation")
    }
  }
  
  if(nrow(get_annotationRows(annotData, (iStart-1))) != 0 ){
      print(paste0(iStart, " ",iEnd))
      warning("Start within gene, exon not considered in calculation")
  }
  
  #selects rows of exon that are in given interval completely
  mask <- sapply(1:nrow(exons), function(x){
    within_annotationRow(exons[x,], start = iStart, end = iEnd)
  }, simplify = "array")
  exons <- exons[mask,]
  
  #generates length of segments if necessary
  if(!"lengthSegm" %in% names(exons)){
    exons <- exons %>% mutate(lengthSegm = abs(start-end)+1)
  }
  return(sum(exons$lengthSegm, na.rm = T)/(abs(iEnd-iStart)+1))
}

#breaks if more that one haplotype boundary within one gene
# haplotype_limitCorrection <- function(index, genes, maxPosition){
#   #FUNCTION CORRECTS A HAPLOTYPE LIMIT IF IT'S WITHIN A GENE
#   # :param index: integer, limit of haploblock which will be corrected if neccessary
#   # :param genes: data frame that contains the columns "start" and "end" 
#   # 
#   # :returns: index or random sample of position befor or after gene (or other row of "genes"), if index lies within one
#   
#   #gets rows of genes which contain index
#   gene <- get_annotationRows(genes, index)
#   if(nrow(gene) == 1){
#     if(index == gene$start | index == gene$end){
#       return(index)
#     }
#     else if(index != 1 && index != maxPosition){
#       corr <- sample(c(gene$start-1, gene$end), 1)
#       print(paste0("Correction of Haplotype Limits: ", index, " -> ", corr))
#       print(get_annotationRows(genes, index))
#     return(corr)
#     }else{
#       warning("Sequence begins or ends with gene, which breaks limit correction function")
#       return(NULL)
#     }
#   }
#   else if(nrow(gene) > 1){
#     warning(paste0("gene not unique at index ", index))
#     return(NULL)
#   }else{return(index)}
# }
