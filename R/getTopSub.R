#' Get sites with more than one-correlated substitutions other than conserved amino acids at that position
#'
#' @param selMat A subset matrix of original multiple sequence alignment with significant correlations identified with 'getCorSites' function
#'
#' @return A data frame with two columns i.e. Pos1 and Pos2 which is a filtered subset of many correlated substitutions based on frequency of substitution after consensus of a column in multiple sequence alignment
#' @export
#' @import seqinr
#' @examples
#' selMatLoc <- system.file("extdata", "selMat.rda", package = "aaSEA")
#' selMat <- readRDS(selMatLoc)
#' getTopSub(selMat = selMat)



getTopSub <- function(selMat){
  conserved <-consensus(selMat)
  posNames <- colnames(selMat)
  names(conserved) <- posNames
  positions <- colnames(selMat)
  lst <- list()
  cntr = 1
  for (i in 1:nrow(selMat)){
    seq = selMat[i,];
    for (j in positions){
      if(seq[j] != conserved[j]){
        for(k in positions){
          if(k > j && seq[k] != conserved[k]){
            pos1 <- paste(conserved[j], j,seq[j], "  ", conserved[k],k, seq[k], sep = "")
            lst[[cntr]] = pos1
            cntr = cntr+1
          }
        }
      }
    }
  }
  mydf <- as.data.frame(unlist(lst))
  colnames(mydf) <- "CorSub"
  
  # tyr to get substitution pairs discovered more than once
  filt <- subset(mydf, duplicated(mydf$CorSub) | duplicated(mydf$CorSub, fromLast = TRUE))
  uni <- unique(filt) # Unique of duplicated
  res <- as.data.frame(uni[!grepl("-", uni$CorSub),]) #exclude sites other than substitutions
  colnames(res) <- "subPair"
  resDf <- data.frame(do.call('rbind', strsplit(as.character(res$subPair),' ',fixed=TRUE)))
  resDf$X2 <- NULL
  colnames(resDf) <- c("Pos1", "Pos2")
  # could report substitutions whihc occured only once by inverting above procedure
  return(resDf)
}



