#' Get amino acid substitutions from multiple sequence alignment
#'
#' @param fileLoc exact location of multiple sequence alignment file in "FASTA" format
#'
#' @return Returns a list of two data frames 1.Single substitutions 2.Multiple substitutions
#' @export
#' @import Bios2cor seqinr
#' @examples
#' file = system.file("extdata", "linB_Prot_ali.fasta", package = "aaSEA")
#' getAASub(fileLoc = file)

getAASub <- function(fileLoc){
  res <- list()
  myali <- read.alignment(file = fileLoc
                          , format = "fasta", forceToLower = TRUE)
  mymat <- as.matrix(myali)
  mymat <- toupper(mymat)
  mycon <- consensus(mymat)
  for (i in 1 : ncol(mymat)){
    uniAA <- unique(mymat[,i])
    uniAA <- uniAA[! uniAA %in% c("-","X", mycon[i])]
    if (length(uniAA) >= 1)
    {
      res[[i]] <- paste(mycon[i],i,uniAA,sep = "")
    }
  } 
  resDF <- data.frame(unlist(res))
  colnames(resDF) <- "substitution"
  rownames(resDF) <- NULL
  
  # now process results
  resDF$wt <- substr(resDF$substitution, 1, 1)
  resDF$mu <- substr(resDF$substitution, nchar(as.character(resDF$substitution)), nchar(as.character(resDF$substitution)))
  resDF$site <- gsub(pattern = "[A-Z, -]", replacement = "", resDF$substitution)
  multiSub <- subset(resDF, duplicated(resDF$site) | duplicated(resDF$site, fromLast = TRUE))
  dups <- unique(multiSub$site)
  singleSub <- subset(resDF, !(resDF$site %in% dups))
  
  return(list("singleSub" = singleSub , 
              "multiSub" = multiSub,
              "All" = resDF))
}