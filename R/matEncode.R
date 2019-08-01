#' Function to encode Correlated columns of alignment matrix with desired properties
#'
#' @param aliMat amino acid multiple sequence alignment in the form of a matrix
#' @param pIndex Amino acid property index to be encoded. It is the row number in the property data frame
#' @param propDf The amino acid property to analyse. It is row number in propertyDF data frame 
#'
#' @return A matrix of input dimensions with amino acid alphabets replaced by amino acid properties of choice
#' @export
#'
#' @examples
#' aliMatLoc <- system.file("extdata", "aliMat.rda", package = "aaSEA")
#' aliMat <- readRDS(aliMatLoc)
#' matEncode(aliMat = aliMat, pIndex = 1, propDf = "Cruciani" )



matEncode <- function(aliMat, pIndex, propDf){
  mycm <- matrix(data = 0, nrow = nrow(aliMat), ncol = ncol(aliMat))
  pIndex <- pIndex
  for (s in 1:nrow(aliMat)){
    mySeq <- aliMat[s,]
    for (i in 1:length(mySeq)){
      for (j in colnames(propDf)){
        if( j == toupper(mySeq[i])){
          mycm[s,i] <- as.numeric(propDf[pIndex, j])
          break
        } 
        else 
        {
          mycm[s,i] <- 0
        }
      }
    }
  }
  return(mycm)  
}