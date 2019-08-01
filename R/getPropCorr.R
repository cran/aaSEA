#' Get amino acid property wise correlations of co-evolving columns of a multiple sequence alignment
#'
#' @param selMat A subset matrix of original multiple sequence alignment with significant correlations identified with 'getCorSites' function
#' @param propertyDF One of the amino acid property data frames. viz. Cruciani, Fasgai, Kidera, AAindex. Default is Cruciani properties
#' @param propertyIndex Specific property row number from the data frame of propertyDF
#'
#' @return A data frame of four columns viz. Pos1, Pos2, Cor and p Value. Results are filtered to find position pairs with correlations above 0.8 and below -0.8
#' @export
#' @importFrom Hmisc rcorr
#' @importFrom stats complete.cases
#' @examples
#' selMatLoc <- system.file("extdata", "selMat.rda", package = "aaSEA")
#' selMat <- readRDS(selMatLoc)
#' getPropCorr(selMat = selMat, propertyDF = "Cruciani", propertyIndex = 1) 


getPropCorr <- function(selMat,
                        propertyDF = "Cruciani",
                        propertyIndex = 1){
  
  if(propertyDF == "Cruciani"){
    propDF <- Cruciani
  } else if(propertyDF == "Fasgai"){
    propDF <- Fasgai
  } else if(propertyDF == "Kidera"){
    propDF <- Kidera
  } else {
    propDF <- AAindex
  }
  # encode with desired properties
  cCod <- matEncode(aliMat = selMat, pIndex = propertyIndex, propDf = propDF)
  cCod
  # caliculate correlations and p value
  
  cp <- rcorr(as.matrix(cCod))
  colnames(cp$r) <- colnames(selMat)
  rownames(cp$r) <- colnames(selMat)
  colnames(cp$P) <- colnames(selMat)
  rownames(cp$P) <- colnames(selMat)
  df <- mat2df(cp$r, cp$P)
  df <- df[complete.cases(df),]
  pfdf <- df[df$p <= 0.005,]
  cfdf <- subset(pfdf, pfdf$cor >= 0.8 | pfdf$cor <= -0.8)
  return(cfdf)
}
