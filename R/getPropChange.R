#' Get wild type and substituted amino acid properties and associated property changes
#'
#' @param subFile A data frame of single or multiple substitutions obtained using 'getAASub' function
#' @param propertyDF Choose one of Cruciani, Fasgai, Kidera or AAindex based amino acid properties
#' @param propertyIndex The amino acid property to analyse. It is row number in propertyDF data frame
#'
#' @return A substitution data frame with three additional columns i.e. wt.Prop, mu.Prop and Delta.Prop 
#' @export
#'
#' @examples
# load multiple substitutions
#' ssFileLoc <- system.file("extdata", "singleSub.rda", package = "aaSEA")
#' singleSubFile <- readRDS(ssFileLoc)
# load multiple substitutions
#' msFileLoc <- system.file("extdata", "multiSub.rda", package = "aaSEA")
#' multiSubFile <- readRDS(msFileLoc)
#' getPropChange(subFile = singleSubFile, propertyDF = "Cruciani", propertyIndex = 1)
#' getPropChange(subFile = multiSubFile, propertyDF = "Cruciani", propertyIndex = 1)


getPropChange <- function(subFile, 
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
  
  subFile$wt.Prop <- c()
  subFile$mu.Prop <- c()
  subFile$Delta.Prop <- c()
  # caliculate wt Property
  for(i in 1:length(subFile$wt)){
    for (j in colnames(propDF)){
      if (subFile$wt[i]==j){
        subFile$wt.Prop[i]<-propDF[propertyIndex,j]
      }
    }
  }
  # caliculate Mu property
  for(m in 1:length(subFile$mu)){
    for (j in colnames(propDF)){
      if (subFile$mu[m]==j){
        subFile$mu.Prop[m]<-propDF[propertyIndex,j]
      }
    }
  }
  subFile$Delta.Prop <- subFile$mu.Prop - subFile$wt.Prop
  return(subFile)
}