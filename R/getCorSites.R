#' Get correlated sites with substitutions
#'
#' @param fileLoc exact location of multiple sequence alignment file in "FASTA" format
#' @param corMethod One of the methods to compute correlated sites viz. 'mip', 'elsc', 'mcbasc' and 'omes'. Default is 'mcbasc'.
#'
#' @return A subset alignment matrix of original multiple sequence alignment with significant correlations.
#' @export
#' 
#' @examples
#' file = system.file("extdata", "linB_toy_ali.fasta", package = "aaSEA")
#' getCorSites(fileLoc = file,corMethod="mcbasc")

getCorSites <- function(fileLoc, 
                        corMethod = "mcbasc"){
# caliculate Correlated sites
  resDF <- data.frame()
# read alignment for cor site selection
  myali <- import.fasta(fileLoc, aa.to.upper = TRUE, gap.to.dash = TRUE)
# read alignment for matrix generation
  myali1 <- read.alignment(fileLoc, format = "fasta", forceToLower = FALSE)
  mymat <- toupper(as.matrix(myali1))
# Calculate correlated sites
    if(corMethod == "mip"){
    mi <- mip(align = myali, diag = 0, gap_val = 0.5, z_score = TRUE)
    resDF <- mat2df(mi$gross, mi$normalized)
    colnames(resDF) <- c("Pos1", "Pos2", "MutualInfo", "Normalized_MI")
  } else if(corMethod == "elsc"){
    el <- elsc(align = myali, gap_val = 0.5, z_score = TRUE)
    resDF <- mat2df(el$gross, el$normalized)
    colnames(resDF) <- c("Pos1", "Pos2", "ELSC", "Normalized_ELSC")
  } else if(corMethod == "mcbasc"){
    mb <- mcbasc(align = myali, gap_val = 0.5, z_score = TRUE)
    resDF <- mat2df(mb$gross, mb$normalized)
    colnames(resDF) <- c("Pos1", "Pos2", "McBASC", "Normalized_McBASC")
  } else {
    om <- omes(align = myali, gap_val = 0.5, z_score = TRUE)
    resDF <- mat2df(om$gross, om$normalized)
    colnames(resDF) <- c("Pos1", "Pos2", "OMES", "Normalized_OMES")
  }
  resDF$Pos1 <- gsub(pattern = "[A-Z, a-z]",replacement = "", resDF$Pos1)
  resDF$Pos2 <- gsub(pattern = "[A-Z, a-z]",replacement = "", resDF$Pos2)
  
  fresDF <- subset(resDF, resDF[,4] > mean(resDF[,4]))
  fresDF <- fresDF[order(fresDF[,1]),]
# Collect sites for selelection
  p3 <- c(fresDF$Pos1,fresDF$Pos2)
  p3 <- sort(as.numeric(unique(p3)))
  sel <- mymat[,p3]
  colnames(sel) <- p3
  return(sel)
}