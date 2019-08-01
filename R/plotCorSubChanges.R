#' Plot co-evolving sites with selected property correlations 
#'
#' @param corSitePropChange A data frame of coevolving sites and associated property changes obtained by "getPropCorr" function with selected MSA matrix (selMat) and desired property selected from property data frame and property index. 
#'
#' @return Returns an interactive heat map of significant sites with selected property correlations 
#' @export
#' @import plotly
#' @importFrom magrittr %>% 
#' @examples
#' fileLocation <- system.file("extdata", "corSitePropChangeDF.rda", package = "aaSEA")
#' corSitePropChange <- readRDS(fileLocation)
#' plotCorSubChanges(corSitePropChange = corSitePropChange)



plotCorSubChanges <- function(corSitePropChange){
  mydf <- corSitePropChange
  plotly::plot_ly(
    x = as.character(mydf$Pos1),
    y = as.character(mydf$Pos2),
    z = round(mydf$cor, 2),
    type = "heatmap",
    zauto = FALSE
  ) %>% layout(title = "Selected property correlations between coevolving sites")
}
