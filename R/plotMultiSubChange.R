#' Plots heat map of multiple substitution associated changes per site
#'
#' @param multiSubChangeDF A data frame of multiple amino acid substitutions per site calculated with 'getAASub' and associated property changes obtained by calling 'getPropChange' function
#'
#' @return An interactive heat map of multiple substitution assocaiated changes per site
#' @export
#' @import plotly
#' @importFrom magrittr %>% 
#' @examples
#' multiSubChangeLoc <- system.file("extdata", "multiSubChange.rda", package = "aaSEA")
#' multiSubChange <- readRDS(multiSubChangeLoc)
#' plotMultiSubChange(multiSubChangeDF = multiSubChange )

plotMultiSubChange <- function(multiSubChangeDF){
 # library(plotly)
  msc <- multiSubChangeDF
  msco <- subset(msc, msc$wt != "-")
  msco <- subset(msco, msco$mu != "-")
  msco$xlab <- paste(msco$wt, msco$site, sep = "")
  
  # simple Heat map route
  plotly::plot_ly(
    x = as.character(msco$xlab),
    y = msco$mu,
    z = round(msc$Delta.Prop, 2),
    type = "heatmap",
    zauto = FALSE
  ) %>% layout(title = "Multiple substitutions and associated property changes")
}