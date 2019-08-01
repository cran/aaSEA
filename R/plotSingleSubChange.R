#' Plots single substitution change histogram
#'
#' @param singleSubChangeDF A data frame of single amino acid substitutions per site calculated with 'getAASub' and associated property changes obtained by calling 'getPropChange' function
#'
#' @return An interactive histogram representing amino acid substitution associaited change
#' @export
#' @import plotly
#' @importFrom magrittr %>% 
#' @examples
#' singleSubChangeLoc <- system.file("extdata", "singleSubChange.rda", package = "aaSEA")
#' singleSubChange <- readRDS(singleSubChangeLoc)
#' plotSingleSubChange(singleSubChangeDF = singleSubChange)



plotSingleSubChange <- function(singleSubChangeDF){
  #library(plotly)
  ss <- singleSubChangeDF
  ss <- subset(ss, ss$wt != "-")
  ss <- subset(ss, ss$mu != "-")
  sso <- ss[order(ss$Delta.Prop), ]
  sso$substitution <- factor(sso$substitution,
                             levels = unique(sso$substitution)[order(sso$Delta.Prop,
                                                                     decreasing = TRUE)])
  m <-
    list(l = 50,
         r = 20,
         b = 50,
         t = 30) # l = left; r = right; t = top; b = bottom
  plotly::plot_ly(
    sso,
    x = ~ substitution,
    y = ~ Delta.Prop,
    type = "bar",
    color = ~ mu
  ) %>% layout(xaxis = list(tickangle = 45),margin = m,title = "Single substitutions and associated property changes")
}