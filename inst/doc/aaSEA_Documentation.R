## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(aaSEA)
library(plotly)
library(networkD3)
library(knitr)

## ------------------------------------------------------------------------
myMSA <- system.file("extdata", "linB_Prot_ali.fasta", package = "aaSEA") 
myMSA

## ------------------------------------------------------------------------
sub <- getAASub(myMSA)
str(sub)

## ------------------------------------------------------------------------
head(sub$singleSub)

## ------------------------------------------------------------------------
head(sub$multiSub)

## ------------------------------------------------------------------------
propChanges <- getPropChange(subFile = sub$singleSub, propertyDF = "Cruciani", propertyIndex = 1)
head(propChanges)

## ----fig1, fig.height = 5, fig.width = 9, warning=FALSE------------------
plotSingleSubChange(singleSubChangeDF = propChanges)

## ----fig2, fig.height=5, fig.width= 9------------------------------------
multiPropChange <- getPropChange(subFile = sub$multiSub, propertyDF = "cruciani", propertyIndex = 2)
plotMultiSubChange(multiSubChange = multiPropChange)

## ----message=FALSE, warning=FALSE, paged.print=FALSE---------------------
corSites <- getCorSites(fileLoc = myMSA, corMethod = "mcbasc") # can use any of the "omes", "elsc", "mip"
dim(corSites)

## ------------------------------------------------------------------------
corPairs <- getTopSub(corSites)
head(corPairs)

## ----fig3, fig.height=5, fig.width= 9------------------------------------
  plotCorNet(corPairs)

## ----fig.width=9, paged.print=FALSE--------------------------------------
coPropChanges = getCorPropChange( corSubFile = corPairs, propertyDF = "Cruciani", propertyIndex = 1)
head(coPropChanges)

## ----fig4, fig.height=5, fig.width= 9------------------------------------
propCor <- getPropCorr(selMat = corSites, propertyDF = "curciani",propertyIndex = 1)
head(propCor)
plotCorSubChanges(propCor)

