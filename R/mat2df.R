#function to flatten correlation matrix
mat2df <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    Pos1 = rownames(cormat)[row(cormat)[ut]],
    Pos2 = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}