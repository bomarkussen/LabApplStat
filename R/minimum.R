#' Minimum between factors
#' 
#' @description 
#' `minimum` finds the minimum of two factors, i.e. the finest factors that is coarser than both of the factors.
#' 
#' @keywords manip
#' 
#' @param x vector that will be interpreted as a factor.
#' @param y vector that will be interpreted as a factor.
#' @param concatenate.names boolean. If \code{TRUE} then the levels of the minimum are constructed as the concatenation of the levels for \code{x} and \code{y}. If \code{FALSE} then the levels of the minimum are given as numbers. Defaults to \code{TRUE}.
#' 
#' @return A factor with the minimum.
#' 
#' @author Bo Markussen
#'
#' @seealso
#' \code{\link{compare}}
#' 
#' @example 
#' x <- rep(c("boy","girl","adult","adult"),4)
#' y <- rep(c("child","child","man","woman"),4)
#' minimum(x,y)
#' minimum(x,y,FALSE)
#' 
#' @export
minimum <- function(x,y,concatenate.names=TRUE) {
  x <- as.character(factor(x))
  y <- as.character(factor(y))
  join <- 1:2  # initialize dummy so initial round is taken
  while (length(join)>1) {
    tmp <- table(x,y)
    for (ii in unique(x)) {
      join <- colnames(tmp)[tmp[ii,]>0]
      if (length(join)>1) y[is.element(y,join)] <- paste(join,collapse=".")
      if (length(join)>1) break
    }
    tmp <- table(x,y)
    for (ii in unique(y)) {
      join <- rownames(tmp)[tmp[,ii]>0]
      if (length(join)>1) x[is.element(x,join)] <- paste(join,collapse=".")
      if (length(join)>1) break
    }
  }
  if (concatenate.names) {
    return(factor(paste(x,y,sep=":")))
  } else {
    return(as.numeric(factor(paste(x,y,sep=":"))))
  }
}
