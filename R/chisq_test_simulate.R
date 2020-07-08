#' Simulate Chi-squared tests with conditioning
#'
#' @description
#' `chisq.test.simulate` simulates the chi-squared test for a 2-way contingency tabel.
#'
<<<<<<< HEAD
#' @keywords htest
#' 
=======
>>>>>>> 0c9bab0686cb3f7b4ca3129c4adb4aee148c9624
#' @param x matrix with the contingency table
#' @param conditioning character string specifying the simulation scenario. Defaults to \code{"total"}. Other possible scenarios are \code{"row"}, \code{"col"}, and \code{"both"}.
#' @param x0 matrix specifying the null distribution. Defaults to \code{NULL}, in which case the null is estimated from the observed data \code{x}.
#' @param B integer specifying the number of replicates used in the Monte Carlo test. Defaults to 10000.
#'
#' @return An object of class \code{"htest"}.
#'
#' @details
<<<<<<< HEAD
#' Using \code{conditioning="both"} corresponds to selecting \code{simulate.p.value=TRUE} in \code{\link{chisq.test}}. However, conditioning on both row and column marginals appears to be rarely justified in real data. Instead \code{conditioning="total"} is the correct choice for testing independence. Similarly, \code{conditioning="row"} is recommended when the row marginals e.g. are fixed by experimental design.
=======
#' Using \code{conditioning="both"} corresponds to selecting \code{simulate.p.value=TRUE} in \code{\link{chisq.test}}. However, conditioning on both row and column marginals appear to be rarely justified in real data. Instead \code{conditioning="total"} is the correct choice for testing independence. Similarly, \code{conditioning="row"} is recommended when the row marginals e.g. are fixed by experimental design.
>>>>>>> 0c9bab0686cb3f7b4ca3129c4adb4aee148c9624
#' The option \code{x0} has no effect when conditioning on both row and column marginals.
#'
#' @note The code has not been optimized for speed, and might be slow.
#'
#' @author Bo Markussen
#'
<<<<<<< HEAD
#' @seealso \code{\link{chisq.test}}
#'
#' @examples
=======
#' @seealso
#' \code{\link{chisq.test}}
#'
#' @example
>>>>>>> 0c9bab0686cb3f7b4ca3129c4adb4aee148c9624
#' # The Avadex dataset
#' Xobs <- matrix(c(2,3,6,40),2,2)
#' rownames(Xobs) <- c("Avadex +","Avadex -")
#' colnames(Xobs) <- c("Tumor +","Tumor -")
#'
#' # In this example only the rows appear to be fixed by experimental design.
#' # As is seen below, conditioning also on the columns is misleading conservative.
#' chisq.test.simulate(Xobs,"both")
#' chisq.test.simulate(Xobs,"row")
#' chisq.test.simulate(Xobs,"total")
#'
#' # Conditioning both on row and column marginals is simular to chisq.test().
#' chisq.test(Xobs,simulate.p.value=TRUE)
<<<<<<< HEAD
#' 
#' @export
=======


>>>>>>> 0c9bab0686cb3f7b4ca3129c4adb4aee148c9624
chisq.test.simulate <- function(x,conditioning="total",x0=NULL,B=1e4) {
  # if the null isn't stated, then take null from data
  if (is.null(x0)) {
    x0 <- x
  } else {
    if (any(dim(x)!=dim(x0))) stop("Dimensions of observed data and specification of null must be the same")
  }

  # set-up test statistic
  X2 <- function(y) {
    expected <- outer(rowSums(y),colSums(y))/sum(y)
    return(sum(((y-expected)^2)/(expected+(expected==0))))
  }

  # set-up simulation methods under the null hypothesis
  simulate.given.both.marginals <- function() {
    # remark: use marginals from observed data
    r2dtable(B,rowSums(x),colSums(x))
  }
  simulate.given.row.marginals <- function() {
    tmp <- matrix(t(sapply(rowSums(x),function(y){rmultinom(B,y,colSums(x0))})),nrow(x)*ncol(x),B)
    return(lapply(split(tmp,col(tmp)),function(y){matrix(y,nrow(x),ncol(x))}))
  }
  simulate.given.col.marginals <- function() {
    tmp <- matrix(t(sapply(colSums(x),function(y){rmultinom(B,y,rowSums(x0))})),ncol(x)*nrow(x),B)
    return(lapply(split(tmp,col(tmp)),function(y){t(matrix(y,ncol(x),nrow(x)))}))
  }
  simulate.given.total.sum <- function() {
    tmp <- rmultinom(B,sum(x),outer(rowSums(x0),colSums(x0)))
    return(lapply(split(tmp,col(tmp)),function(y){matrix(y,nrow(x),ncol(x))}))
  }

  # Expected counts given the observed marginals
  E <- outer(rowSums(x),colSums(x))/sum(x)

  # Observed test statistic
  X2.obs <- sum(((x-E)^2)/(E+(E==0)))
  names(X2.obs) <- "X-squared"

  # use selected simulation method
  method <- pmatch(conditioning,c("both","rows","cols","total"),4)
  if (method==1) {
    pval <- mean(unlist(lapply(simulate.given.both.marginals(),X2)) >= X2.obs)
    METHOD <- "Chi-squared test for given row and column marginals"
  }
  if (method==2) {
    pval <- mean(unlist(lapply(simulate.given.row.marginals(),X2)) >= X2.obs)
    METHOD <- "Chi-squared test for given row marginals"
  }
  if (method==3) {
    pval <- mean(unlist(lapply(simulate.given.col.marginals(),X2)) >= X2.obs)
    METHOD <- "Chi-squared test for given column marginals"
  }
  if (method==4) {
    pval <- mean(unlist(lapply(simulate.given.total.sum(),X2)) >= X2.obs)
    METHOD <- "Chi-squared test for given total sum"
  }

  # Standard error on p-value
  SE <- sqrt(pval*(1-pval)/B)
  names(SE) <- "Standard error (p-value)"

  # return result
  structure(list(statistic = X2.obs,
                 p.value = pval,
                 estimate = SE,
                 method = paste(METHOD,"(based on",B,"replicates)"),
                 data.name = deparse(substitute(x)),
                 observed = x,
                 expected = E,
                 residuals = (x - E)/sqrt(E+(E==0))), class = "htest")
}
