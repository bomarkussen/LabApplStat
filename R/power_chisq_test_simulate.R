#' Simulate power of Chi-squared tests with conditioning
#' 
#' @description 
#' `power.chisq.test.simulate` simulates power for tests for 2-way contingency tables based on the Pearson Chi-squared test statistics by simulation under 4 different conditioning scenarios.
#' 
#' @keywords htest
#' 
#' @param x matrix specifying the alternative distribution of the contingency table.
#' @param conditioning character string specifying the simulation scenario. Defaults to \code{"total"}. Other possible scenarios are \code{"row"}, \code{"col"}, and \code{"both"}.
#' @param x0 matrix specifying the null distribution. Defaults to \code{NULL}, in which case the null is estimated from the alternative \code{x}.
#' @param sig.level significance level used in test. Defaults to 0.05.
#' @param B integer specifying the number of replicates used in the Monte Carlo test. Defaults to 10000.
#' 
#' @return An object of class \code{"power.htest"}.
#' 
#' @details 
#' Using \code{conditioning="both"} corresponds to selecting \code{simulate.p.value=TRUE} in \code{\link{chisq.test}}. However, conditioning on both row and column marginals appears to be rarely justified in real data. Instead \code{conditioning="total"} is the correct choice for testing independence. Similarly, \code{conditioning="row"} is recommended when the row marginals e.g. are fixed by experimental design.
#' Both the alternative and the null are simulated under the parametric scenario estimated from the data matrix \code{x}. This possibly induces a discrepancy with \code{\link{chisq.test.simulate}}, where the null also is simulated from the specific data instance. Thus, the problem is that the null distribution depends on the model parameters.
#' 
#' @note The code has not been optimized for speed, and might be slow.
#' 
#' @author Bo Markussen
#' 
#' @seealso \code{\link{chisq.test.simulate}}
#' 
#' @examples 
#' # The Avadex dataset
#' Xobs <- matrix(c(2,3,6,40),2,2)
#' rownames(Xobs) <- c("Avadex +","Avadex -")
#' colnames(Xobs) <- c("Tumor +","Tumor -")
#' 
#' # In this example only the rows appear to be fixed by experimental design.
#' power.chisq.test.simulate(Xobs,"row")
#' power.chisq.test.simulate(Xobs,"total")
#' 
#' @export
power.chisq.test.simulate <- function(x,conditioning="total",x0=NULL,sig.level=0.05,B=1e4) {
  # if the null isn't stated, then take null from the alternative
  if (is.null(x0)) {
    x0 <- x
  } else {
    if (any(dim(x)!=dim(x0))) stop("Dimensions of alternative and of null must be the same")
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
  
  # set-up simulation methods under the alternative
  alternative.given.row.marginals <- function() {
    tmp <- matrix(t(apply(x,1,function(y){rmultinom(B,sum(y),y)})),nrow(x)*ncol(x),B)
    return(lapply(split(tmp,col(tmp)),function(y){matrix(y,nrow(x),ncol(x))}))
  }
  alternative.given.col.marginals <- function() {
    tmp <- matrix(t(apply(x,1,function(y){rmultinom(B,sum(y),y)})),ncol(x)*nrow(x),B)
    return(lapply(split(tmp,col(tmp)),function(y){t(matrix(y,ncol(x),nrow(x)))}))
  }
  alternative.given.total.sum <- function() {
    tmp <- rmultinom(B,sum(x),as.vector(x))
    return(lapply(split(tmp,col(tmp)),function(y){matrix(y,nrow(x),ncol(x))}))
  }
  
  # use selected simulation method
  method <- pmatch(conditioning,c("both","rows","cols","total"),4)
  if (method==1) {
    stop("Power computation not yet implemented for conditioning on both marginals")
    HA <- simulate.given.both.marginals()
    crit <- quantile(unlist(lapply(simulate.given.both.marginals(),X2)),1-sig.level,type=1)
    pwr <- mean(unlist(lapply(HA,X2)) > crit)
    METHOD <- "Chi-squared test for given row and column marginals"
  }
  if (method==2) {
    HA <- alternative.given.row.marginals()
    crit <- quantile(unlist(lapply(simulate.given.row.marginals(),X2)),1-sig.level,type=1)
    pwr <- mean(unlist(lapply(HA,X2)) > crit)
    METHOD <- "Chi-squared test for given row marginals"
  }
  if (method==3) {
    HA <- alternative.given.col.marginals()
    crit <- quantile(unlist(lapply(simulate.given.col.marginals(),X2)),1-sig.level,type=1)
    pwr <- mean(unlist(lapply(HA,X2)) > crit)
    METHOD <- "Chi-squared test for given column marginals"
  }
  if (method==4) {
    HA <- alternative.given.total.sum()
    crit <- quantile(unlist(lapply(simulate.given.total.sum(),X2)),1-sig.level,type=1)
    pwr <- mean(unlist(lapply(HA,X2)) > crit)
    METHOD <- "Chi-squared test for given total sum"
  }
  
  # Standard error on power
  SE <- sqrt(pwr*(1-pwr)/B)
  names(SE) <- "Standard error (power)"
  
  # return result
  structure(list(N = sum(x),
                 sig.level = sig.level,
                 power = pwr,
                 SE.power = SE,
                 data.name = deparse(substitute(x)),
                 method = paste(METHOD,"(based on",B,"replicates)")),
            class = "power.htest")
}
