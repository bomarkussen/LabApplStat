#' @title Make emmeans object for an expected dose
#' 
#' @description 
#' Solves linear equations in continuous explanatory variables in order to find the expected dose. A typical application could be to find LD50, i.e. the lethal dose killing 50 percent of the population, from a probit analysis fitted by \code{\link{glm}}. The associated variance-covariance matrix is found using the Delta method.
#' 
#' @param object An object that can be given to \code{\link[emmeans]{emmeans}}. Typically a model fitted by \code{\link{glm}}.
#' @param specs As for \code{\link[emmeans]{emmeans}}. Typically as one-sided \code{\link{formula}}.
#' @param left A list specifying the left end point of the linear span of continuous variables in which to measure the ED values.
#' @param right A list specifying the right end point of the linear span of continuous variables in which to measure the ED values.
#' @param tran Possible transformation of the scale of the ED values. If given then backtransformation can be done using the technology of the \code{\link[emmeans]{emmeans}}. The default value \code{tran=NULL} corresponds to no transformation.
#' @param p Numeric vector given the targeted predictions. Typically probabilities, where the default value \code{p=0.5} corresponds to ED50.
#' @param p.name The name of the variable containing \code{p}. If \code{p} contains more than one value, then this will also appear in \code{@misc$by.vars} in the \code{emmGrid} object. Defaults to \code{p.name="probability"}.
#' 
#' @author Bo Markussen
#' 
#' @export 
emmeans_ED <- function(object,specs,left,right,tran=NULL,p=0.5,p.name="probability") {
  # Find reference grids
  suppressMessages(em0 <- emmeans(object,specs,at=left))
  suppressMessages(em1 <- emmeans(object,specs,at=right))
  # Find intercept
  if (is.null(em0@misc$tran))      intercept <- p
  if (is.list(em0@misc$tran))      intercept <- em0@misc$tran$linkfun(p)
  if (is.character(em0@misc$tran)) intercept <- stats::make.link(em0@misc$tran)$linkfun(p)
  # Find variables and levels
  pri.vars <- em0@misc$pri.vars
  by.vars  <- em0@misc$by.vars
  xlev     <- em0@levels
  # Find grid
  if (is.element(".wgt.",c(pri.vars,by.vars))) {
    grid <- em0@grid
  } else {
    grid <- em0@grid[,!is.element(names(em0@grid),".wgt.")]
  }
  # Extend if length(p)>1
  if (length(p)>1) {
    if (is.element(p.name,c(pri.vars,by.vars))) stop("p.name must not appear as a variable in object")
    tmp     <- data.frame(p=factor(p)); names(tmp) <- p.name
    grid    <- merge(grid,tmp)
    by.vars <- c(by.vars,p.name)
    tmp     <- list(p=as.character(p)); names(tmp) <- p.name
    xlev    <- c(xlev,tmp)
  }
  # Find linear transformations
  A <- em0@linfct
  B <- em1@linfct - em0@linfct
  # Find indices of estimable LD's
  nbasis0 <- em0@nbasis; if (all(is.na(nbasis0))) nbasis0 <- rep(0,length(em0@bhat)) 
  nbasis1 <- em1@nbasis; if (all(is.na(nbasis1))) nbasis1 <- rep(0,length(em1@bhat)) 
  ii <- apply(cbind(em0@linfct %*% nbasis0,em1@linfct %*% nbasis1),1,function(x) {all(x==0)})
  # Find estimates
  hat.theta <- em0@bhat
  jj    <- !is.na(hat.theta)
  alpha <- c(A[,jj,drop=FALSE]%*%hat.theta[jj])
  beta  <- c(B[,jj,drop=FALSE]%*%hat.theta[jj])
  if (length(p)>1) {
    f <- c(-alpha/beta,1/beta)
    f[rep(!ii,length(p))] <- NA
  } else {
    f <- (intercept-alpha)/beta
    f[!ii] <- NA
  }
  # Find variance matrix
  if (length(p)>1) {
    H <- rbind(cbind(diag(-1/beta[ii],sum(ii)),
                     diag(alpha[ii]/(beta[ii]^2),sum(ii))),
               cbind(diag(0,sum(ii)),
                     diag(-1/(beta[ii]^2),sum(ii))))%*%rbind(A[ii,jj],B[ii,jj])
  } else {
    H <- cbind(diag(-1/beta[ii],sum(ii)),diag((alpha[ii]-intercept)/(beta[ii]^2),sum(ii)))%*%rbind(A[ii,jj],B[ii,jj])
  }
  # Make emmGrid object by altering em0
  em <- em0
  em@model.info <- list(call = match.call(), xlev = xlev)
  if (length(p)>1) em@roles$predictors <- c(em@roles$predictors,p.name)
  em@grid   <- grid
  em@levels <- xlev
  em@bhat   <- f
  em@V      <- H%*%(em@V)%*%t(H)
  if (length(p)>1) {
    em@linfct <- matrix(aperm(outer(diag(length(alpha)),
                                    matrix(c(rep(1,length(p)),intercept),length(p),2)),
                              c(1,3,2,4)),
                        length(alpha)*length(p),length(alpha)*2)
  } else {
    em@linfct <- diag(length(alpha))
  }
  if (sum(ii)==0) {
    em@nbasis <- matrix(NA,1,1)
  } else {
    em@nbasis <- matrix(rep(as.numeric(!ii),min(2,length(p))),
                        length(ii)*min(2,length(p)),1)
  }
  em@dffun <- function(k,dfargs) {Inf}
  em@dfargs <- list()
  em@misc <- list(estName = "estimate", estType = "prediction", 
                  infer = c(TRUE, FALSE), level = 0.95, adjust = "none", 
                  famSize = nrow(em@linfct), avgd.over = character(0), 
                  pri.vars = pri.vars, by.vars = by.vars,
                  methDesc = "emmobj")
  if (!is.null(tran)) em <- stats::update(em,tran=tran)
  # return result
  return(em)
}
