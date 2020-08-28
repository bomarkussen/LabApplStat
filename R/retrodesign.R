#' @title Post hoc power computation
#' 
#' @description 
#' `retrodesign` performs post hoc power analysis as proposed by Gelman & Carlin.
#' 
#' @keywords htest
#' 
#' @param A numeric vector with effect size.
#' @param SE numeric with standard error (default=1).
#' @param sig.level significance level (default=0.05).
#' @param df degrees of freedom (default=infinity).
#' @param B number of simulation replications used for typeM computation (default=10000).
#' 
#' @return A data frame with variables
#' \describe{
#'   \item{\code{power}}{power of test}
#'   \item{\code{typeS}}{risk of Type S error}
#'   \item{\code{exaggeration}}{average exaggeration ratio}
#' }
#' 
#' @references Gelman & Carlin: "Beyond Power Calculations: Assessing Type S (Sign) and Type M (Magnitude) Errors", Perspectives on Psychological Science, 1-11, 2014.
#'
#' @note R code adapted from Gelman & Carlin (2014).
#' 
#' @author Bo Markussen  
#' 
#' @examples 
#' # Basic design
#' mydata <- retrodesign(c(seq(0,1,0.01),seq(1,10,0.1),100),1)
#' 
#' # Type S error
#' plot(typeS~power,xlab="Power",ylab="Type S error",type="l",data=mydata)
#' 
#' # Type M error
#' plot(exaggeration~power,xlab="Power",ylab="Exaggeration ratio",ylim=c(0,12),type="l",
#' data=mydata)
#' abline(1,0,lty=2)
#' 
#' @export 
retrodesign <- function(A, SE=1, sig.level=0.05, df=Inf, B=10000){
  # A         = effect size (may be a vector)
  # SE        = standard error (default=1)
  # sig.level = significance level (default=5%)
  # df        = degrees of freedom (default=infinity)
  # B         = number of simulation replications used for typeM computation (default=10000)
  if (any(A<0)) stop("effect size must be non-negative")
  z <- qt(1-sig.level/2, df)
  p.hi <- 1 - pt(z-A/SE, df)
  p.lo <- pt(-z-A/SE, df)
  power <- p.hi + p.lo
  typeS <- p.lo/power
  mysims <- rt(B,df)
  typeM <- vapply(A,function(a) {estimate <- a+SE*mysims; mean(abs(estimate[abs(estimate)>SE*z]))/a},FUN.VALUE=1)
  return(data.frame(power=power, typeS=typeS, exaggeration=typeM))
}

