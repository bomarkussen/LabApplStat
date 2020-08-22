#' Design diagram for a linear model
#' 
#' @description 
#' `DD` computes the Design Diagram for a linear model.
#' 
#' @keywords manip design
#' 
#' @param fixed formula with fixed effects. A response may the specified, but this optional.
#' @param random formula with random effects. Defaults to \code{NULL} meaning that there are no other random effects than the residual, which is added to all designs.
#' @param data data frame with the explanatory variables and the response (if specified).
#' @param eps threshold for deeming singular values to be "zero". Defaults to 1e-12.
#' 
#' @return An object of class \code{\link{designDiagram-class}}
#' 
#' @author Bo Markussen
#' 
#' @seealso \code{\link{minimum}}, \code{\link{plot.designDiagram}}
#' 
#' @examples
#' # 3-way ANOVA
#' x <- factor(rep(rep(1:4,times=4),times=4))
#' y <- factor(rep(rep(1:4,times=4),each=4))
#' z <- factor(rep(rep(1:4,each=4),each=4))
#' myDD <- DD(~x*y*z,data=data.frame(x=x,y=y,z=z))
#' summary(myDD)
#' 
#' #Making the factor diagram closed under minima
#' mydata <- data.frame(age=rep(c("boy","girl","adult","adult"),4),
#'                      gender=rep(c("child","child","man","woman"),4))
#' myDD <- DD(~0+age+gender,data=mydata)
#' plot(myDD)
#' 
#' # Example of collinearity
#' mydata <- data.frame(age=rnorm(102),edu=rnorm(102),sex=factor(rep(c(1,2),51)))
#' mydata <- transform(mydata,exper=age-edu+0.1*rnorm(102))
#' mydata <- transform(mydata,wage=2*edu+2*exper+rnorm(102))
#' summary(myDD <- DD(wage~sex*(age+exper+edu),data=mydata))
#' 
#' # growth of rats
#' antibiotica <- factor(rep(c(0,40),each=6))
#' vitamin <- factor(rep(rep(c(0,5),each=3),2))
#' growth <- c(1.30,1.19,1.08,1.26,1.21,1.19,1.05,1.00,1.05,1.52,1.56,1.55)
#' mydata <- data.frame(antibiotica=antibiotica,vitamin=vitamin,growth=growth)
#' plot(DD(growth~antibiotica*vitamin,data=mydata),"MSS")
#' 
#' \dontrun{
#'   # ANCOVA: Non-orthogonal design
#'   library(isdals)
#'   data(birthweight)
#'   plot(DD(weight~sex*I(age-42),data=birthweight),"MSS")
#'   plot(DD(weight~I(age-42)+sex:I(age-42)+sex,data=birthweight),"MSS")
#' }
#' 
#' @export
DD <- function(fixed,random=NULL,data,eps=1e-12) {
  # sanity check
  if (class(fixed)!="formula") stop("fixed-argument must be a formula")
  if ((!is.null(random)) && (class(random)!="formula")) stop("random-argument must be a formula")
  if (!is.data.frame(data)) stop("data-argument must be a data frame")
  
  # -------------------------
  # Initialize
  # -------------------------
  
  # find terms in the design and place square brackets around random terms
  myterms <- attr(terms(fixed,keep.order=TRUE),"term.labels")
  if (attr(terms(fixed),"intercept")==1) myterms <- c("1",myterms)
  if (is.null(random)) {myterms.random <- NULL} else {
    myterms.random <- setdiff(attr(terms(random,keep.order=TRUE),"term.labels"),myterms)
    if (length(myterms.random)>0) {
      myterms.random <- paste("[",myterms.random,"]",sep="")
      myterms <- c(myterms,myterms.random)
    }
  }
  M <- length(myterms)

  # find number of complete observations
  N <- nrow(model.frame(formula(paste0("~",paste(sub("[","",sub("]","",myterms,fixed=TRUE),fixed=TRUE),collapse="+"))),data))
  
  # extract response variable if it is present
  y <- model.response(model.frame(fixed,data))
  
  # make design matrices for each of the terms
  mydesigns <- vector("list",M)
  Nparm <- rep(0,M)
  for (i in 1:M) {
    # removed square brackets if present
    tmp <- sub("[","",sub("]","",myterms[i],fixed=TRUE),fixed=TRUE)
    # make design matrix
    A <- model.matrix(as.formula(paste0("~0+",tmp)),data=data)
    # remove non-needed attributed
    attr(A,"assign") <- NULL
    attr(A,"contrasts") <- NULL
    # remove zero columns
    A <- A[,!apply(A,2,function(x){all(x==0)}),drop=FALSE]
    # reduce to full rank (relative to eps) basis
    tmp <- svd(A)
    mydesigns[[i]] <- tmp$u[,tmp$d>eps,drop=FALSE]
    # extract number of parameters
    Nparm[i] <- ncol(mydesigns[[i]])
  }
  
  # Initialize matrix of relations
  relations <- matrix(NA,M,M)
  diag(relations) <- "="

  # --------------------------------------------------------------------------
  # Complete relations:
  #   1) fills in matrix of relations
  #   2) remove duplicate variables
  #   3) extend with missing minima
  #   Uses and modifies the variables: M, myterms, mydesigns, Nparm, relations
  # --------------------------------------------------------------------------
  
  while (any(is.na(relations))) {
    # Find row and column of non-decided relation
    j <- min(which(is.na(relations)))
    i <- 1+((j-1) %% M)
    j <- 1+((j-1) %/% M)
    # Find common subspace, i.e. the null space for the union 
    tmp <- svd(cbind(mydesigns[[i]],mydesigns[[j]]))
    NullDim <- min(c(sum(tmp$d<=eps)+max(0,dim(tmp$v)[1]-dim(tmp$u)[1]),Nparm[i],Nparm[j]))
    # Look for ordering
    if ((NullDim==0) | (NullDim==Nparm[i]) | (NullDim==Nparm[j])) {
      # Are designs linearly independent?
      if (NullDim==0) {
        relations[i,j] <- relations[j,i] <- "0"
      } else {
        # Is there an order between the designs?
        if (Nparm[i] < Nparm[j]) {relations[i,j] <- "<"; relations[j,i] <- ">"}
        if (Nparm[i] > Nparm[j]) {relations[i,j] <- ">"; relations[j,i] <- "<"}
        if (Nparm[i]==Nparm[j]) {
          if ((is.element(myterms[i],myterms.random)) & (!is.element(myterms[j],myterms.random))) {
            warning(paste(myterms[j],"and",myterms[i],"are (almost) identical. The term with random effect is removed from the design."))
            myterms.random <- setdiff(myterms.random,myterms[i])
            k <- i
          }
          if ((!is.element(myterms[i],myterms.random)) & (is.element(myterms[j],myterms.random))) {
            warning(paste(myterms[i],"and",myterms[j],"are (almost) identical. The term with random effect is removed from the design."))
            myterms.random <- setdiff(myterms.random,myterms[j])
            k <- j
          }
          if ((!is.element(myterms[i],myterms.random)) & (!is.element(myterms[j],myterms.random))) {
            warning(paste(myterms[i],"and",myterms[j],"are (almost) identical. The latter is removed from the design."))
            k <- j
          }
          M <- M-1
          myterms <- myterms[-k]
          mydesigns <- mydesigns[-k]
          Nparm <- Nparm[-k]
          relations <- relations[-k,-k,drop=FALSE]            
        }
      }
    } else {
      # Design are neither linearly independent nor ordered
      # Parametrize common space (ie. minimum) (using design with lowest index)
      if (i<j) {
        NullSpace <- mydesigns[[i]]%*%tmp$v[1:Nparm[i],tmp$d<=eps,drop=FALSE]
      } else {
        NullSpace <- mydesigns[[j]]%*%tmp$v[(Nparm[i]+1):(Nparm[i]+Nparm[j]),tmp$d<=eps,drop=FALSE]
      }
      tmp <- svd(NullSpace,nv=0)
      NullSpace <- tmp$u[,tmp$d>eps,drop=FALSE]
      # Look for minimum among the existing terms
      for (k in setdiff((1:M),c(i,j))) {
        if ((NullDim==Nparm[k]) && (NullDim==sum(svd(cbind(NullSpace,mydesigns[[k]]))$d>eps))) {
          relations[i,j] <- relations[j,i] <- myterms[k]
          break
        }
      }
      # Has minimum been found among the existing terms?
      if (is.na(relations[i,j])) {
        # Make name for new minimum
        tmp <-  paste0("min{",sub("[","",sub("]","",myterms[min(c(i,j))],fixed=TRUE),fixed=TRUE),
                       ",",sub("[","",sub("]","",myterms[max(c(i,j))],fixed=TRUE),fixed=TRUE),"}")
        # Minimum is random if and only if both variables are random
        if ((is.element(myterms[min(c(i,j))],myterms.random)) & (is.element(myterms[max(c(i,j))],myterms.random))) {
          tmp <- paste0("[",tmp,"]")
          myterms.random <- c(myterms.random,tmp)
        }
        # Add design of minimum.
        M <- M+1
        myterms   <- c(myterms,tmp)
        Nparm     <- c(Nparm,NullDim)
        mydesigns <- c(mydesigns,list(NullSpace))
        relations[i,j] <- relations[j,i] <- myterms[M]
        relations <- cbind(rbind(relations,rep(NA,M-1)),rep(NA,M))
        relations[M,M] <- "="
      }
    }
    # End while() loop
  }

  # -------------------------------------------------------
  # Remove underlying designs and find degrees of freedom
  # -------------------------------------------------------
  
  # Remove nested designs from the designs
  if (M>0) for (i in M:1) {
    tmp <- is.element(relations[i,],c(">","->"))
    if (any(tmp)) {
      B <- NULL
      for (j in which(tmp)) B <- cbind(B,mydesigns[[j]])
      tmp <- svd(B,nv=0)
      B <- tmp$u[,tmp$d>eps,drop=FALSE]
      tmp <- svd(mydesigns[[i]]-B%*%t(B)%*%mydesigns[[i]],nv=0)
      mydesigns[[i]] <- tmp$u[,tmp$d>eps,drop=FALSE]
    }
  }
  
  # Compute degrees of freedom
  mydf <- unlist(lapply(mydesigns,function(x){dim(x)[2]}))
  
  
  # -------------------------------------------------------
  # Investigate orthogonality by computing inner products
  # -------------------------------------------------------

  # Compute inner products
  inner <- matrix(NA,M,M)
  for (i in 1:M) for (j in 1:M) {
    inner[i,j] <- round(sum(c(t(mydesigns[[i]])%*%mydesigns[[j]])^2),floor(-log10(eps)))
  }
  
  # Issue warning for non-orthogonal designs
  if (any(inner[upper.tri(inner)]!=0)) warning("Design is non-orthogonal: Sum-of-Squares and p-values may depend on order of terms.")
  
  
  # ----------------------------------------------------
  # Find sequential ordering of the terms
  # Uses and modifies the variables: 
  #   myorder, M, myterms, mydesigns, Nparm, relations
  # ----------------------------------------------------

  # Initialize and find basis
  myorder <- rep(NA,M)

  # Loop through all variables
  for (k in 1:M) {
    # find next variable in a sequential ordering of the variables
    if (k==1) {
      # Note: which.min() selects the index of the first minimum, which hence 
      #       complies with the initial order as much as possible.
      i <- which.min(M*apply(relations==">",1,sum)-apply(relations=="<",1,sum))
    } else {
      i <- which.min(M*is.element(1:M,myorder[1:(k-1)])
                     +M*apply(relations[,-myorder[1:(k-1)],drop=FALSE]==">",1,sum)
                     -apply(relations[,-myorder[1:(k-1)],drop=FALSE]=="<",1,sum))
    }
    myorder[k] <- i
  }

  # Reorder variables
  myterms   <- myterms[myorder]
  Nparm     <- Nparm[myorder]
  relations <- relations[myorder,myorder]
  mydesigns <- mydesigns[myorder]
  

  # ----------------------------------
  # Compute summaries and statistics
  # ----------------------------------

  # Extend with the identity variable
  # TO DO: What is the identity variable is already included!?
  #        Then it will have df=0, but what does this imply??
  myterms.random <- c(myterms.random,"[I]")
  myterms   <- c(myterms,"[I]")
  Nparm     <- c(Nparm,N)
  mydf      <- c(mydf,Nparm[M+1]-sum(mydf))
  relations <- cbind(rbind(relations,matrix(">",1,M)),matrix(c(rep("<",M),"="),M+1,1))
  M <- M+1
  
  # Make ghost relations, where terms with df=0 are removed.
  # This is used to find testable hypothesis, i.e. collapses that dont remove non-trivial minima
  relations.ghost <- relations
  relations.ghost[mydf==0,] <- "df=0"
  relations.ghost[,mydf==0] <- "df=0"
  
  # Identify arrows: Should be run after computation of degrees of freedom
  for (i in 1:M) {
    below <- (relations[i,]==">")
    if (any(below)) {
      below <- below&(!apply(relations[below,,drop=FALSE],2,function(x){any(is.element(x,c(">","->")))}))
      relations[i,below] <- "->"
      relations[below,i] <- "<-"
    }
  }
  
  # Identify arrows in ghost relations
  for (i in 1:M) {
    below <- (relations.ghost[i,]==">")
    if (any(below)) {
      below <- below&(!apply(relations.ghost[below,,drop=FALSE],2,function(x){any(is.element(x,c(">","->")))}))
      relations.ghost[i,below] <- "->"
      relations.ghost[below,i] <- "<-"
    }
  }
  
  # Find projections onto basis functions and make Type-I F-tests
  pvalue <- matrix(NA,M,M)
  SS     <- matrix(NA,M,M)
  MSS    <- matrix(NA,M)
  if (!is.null(y)) {
    # find orthogonal basis
    if (M>1) for (i in 1:(M-1)) {
      for (j in i:(M-1)) {
        if (i==j) {
          A <- mydesigns[[j]]
        } else {
          tmp <- svd(do.call("cbind",mydesigns[i:(j-1)]),nv=0)
          B <- tmp$u[,tmp$d>eps,drop=FALSE]
          A <- mydesigns[[j]]-B%*%t(B)%*%mydesigns[[j]]
          tmp <- svd(A,nv=0)
          A <- tmp$u[,tmp$d>eps,drop=FALSE]
        }
        SS[i,j] <- sum(c(y%*%A)^2)
      }
    }
    # Residual SS
    SS[M,M] <- sum(y^2)
    SS[1:(M-1),M] <- SS[M,M]-rowSums(SS[1:(M-1),1:(M-1)],na.rm=TRUE)

    # Compute MSS
    MSS <- SS*matrix(rep(ifelse(mydf>0,1/mydf,0),each=M),M,M)

    # F-tests only for terms with positive degrees of freedom, which have
    # only one term nested within them.
    for (i in which((mydf>0) & apply(relations.ghost,1,function(x){sum(x=="<-")})==1)) {
      # Find candidates for denominator: random effects with positive degrees of freedom, which are "<" or "<-"
      j <- (is.element(myterms,myterms.random) & (mydf>0) & is.element(relations.ghost[i,],c("<","<-")))
      # Remove candidates that are nested within some of the other candidates
      if (sum(j)>0) j[j] <- j[j] & (!apply(relations.ghost[j,j,drop=FALSE],2,function(x){any(is.element(x,c("<","<-")))}))
      # Make F-test if only one(!) candidate is left
      if (sum(j)==1) {
        j <- which(j)
        # allocate p-value to an "<-" with df>0
        k <- which((mydf>0) & (relations[i,]=="<-"))
        if (length(k)==0) k <- which(relations[i,]=="<-")[1]
        pvalue[i,k] <- 1-pf(MSS[1,i]/MSS[1,j],mydf[i],mydf[j])
      }
    }
  }

  # return result
  names(myterms) <- names(Nparm) <- names(mydf) <- 
    colnames(SS) <- colnames(MSS) <- rownames(relations) <- colnames(relations) <- 
    rownames(pvalue) <- colnames(pvalue) <- myterms
  rownames(SS) <- rownames(MSS) <- c("-",myterms[-M])
  rownames(inner) <- colnames(inner) <- names(mydesigns) <- myterms[-M]
  res <- list(terms=myterms,random.terms=myterms.random,Nparm=Nparm,df=mydf,
              SS=SS,MSS=MSS,relations=relations,pvalue=pvalue,
              inner=inner,response=!is.null(y))
  class(res) <- "designDiagram"
  return(res)
}
