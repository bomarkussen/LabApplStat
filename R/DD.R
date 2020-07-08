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
#' @param threshold threshold for removing (approximative) collinearities in the design. Defaults to 0.1.
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
#' #plot(myDD)
#' 
#' # Example of collinearity
#' mydata <- data.frame(age=rnorm(102),edu=rnorm(102),sex=factor(rep(c(1,2),51)))
#' mydata <- transform(mydata,exper=age-edu)
#' summary(myDD <- DD(~sex*(age+exper+edu),data=mydata))
#' 
#' # growth of rats
#' antibiotica <- factor(rep(c(0,40),each=6))
#' vitamin <- factor(rep(rep(c(0,5),each=3),2))
#' growth <- c(1.30,1.19,1.08,1.26,1.21,1.19,1.05,1.00,1.05,1.52,1.56,1.55)
#' mydata <- data.frame(antibiotica=antibiotica,vitamin=vitamin,growth=growth)
#' #plot(DD(growth~antibiotica*vitamin,data=mydata),"MSS")
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
DD <- function(fixed,random=NULL,data,threshold=0.1,eps=1e-12) {
  # sanity check
  if (class(fixed)!="formula") stop("fixed-argument must be a formula")
  if ((!is.null(random)) && (class(random)!="formula")) stop("random-argument must be a formula")
  if (!is.data.frame(data)) stop("data-argument must be a data frame")
  
  # -------------------------
  # Initialize
  # -------------------------
  
  # find terms in the design and place square brakets around random terms
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

  # -------------------------------------------------------------------------
  # Find basis:
  #   1) find sequential ordering of the variables
  #   2) make associated basis for independent designs "from right to left"
  #   3) remove approximate collinearities (relative to threshold-argument)
  #   Uses and modifies the variables: myorder, mybasis
  #   Also uses the variables: M, myterms, mydesigns, Nparm, relations
  # -------------------------------------------------------------------------

  # Initialize and find basis
  mybasis <- vector("list",M)
  myorder <- rep(NA,M)

  # Initialize basis for lower order variables
  B <- matrix(0,N,0)
  
  # Loop through all variables
  for (k in 1:M) {
    # find next variable in a sequential ordering of the variables
    if (k==1) {
      # Note: which.min() selects the index of the first minimum, which hence complies with the 
      #       initial order as much as possible.
      i <- which.min(M*apply(relations==">",1,sum)-apply(relations=="<",1,sum))
    } else {
      i <- which.min(M*is.element(1:M,myorder[1:(k-1)])
                     +M*apply(relations[,-myorder[1:(k-1)],drop=FALSE]==">",1,sum)
                     -apply(relations[,-myorder[1:(k-1)],drop=FALSE]=="<",1,sum))
    }
    # find orthogonal basis
    A <- mydesigns[[i]]-B%*%t(B)%*%mydesigns[[i]]
    tmp <- svd(A)
    mybasis[[i]] <- tmp$u[,tmp$d>threshold,drop=FALSE]
    # update myorder
    myorder[k] <- i
    # updata basis of lower order variables
    B <- cbind(B,mybasis[[i]])
  }

  # ----------------------------------
  # Compute summaries and statistics
  # ----------------------------------

  # Reorder variables
  myterms   <- myterms[myorder]
  Nparm     <- Nparm[myorder]
  relations <- relations[myorder,myorder]
  mybasis   <- mybasis[myorder]
  mydesigns <- mydesigns[myorder]
  
  # compute degrees of freedom
  mydf <- unlist(lapply(mybasis,function(x){dim(x)[2]}))

  # compute number of collinearities
  mycol <- rep(NA,M)
  for (i in 1:M) {
    mycol[i] <- Nparm[i]-sum(mydf[relations[i,]==">"])-mydf[i]
  }
  
  # Extend with the identity variable
  myterms.random <- c(myterms.random,"[I]")
  myterms   <- c(myterms,"[I]")
  Nparm     <- c(Nparm,N)
  mydf      <- c(mydf,Nparm[M+1]-sum(mydf))
  mycol     <- c(mycol,0)
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
  
  # Find projections onto basis functions and make F-tests
  mycoef <- vector("list",M-1)
  pvalue <- matrix(NA,M,M)
  SS     <- rep(NA,M)
  MSS    <- rep(NA,M)
  if (!is.null(y)) {
    # projections onto basis for non-residual terms
    if (M>1) for (i in 1:(M-1)) mycoef[[i]] <- c(y%*%mybasis[[i]])
    # projections onto basis for residual term: REMOVED, but could be reintroduced for purpose of model validation
    #mycoef[[M]] <- numeric(0)
    #if (mydf[M]>0) {
    #  # find basis for residual term
    #  A <- mybasis[[1]]
    #  if (M>2) for (i in 2:(M-1)) A <- cbind(A,mybasis[[i]])
    #  tmp <- svd(A)
    #  A <- tmp$u[,tmp$d>eps]
    #  # find projections
    #  mycoef[[M]] <- c(y%*%svd(cbind(diag(1,nrow=length(y)),A))$u[,dim(A)[2]+1:mydf[M]])
    #}

    # Compute SS and MSS
    SS[M] <- sum(y^2)
    if (M>1) {
      SS[1:(M-1)]  <- unlist(lapply(mycoef,function(x){sum(x^2)}))
      SS[M] <- SS[M] - sum(SS[1:(M-1)])
    }
    MSS[1:M]     <- SS/mydf
    MSS[mydf==0] <- 0

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
        pvalue[i,k] <- 1-pf(MSS[i]/MSS[j],mydf[i],mydf[j])
      }
    }
  }

  # Investigate orthogonality:
  # 1. Remove underlying nested designs from designs. Done top-down!
  if (M>1) for (i in (M-1):1) {
    if (mydf[i]==0) {
      mydesigns[[i]] <- matrix(0,nrow(data),0)
    } else {
      tmp <- is.element(relations[i,],c(">","->"))
      if (any(tmp)) {
        B <- NULL
        for (j in which(tmp)) B <- cbind(B,mydesigns[[j]])
        tmp <- svd(B)
        B <- tmp$u[,tmp$d>eps]
        mydesigns[[i]] <- svd(mydesigns[[i]]-B%*%t(B)%*%mydesigns[[i]],nu=mydf[i],nv=0)$u
      }
    }
  }
  # 2. Compute inner products
  inner <- matrix(NA,M-1,M-1)
  for (i in 1:(M-1)) for (j in 1:(M-1)) {
    inner[i,j] <- round(sum(c(t(mydesigns[[i]])%*%mydesigns[[j]])^2),6)
  }
  # 3. issue warning for non-orthogonal designs
  if (any(inner[upper.tri(inner)]!=0)) warning("Design is non-orthogonal: Sum-of-Squares and p-values may depend on order of terms.")
  
  # return result
  names(myterms) <- names(Nparm) <- names(mydf) <- names(mycol) <- names(SS) <- names(MSS) <-
    rownames(relations) <- colnames(relations) <- rownames(pvalue) <- colnames(pvalue) <- myterms
  rownames(inner) <- colnames(inner) <- names(mycoef) <- names(mybasis) <- names(mydesigns) <- myterms[-M]
  res <- list(terms=myterms,random.terms=myterms.random,Nparm=Nparm,df=mydf,collinearities=mycol,
              SS=SS,MSS=MSS,relations=relations,pvalue=pvalue,
              inner=inner,response=!is.null(y))
  class(res) <- "designDiagram"
  return(res)
}
