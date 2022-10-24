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
#' @param keep formula which effects that will not be removed in the collinarity analysis. Defaults to \code{~1} meaning that the intercept will be kept if it is present.
#' @param center boolean deciding whether to centralize numerical predictors when an intercept is present. Defaults to \code{FALSE}.
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
#' myDD <- DD(growth~antibiotica*vitamin,data=mydata)
#' plot(myDD,"MSS")
#' plot(myDD,"I2")
#' 
#' \donttest{
#'   # ANCOVA: Non-orthogonal design
#'   library(isdals)
#'   data(birthweight)
#'   plot(DD(weight~sex*I(age-42),data=birthweight),"MSS")
#'   plot(DD(weight~I(age-42)+sex:I(age-42)+sex,data=birthweight),"MSS")
#' }
#' 
#' @export
DD <- function(fixed,random=NULL,data,keep=~1,center=FALSE,eps=1e-12) {
  # sanity check
  if (!inherits(fixed,"formula")) stop("fixed-argument must be a formula")
  if ((!is.null(random)) && (!inherits(random,"formula"))) stop("random-argument must be a formula")
  if (!is.data.frame(data)) stop("data-argument must be a data frame")
  if (!inherits(keep,"formula"))  stop("keep-argument must be a formula")
  
  # -------------------------
  # Initialize
  # -------------------------

  # find model frame and number of rows
  N.all <- nrow(data)
  if (is.null(random)) {
    data <- model.frame(fixed,data)
  } else {
    data <- model.frame(update(fixed,as.formula(paste(as.character(random),collapse=".+"))),data)
  }
  N <- nrow(data)
  if (N!=N.all) warning(paste("Removed",N.all-N,"rows with missing values in response or explanatory variables"))

  # TO DO: add option "keep.order=TRUE" to terms() call below?
  
  # find terms in the design and place square brackets around random terms
  myterms <- attr(terms(fixed,data=data),"term.labels")
  if (attr(terms(fixed,data=data),"intercept")==1) {
    myterms <- c("1",myterms)
  } else {
    # if no intercept, then centralization is switched off!
    center <- FALSE
  }
  if (is.null(random)) {myterms.random <- NULL} else {
    # terms in random-option will be treated as random
    # but these may also be given in fixed-option to specify the order
    myterms.random <- attr(terms(random,data=data),"term.labels")
    if (length(myterms.random)>0) {
      myterms <- c(myterms,setdiff(myterms.random,myterms))
      i <- is.element(myterms,myterms.random)
      myterms[i] <- paste("[",myterms[i],"]",sep="")
      myterms.random <- paste("[",myterms.random,"]",sep="")
    }
  }
  M <- length(myterms)

  # stop if no terms
  # TO DO: Should we allow for M=0 below? If so, then updates are required
  if (M==0) stop("There should be at least one term (including the intercept)")

  # for each term
  #   1) check if it is a categorical factor
  #   2) find design matrices in the original parametrization
  #   4) make orthonormal design matrices
  myisfactor <- rep(TRUE,M)
  mydesigns  <- vector("list",M)
  mybasis    <- vector("list",M)
  Nparm      <- rep(0,M)
  for (i in 1:M) {
    # removed square brackets if present
    tmp <- sub("[","",sub("]","",myterms[i],fixed=TRUE),fixed=TRUE)
    # make design matrix
    A <- model.matrix(as.formula(paste0("~0+",tmp)),data=data)
    # check if is is a categorical factor
    myisfactor[i] <- (length(setdiff(attr(terms(as.formula(paste0("~0+",tmp))),"term.labels"),
                                     names(attr(A,"contrasts"))))==0)
    # remove zero columns
    A <- A[,!apply(A,2,function(x){all(x==0)}),drop=FALSE]
    # centralize numerical variables?
    if (center && (tmp!="1") && (!myisfactor[i])) {
      A <- A - mean(c(A),na.rm=TRUE)
    }
    # remove non-needed attributed
    attr(A,"assign") <- NULL
    attr(A,"contrasts") <- NULL
    # reduce to full rank (relative to eps) basis
    if (ncol(A)>0) tmp <- svd(A)
    if ((ncol(A)>0) && (sum(tmp$d>eps)>0)) {
      mydesigns[[i]] <- tmp$u[,tmp$d>eps,drop=FALSE]%*%diag(tmp$d[tmp$d>eps],nrow=sum(tmp$d>eps))%*%t(tmp$v[,tmp$d>eps,drop=FALSE])
      mybasis[[i]]   <- tmp$u[,tmp$d>eps,drop=FALSE]
    } else {
      N <- nrow(data)  # Hack: Check this!
      mydesigns[[i]] <- matrix(0,N,0)
      mybasis[[i]] <- matrix(0,N,0)
    }
    # extract number of parameters
    Nparm[i] <- ncol(mybasis[[i]])
  }

  
  # extract response variable if it is present
  y <- model.response(data)
  
  
  # --------------------------------------------------------------------------
  # Complete relations:
  #   1) fills in matrix of relations
  #   2) remove duplicate variables
  #   3) extend with missing minima
  #   Uses and modifies the variables: 
  #     M, myterms, mydesigns, mybasis, Nparm, relations
  # --------------------------------------------------------------------------
  
  # Initialize matrix of relations
  relations <- matrix(NA,M,M)
  diag(relations) <- "="

  # Complete relations
  while (any(is.na(relations))) {
    # Find row and column of non-decided relation
    j <- min(which(is.na(relations)))
    i <- 1+((j-1) %% M)
    j <- 1+((j-1) %/% M)
    # Find common subspace, i.e. the null space for the union 
    tmp <- svd(cbind(mybasis[[i]],mybasis[[j]]))
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
          myterms   <- myterms[-k]
          mydesigns <- mydesigns[-k]
          mybasis   <- mybasis[-k]
          Nparm     <- Nparm[-k]
          relations <- relations[-k,-k,drop=FALSE]            
        }
      }
    } else {
      # Design are neither linearly independent nor ordered
      # Parametrize common space (ie. minimum) (using design with lowest index)
      if (i<j) {
        NullSpace <- mybasis[[i]]%*%tmp$v[1:Nparm[i],tmp$d<=eps,drop=FALSE]
      } else {
        NullSpace <- mybasis[[j]]%*%tmp$v[(Nparm[i]+1):(Nparm[i]+Nparm[j]),tmp$d<=eps,drop=FALSE]
      }
      tmp <- svd(NullSpace,nv=0)
      NullSpace <- tmp$u[,tmp$d>eps,drop=FALSE]
      # Look for minimum among the existing terms
      for (k in setdiff((1:M),c(i,j))) {
        if ((NullDim==Nparm[k]) && (NullDim==sum(svd(cbind(NullSpace,mybasis[[k]]),nu=0,nv=0)$d>eps))) {
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
        myterms    <- c(myterms,tmp)
        Nparm      <- c(Nparm,NullDim)
        myisfactor <- c(myisfactor,myisfactor[i] & myisfactor[j])
        if (myisfactor[i] && myisfactor[j]) {
          # Experimental: If both are factors, then use categorical minimum
          tmp <- eval(substitute(with(df,LabApplStat::minimum(x,y)),
                                 list(df=data,x=as.name(sub("[","",sub("]","",myterms[i],fixed=TRUE),fixed=TRUE)),
                                      y=as.name(sub("[","",sub("]","",myterms[j],fixed=TRUE),fixed=TRUE)))))
          tmp <- svd(model.matrix(~0+tmp))
          mydesigns <- c(mydesigns,list(tmp$u[,tmp$d>eps,drop=FALSE]%*%diag(tmp$d[tmp$d>eps],nrow=sum(tmp$d>eps))%*%t(tmp$v[,tmp$d>eps,drop=FALSE])))
        } else {
          # otherwise use orthonormal representation of nullspace
          mydesigns <- c(mydesigns,list(NullSpace))
        }
        mybasis        <- c(mybasis,list(NullSpace))
        relations[i,j] <- relations[j,i] <- myterms[M]
        relations      <- cbind(rbind(relations,rep(NA,M-1)),rep(NA,M))
        relations[M,M] <- "="
      }
    }
    # End while() loop
  }

  # ----------------------------------------------------
  # Find sequential ordering of the terms
  # Uses and modifies the variables: 
  #   myorder, M, myterms, mydesigns, mybasis, Nparm, relations
  # ----------------------------------------------------
  
  # Initialize and find basis
  myorder <- rep(NA,M)
  
  # Loop through all variables
  for (k in 1:M) {
    # find next variable in a sequential ordering of the variables
    if (k==1) {
      # Note: which.min() selects the index of the first minimum, which hence 
      #       complies with the initial order as much as possible.
      #      i <- which.min(M*apply(relations==">",1,sum)-apply(relations=="<",1,sum))
      i <- which.min(1*apply(relations==">",1,any))
    } else {
      #      i <- which.min(M*is.element(1:M,myorder[1:(k-1)])+
      #                     M*apply(relations[,-myorder[1:(k-1)],drop=FALSE]==">",1,sum)-
      #                     apply(relations[,-myorder[1:(k-1)],drop=FALSE]=="<",1,sum))
      i <- which.min(2*is.element(1:M,myorder[1:(k-1)])+
                       apply(relations[,-myorder[1:(k-1)],drop=FALSE]==">",1,any))
    }
    myorder[k] <- i
  }
  
  # Reorder variables
  myterms   <- myterms[myorder]
  Nparm     <- Nparm[myorder]
  relations <- relations[myorder,myorder]
  mydesigns <- mydesigns[myorder]
  mybasis   <- mybasis[myorder]

  
  # -------------------------------------------------------
  # Remove underlying designs and find degrees of freedom
  # -------------------------------------------------------

  # Make orthonormal basis, where nested designs are removed from the designs
  # Note that these are not necessarily pairwise orthogonal, e.g. in case of collinearity
  # We also update the associated singular values
  if (M>0) for (i in M:1) if (ncol(mybasis[[i]])>0) {
    tmp <- (relations[i,]==">")
    if (any(tmp)) {
      # combine preceding subspaces
      tmp <- svd(do.call("cbind",mybasis[which(tmp)]),nv=0)
      B <- tmp$u[,tmp$d>eps,drop=FALSE]
      # update orthonormal designs
      tmp <- svd(mybasis[[i]] - B%*%(t(B)%*%mybasis[[i]]),nv=0)
      mybasis[[i]] <- tmp$u[,tmp$d>eps,drop=FALSE]
    }
  }
  
  # Compute degrees of freedom
  mydf <- unlist(lapply(mybasis,ncol))
  
  # -------------------------------------------------------
  # Investigate orthogonality by computing inner products
  # -------------------------------------------------------

  # Compute inner products
  inner <- matrix(NA,M,M)
  for (i in 1:M) for (j in 1:M) {
    inner[i,j] <- round(sum(c(t(mybasis[[i]])%*%mybasis[[j]])^2),floor(-log10(eps)))
  }
  
  # Issue warning for non-orthogonal designs
  if (any(inner[upper.tri(inner)]!=0)) warning("Design is non-orthogonal: Sum-of-Squares and p-values may depend on order of terms.")

  
  # -----------------------------------
  # Extend with the identity variable
  # -----------------------------------
  
  # TO DO: What is the identity variable is already included!?
  #        Then [I] will have df=0, but what does this imply??
  myterms.random <- c(myterms.random,"[I]")
  myterms        <- c(myterms,"[I]")
  Nparm          <- c(Nparm,N)
  mydf           <- c(mydf,Nparm[M+1]-sum(svd(do.call("cbind",mybasis),nu=0,nv=0)$d>eps))
  # TO DO: Following line is only necessarily correct when mydf("[I]")>0
  #        However, relations.ghost take care of this
  relations <- cbind(rbind(relations,matrix(">",1,M)),matrix(c(rep("<",M),"="),M+1,1))
  M <- M+1
  
  # Make ghost relations, where terms with df=0 are removed.
  # This is used to find testable hypothesis, i.e. collapses that don't remove non-trivial minima
  relations.ghost <- relations
  relations.ghost[mydf==0,] <- "df=0"
  relations.ghost[,mydf==0] <- "df=0"
  
  # -----------------
  # Identify arrows
  # -----------------
  
  # Should be run after computation of degrees of freedom
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

  # ----------------------------------
  # Compute summaries and statistics
  # ----------------------------------

  # find terms to be removed in the collinearity analysis
  # NB: placed here to ensure same ordering as in myterms
  myterms.remove <- setdiff(myterms,attr(terms(keep),"term.labels"))
  if (attr(terms(keep),"intercept")==1) myterms.remove <- setdiff(myterms.remove,"1")
  myterms.remove <- setdiff(myterms.remove,"[I]")
  
  # Find projections onto basis functions and make Type-I F-tests
  # New: Introduce tau2 to contain cummulated variance estimates
  pvalue <- matrix(NA,M,M)
  ii     <- 1+length(myterms.remove)
  SS     <- matrix(NA,ii,M)
  MSS    <- matrix(NA,ii,M)
  df.tmp <- matrix(NA,ii,M)
  sigma2 <- rep(NA,M); names(sigma2) <- myterms.random
  ranef_mat <- diag(1,nrow=length(myterms.random))
  rownames(ranef_mat) <- colnames(ranef_mat) <- myterms.random
  
  if (!is.null(y)) {
    # find orthogonal basis
    if (M>1) for (i in 1:ii) {
      for (j in setdiff(1:(M-1),match(myterms.remove[0:(i-1)],myterms))) {
        # terms to remove 
        k <- setdiff(0:(j-1),c(0,match(myterms.remove[0:(i-1)],myterms)))
        # find basis for orthogonal complement
        if (length(k)==0) {
          A <- mydesigns[[j]]
        } else {
          tmp <- svd(do.call("cbind",mybasis[k]),nv=0)
          B <- tmp$u[,tmp$d>eps,drop=FALSE]
          A <- mydesigns[[j]]-B%*%(t(B)%*%mydesigns[[j]])
        }
        if (ncol(A)>0) {
          tmp <- svd(A,nv=0)
          A <- tmp$u[,tmp$d>eps,drop=FALSE]
          if ((i==1) && is.element(myterms[j],myterms.random)) {
            # update ranef_mat matrix using the orthogonalized designs
            for (k in setdiff(myterms.random,"[I]")) {
              tmp <- t(A)%*%mydesigns[[match(k,myterms)]]
              ranef_mat[myterms[j],k] <- sum(diag(tmp%*%t(tmp))) / mydf[match(k,myterms)]
            }
            ranef_mat[myterms[j],"[I]"] <- sum(diag(t(A)%*%A)) / mydf[match(k,myterms)]
          }
        } else {
          A <- matrix(0,N,0)
        }
        # compute sum-of-squares
        SS[i,j] <- sum(c(y%*%A)^2)
        df.tmp[i,j] <- ncol(A)
      }
    }
    # Residual SS
    SS[,M]     <- sum(y^2) - rowSums(SS[,1:(M-1),drop=FALSE],na.rm=TRUE)
    df.tmp[,M] <- length(y) - rowSums(df.tmp[,1:(M-1),drop=FALSE],na.rm=TRUE)
    
    # Compute MSS
    MSS <- SS/df.tmp

    # Compute sigma2 for random effects
    sigma2 <- solve(ranef_mat,MSS[1,match(myterms.random,myterms)])
    
    # F-tests only for terms with positive degrees of freedom, which have
    # at least one term nested within them.
    for (i in which((mydf>0) & apply(relations.ghost,1,function(x){sum(x=="<-")})>0)) {
      # Find candidates for denominator: random effects with positive degrees of freedom, which are "<" or "<-"
      j <- (is.element(myterms,myterms.random) & (mydf>0) & is.element(relations.ghost[i,],c("<","<-")))
      # Remove candidates that are nested within some of the other candidates
      if (sum(j)>0) j[j] <- j[j] & (!apply(relations.ghost[j,j,drop=FALSE],2,function(x){any(is.element(x,c("<","<-")))}))
      # Make F-test if only one(!) candidate is left
      if (sum(j)==1) {
        j <- which(j)
        # allocate p-value to an "<-" with df>0, and in which j is nested
        k <- which((mydf>0) & (relations[i,]=="<-") & (is.element(relations[j,],c(">","->","="))))
        # if not possible, then relax conditions df>0
        if (length(k)==0) {
          k <- which((relations[i,]=="<-") & (is.element(relations[j,],c(">","->","="))))
        }
        # in rare cases (!?) with more than one k, simply use the first
        pvalue[i,k[1]] <- 1-pf(MSS[1,i]/MSS[1,j],mydf[i],mydf[j])
      }
    }
  }


  # Compute variance-covariance matrices for each fixed effect 
  # relative to each random effect (for variance=1)
  myvarcov <- vector("list",0)
  for (i in which(!is.element(myterms,myterms.random))) {
    varcov <- vector("list",0)
    for (j in which(is.element(myterms,setdiff(myterms.random,"[I]")))) {
      # compute variance-covariance matrix for term i relative to term j
      tmp <- solve(t(mydesigns[[i]])%*%mydesigns[[i]],t(mydesigns[[i]])%*%mydesigns[[j]])
      varcov <- c(varcov,list(round(tmp%*%t(tmp),-log10(eps))))
    }
    # take care of [I]
    varcov <- c(varcov,list(round(solve(t(mydesigns[[i]])%*%mydesigns[[i]]),-log10(eps))))
    # add names of random effects
    names(varcov) <- myterms.random
    # save results in list
    myvarcov <- c(myvarcov,list(varcov))
  }
  names(myvarcov) <- setdiff(myterms,myterms.random)

    
  # Make data frame with coordinates of nodes in a Sugiyama layout ----
  from <- 1+(which(relations=="<-")-1)%/%length(myterms)
  to   <- 1+(which(relations=="<-")-1)%%length(myterms)
  ii   <- order(from,to,decreasing=TRUE)
  g <- ggraph::create_layout(
    data.frame(from=as.character(from[ii]),to=as.character(to[ii])),
    "sugiyama")
  coordinates           <- data.frame(x=max(g$y)-g$y,y=g$x)
  rownames(coordinates) <- myterms[as.numeric(g$name)]

      
  # Return result ----
  names(myterms) <- names(Nparm) <- names(mydf) <- 
    colnames(SS) <- colnames(MSS) <- rownames(relations) <- colnames(relations) <- 
    rownames(pvalue) <- colnames(pvalue) <- myterms
  rownames(SS) <- rownames(MSS) <- c("-",myterms.remove)
  rownames(inner) <- colnames(inner) <- names(mybasis) <- myterms[-M]
  return(structure(list(response=!is.null(y),terms=myterms,random.terms=myterms.random,
                        relations=relations,inner=inner,
                        Nparm=Nparm,df=mydf,
                        SS=SS,MSS=MSS,pvalue=pvalue,
                        sigma2=sigma2, varcov=myvarcov,
                        coordinates=coordinates),
                   class="designDiagram"))
}
