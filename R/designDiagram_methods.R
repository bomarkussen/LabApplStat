# --------------------------
# Implemented methods:
# 1. print
# 2. summary
# 3. plot
# --------------------------

# --------------------------
# designDiagram print
# --------------------------

print.designDiagram <- function(x,...) {
  # Table of dimensions
  if (any(x$inner[upper.tri(x$inner)]!=0)) {
    cat("Non-orthogonal design with dimensions:\n")
  } else {
    cat("Orthogonal design with dimensions:\n")
  }
  print(rbind(Nparm=x$Nparm,df=x$df,Ncol=x$collinearities))
}

# --------------------------
# designDiagram summary
# --------------------------

summary.designDiagram <- function(object,...) {
  # Table of dimensions
  if (any(object$inner[upper.tri(object$inner)]!=0)) {
    cat("Non-orthogonal design with dimensions:\n")
  } else {
    cat("Orthogonal design with dimensions:\n")
  }
  print(rbind(Nparm=object$Nparm,df=object$df,Ncol=object$collinearities))
  
  # Table of inner products
  if (any(object$inner[upper.tri(object$inner)]!=0)) {
    cat("\n")
    cat("Note: Sum-of-Squares and p-values may depend on order of terms in an non-orthogonal design.\n")
    cat("\n")
    cat("Total inner products between subspaces (used to decide orthogonality):\n")
    print(object$inner)
  }
  
  # Table of relations
  cat("\n")
  cat("Table of relations:\n")
  print(object$relations)
  
  # Additional output if response is present
  cat("\n")
  if (!object$response) {
    cat("No response variable specified.\n")
  } else {
    cat("Orthogonal decomposition of response variable:\n")
    print(rbind(SS=object$SS,MSS=object$MSS))
    cat("\n")
    cat("P-values for F-tests against nested random effect:\n")
    print(signif(object$pvalue,6))
  }
}

# --------------------------
# designDiagram plot
# --------------------------

plot.designDiagram <- function(x,circle="none",pvalue=(circle=="MSS"),kill.intercept=TRUE,
                               diam=80,color=ifelse(circle=="MSS","lightblue","lightgreen"),
                               border=c(0,0.1,0.1,0.2),
                               ...) {
  
  # sanity check
  if (!is.element(circle,c("none","SS","MSS"))) stop("circle-argumente must be either none, SS, or MSS")
  if ((!x$response) & (is.element(circle,c("SS","MSS")))) {
    circle <- "none"
    warning("Sum of Squares unavailable")
  }
  if (length(border)==1) border <- rep(border,4)

  # set-up basic graph
  N <- length(x$terms)
  myedges <- c(sapply(which(x$relations=="<-"),function(z){c(1+(z-1)%/%N,1+(z-1)%%N)}))
  g <- make_graph(myedges,directed=TRUE) %>% set_vertex_attr("label",value=rep(" ",N))
  if (pvalue) E(g)$label <- as.character(signif(x$pvalue[x$relations=="<-"],digits=3))
  lay1 <- layout_with_sugiyama(g,attributes="all")

  # takeout layout and turn the vertical direction: NOT USED ANYMORE
  #tmp     <- get.graph.attribute(lay1$extd_graph)$layout
  #y.max   <- 1+max(tmp[,2])
  #tmp[,2] <- y.max-tmp[,2]
  #lay1$extd_graph <- set_graph_attr(lay1$extd_graph,"layout",tmp)
  #tmp     <- lay1$layout
  #tmp[,2] <- y.max-tmp[,2]
  #lay1$layout <- tmp
  #tmp     <- lay1$layout.dummy
  #tmp[,2] <- y.max-tmp[,2]
  #lay1$layout.dummy <- tmp
  
  # Plain design diagram
  if (circle=="none") {
    plot(lay1$extd_graph,
         vertex.size=c(rep(30,length(lay1$layout)/2),rep(0,length(lay1$layout.dummy)/2)),
         vertex.color=NA,vertex.frame.color=NA,
         rescale=FALSE,
         xlim=c(min(lay1$layout[,1])-border[2],max(lay1$layout[,1])+border[4]),
         ylim=c(min(lay1$layout[,2])-border[1],max(lay1$layout[,2])+border[3]))
  }
  
  # Design diagram with circles
  if (is.element(circle,c("SS","MSS"))) {
    # find circle diameters
    if (circle=="SS") {area <- x$SS} else {area <- x$MSS}
    if (kill.intercept & (is.element("1",names(area)))) area["1"] <- 0
    diam <- diam*sqrt(area/max(area,na.rm=TRUE))
    diam[is.na(diam)] <- 0
    diam[is.nan(diam)] <- 0
    # make graph
    plot(lay1$extd_graph,
         vertex.size=c(diam,rep(0,length(lay1$layout.dummy)/2)),
         vertex.color=color,vertex.frame.color=NA,
         rescale=FALSE,
         xlim=c(min(lay1$layout[,1])-border[2],max(lay1$layout[,1])+border[4]),
         ylim=c(min(lay1$layout[,2])-border[1],max(lay1$layout[,2])+border[3]))
  }
  
  # insert vertex names
  for (i in 1:N) {
    if (x$collinearities[i]==0) {
      text(lay1$layout[i,1],lay1$layout[i,2],cex=1.2,
           substitute(x[b]^a,list(x=names(x$terms)[i],a=x$Nparm[i],b=x$df[i])))
    } else {
      tmp <- x$collinearities[i]
      text(lay1$layout[i,1],lay1$layout[i,2],cex=1.2,
           substitute(x[b-c]^a,list(x=names(x$terms)[i],a=x$Nparm[i],b=x$df[i]+tmp,c=tmp)))
    }
  }
}
