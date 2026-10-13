#' Generate a random contphylo
#' @param n Number of leaves
#' @param model Model of tree: 1 for ultrametric, 2 for rooted, 3 for unrooted.
#' @returns Object of class contphylo
#' @export
rcontphylo = function(n,model=2) {
  r=list()
  if (model==1) {
    dim=n-1
    rn=rnorm(dim)
    r$coordinates=matrix(c(0,rn),nrow=n,ncol=1)
  }
  if (model==2) {
    dim=2*n-2
    rn=rnorm(dim)
    r$coordinates=rbind(c(0,0),matrix(rn,n-1,2))
  }
  if (model==3) {
    dim=2*n-3
    rn=rnorm(dim)
    r$coordinates=rbind(c(0,0),matrix(rn[-1],n-2,2),c(0,tail(rn,1)))
  }
  r$p=sum(dnorm(rn,log=T))
  r$tip.label=paste0('t',1:n)
  r$model=model
  class(r)<-'contphylo'
  return(r)
}

#' Plotting method for contphylo
#' @param x Object of class contphylo
#' @param ... Additional parameters are passed on
#' @return Plot of contphylo
#' @export
plot.contphylo = function(x, ...) {
  stopifnot(inherits(x, "contphylo"))
  p=contphylo2phylo(x)
  plot(p,...)
}

#' Print function for contphylo objects
#' @param x Object of class contphylo
#' @param ... Passed on
#' @return Print out details of contphylo object
#' @export
print.contphylo <- function(x, ...) {
  stopifnot(inherits(x, "contphylo"))
  cat(sprintf('contphylo object with %d leaves\n',nrow(x$coordinates)))
  invisible(x)
}

#' Summary function for contphylo objects
#' @param object Output from ClonalFrame2
#' @param ... Passed on
#' @return Print out details of contphylo object
#' @export
summary.contphylo <- function(object, ...) {
  stopifnot(inherits(object, "contphylo"))
  print(object,...)
  invisible(object)
}

#' Sum of branch lengths equation
#' @param x Object of class contphylo
#' @return Sum of branch lengths
#' @export
sum_bra_len <- function(x) {
  stopifnot(inherits(x, "contphylo"))
  if (x$model==1) x=phylo2contphylo(contphylo2phylo(x),2)
  l=0
  coords=x$coordinates
  for (i in 2:nrow(coords)) {
    l=l+exp(coords[i,2])+max(0,coords[i,1]-2*l)+max(0,-coords[i,1])
  }
  return(l)
}

#' Distance between two trees, based on Mean Square Error of pairwise distances
#' @param t1 First tree
#' @param t2 Second tree
#' @return MSE of pairwise distances
#' @export
treediff = function(t1,t2) {
  d1=cophenetic(t1)
  d2=cophenetic(t2)
  return(mean((d1-d2)^2))
}
