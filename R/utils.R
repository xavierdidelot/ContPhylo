#' Generate a random contphylo
#' @param n Number of leaves

#' @returns Object of class contphylo
#' @param type Type of tree: 1 for ultrametric, 2 for rooted, 3 for unrooted.
#' @export
rcontphylo = function(n,type=2) {
  r=list()
  r$coordinates=rbind(c(0,0),matrix(rnorm((n-1)*2),n-1,2))
  r$tip.label=paste0('t',1:n)
  if (type==1) r$coordinates=r$coordinates[,1,drop=F]
  if (type==3) r$coordinates[n,1]=0
  r$type=type
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
