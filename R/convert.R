#' Convert a contphylo into a phylo
#' @param r Object of class contphylo
#' @returns Object of class phylo
#' @export
contphylo2phylo = function(r) {
  stopifnot(inherits(r, "contphylo"))
  n=nrow(r$coordinates)
  p=rtree(n)#TODO
  p$tip.label=r$tip.label
  return(p)
}

#' Convert a phylo into a contphylo
#' @param p Object of class phylo
#' @returns Object of class contphylo
#' @export
phylo2contphylo = function(p) {
  stopifnot(inherits(p, "phylo"))
  n=Ntip(p)
  coords=rbind(c(0,0),matrix(rnorm((n-1)*2),n-1,2))#TODO
  r=list()
  r$coordinates=coords
  r$tip.label=p$tip.label
  class(r)<-'contphylo'
  return(r)
}

#' Generate a random contphylo
#' @param n Number of leaves
#' @returns Object of class contphylo
#' @export
rcontphylo = function(n) {
  r=list()
  r$coordinates=rbind(c(0,0),matrix(rnorm((n-1)*2),n-1,2))
  r$tip.label=paste0('t',1:n)
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
  plot(c(),c(),xlab='',ylab='',xlim=range(x$coordinates[,1]),ylim=range(x$coordinates[,2]))
  text(x$coordinates[,1],x$coordinates[,2],x$tip.label)
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
