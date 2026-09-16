#' Convert a contphylo into a phylo
#' @param r Object of class contphylo
#' @returns Object of class phylo
#' @export
contphylo2phylo = function(r) {
  stopifnot(inherits(r, "contphylo"))
  coords=exp(r$coordinates)
  n=nrow(coords)
  if (n<2) stop('Need at least two leaves.')
  m=matrix(NA,n*2-1,4)#child1,child2,father,dist_to_father
  for (i in 2:n) {
    if (i==2) graft=1 else graft=sample(c(1:(i-1),n+(1:(i-2))),1)
    graft_dist=rexp(1)#TODO contour traversal based on coords[i,1]
    if (!is.na(m[graft,4]) && graft_dist>m[graft,4]) graft_dist=m[graft,4]/2#TODO remove
    m[i,3]=n+i-1
    m[i,4]=coords[i,2]
    m[n+i-1,1]=graft
    m[n+i-1,2]=i
    m[n+i-1,3]=m[graft,3]
    m[n+i-1,4]=m[graft,4]-graft_dist
    m[graft,3]=n+i-1
    m[graft,4]=graft_dist
  }
  root=which(is.na(m[,3]))
  if (root!=n+1) {
    #swap root and n+1
    row=m[n+1,]
    m[n+1,]=m[root,]
    m[root,]=row
    m[which(m[,1:3]==n+1, arr.ind=T)]=-1
    m[which(m[,1:3]==root,arr.ind=T)]=n+1
    m[which(m[,1:3]==-1. ,arr.ind=T)]=root
  }
  p=list()
  p$Nnode=n-1
  p$edge=cbind(m[,3],1:(n*2-1))
  p$edge=p$edge[-(n+1),]
  p$edge.length=m[-(n+1),4]
  p$tip.label=r$tip.label
  class(p)='phylo'
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
