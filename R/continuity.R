#' Assessment of continuity
#' @param n Number of tips
#' @param epsmax Maximum value of perturbation
#' @param repeats Number of repeats to perform
#' @param showProgress Whether to show a progress bar
#' @returns Continuity assessment
#' @export
continuity = function(n=10,epsmax=0.1,repeats=100,showProgress=T) {
  res=matrix(NA,repeats,4)
  if (showProgress) pb <- txtProgressBar(min=0,max=repeats,style = 3)
  for (i in 1:repeats) {
    if (showProgress) setTxtProgressBar(pb, i)
    set.seed(i)
    p1=rtree(n)
    r1=phylo2contphylo(p1)
    r2=r1
    new=r1$coordinates[-1,]
    eps=runif(1,0,epsmax)
    w=sample.int(nrow(new)*2,1)
    if (w<=nrow(new)) new[w,1]=new[w,1]+runif(1,-eps,eps) else
    new[w-nrow(new),2]=log(abs(exp(new[w-nrow(new),2])+runif(1,-eps,eps)))
    new=rbind(c(0,0),new)
    r2$coordinates=new
    p2=contphylo2phylo(r2)
    d1=treediff(p1,p2)#my own distance
    d2=phangorn::KF.dist(p1,p2)#branch score distance
    d3=distory::dist.multiPhylo(list(p1,p2))#geodesic distance from Billera et al 2001
    res[i,]=c(eps,d1,d2,d3)
    colnames(res)<-c('eps','mydist','scoredist','geodesic')
  }
  pairs(res)
  if (showProgress) close(pb)
  return(res)
}
