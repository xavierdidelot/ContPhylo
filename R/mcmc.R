#' MCMC algorithm
#' @param n Number of tips
#' @param niter Number of MCMC iterations
#' @param thin Thining interval
#' @param showProgress Whether to show a progress bar
#' @returns MCMC result
#' @export
mcmc = function(n=10,niter=1e3,thin=max(1,round(niter/1e3)),showProgress=T) {
  tree=rtree(n)
  p=prior(tree)
  r=phylo2contphylo(tree)
  res=vector('list',niter/thin)
  if (showProgress) pb <- txtProgressBar(min=0,max=niter,style = 3)
  for (i in 1:niter) {
    if (showProgress) setTxtProgressBar(pb, i)
    r2=r
    new=r$coordinates[-1,]
    w=sample.int(nrow(new)*2,1)
    if (w<=nrow(new)) new[w,1]=new[w,1]+rnorm(1,0,1e-3) else
      new[w-nrow(new),2]=log(abs(exp(new[w-nrow(new),2])+rnorm(1,0,1e-3)))
    new=rbind(c(0,0),new)
    r2$coordinates=new
    tree2=contphylo2phylo(r2)
    p2=prior(tree2)
    if (log(runif(1))<p2-p) {
      p=p2
      r=r2
      tree=tree2
    }
    if (i%%thin==0) res[[i/thin]]=tree
  }
  close(pb)
  return(res)
}

#' MCMC algorithm using package adaptMCMC
#' @param n Number of tips
#' @param niter Number of MCMC iterations
#' @param thin Thining interval
#' @param showProgress Whether to show a progress bar
#' @returns MCMC result
#' @export
mcmc2 = function(n=10,niter=1e3,thin=max(1,round(niter/1e3)),showProgress=T) {
  target = function(v) {
    r=list()
    r$coordinates=rbind(c(0,0),matrix(v,ncol=2,byrow=F))
    r$tip.label=paste0('t',1:n)
    class(r)<-'contphylo'
    t=contphylo2phylo(r)
    p=prior(t)
    return(p)
  }
  tree=rtree(n)
  init=as.vector(phylo2contphylo(tree)$coordinates[-1,])
  r=adaptMCMC::MCMC(p=target,n=niter,init=init,adapt=T,acc.rate=0.234,showProgressBar=showProgress)
  res=vector('list',niter/thin)
  l=list()
  l$tip.label=paste0('t',1:n)
  class(l)<-'contphylo'
  for (i in seq(thin,niter,thin)) {
    l$coordinates=rbind(c(0,0),matrix(r$samples[i,],ncol=2,byrow=F))
    res[[i/thin]]=contphylo2phylo(l)
  }
  return(res)
}

#' Perform topology test for a given set of trees (eg output from mcmc())
#' @param trees List of trees
#' @returns Topology test p-value
#' @export
topology_test=function(trees=mcmc(n=3)) {
  n=Ntip(trees[[1]])
  clusters=rep(0,length(trees))
  reps=vector('list',length(trees))
  k=0
  for (i in seq_along(trees)) {
    assigned=FALSE
    if (k>0) for (j in 1:k)
      if (phangorn::RF.dist(trees[[i]], reps[[j]],rooted=T) == 0) {clusters[i]=j;assigned=TRUE;break}
    if (!assigned) {k=k+1;reps[[k]]=trees[[i]];clusters[i]=k}
  }
  counts=table(clusters)
  print(counts)
  expected=round(exp(lfactorial(2*n-2)-(n-1)*log(2)-lfactorial(n-1)))
  if (length(counts)<expected) stop(sprintf('Found %d topologies but was expecting %d.',length(counts),expected))
  t=chisq.test(counts)
  print(t)
  invisible(t)
}
