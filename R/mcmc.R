#' MCMC algorithm
#' @param dat Object of class phyDat containing the genetic data
#' @param model Model of tree: 1 for ultrametric, 2 for rooted, 3 for unrooted
#' @param initTree Initial tree to use as start of MCMC
#' @param niter Number of MCMC iterations
#' @param thin Thining interval
#' @param scale Scale of the branch lengths
#' @param showProgress Whether to show a progress bar
#' @returns MCMC result
#' @export
mcmc = function(dat=phyDatEmpty(5),model=2,initTree,niter=1e3,thin=max(1,round(niter/1e3)),scale=1,showProgress=T) {
  if (missing(initTree)) {
    n=length(dat)
    if (model==1) tree=rcoal(n) else tree=rtree(n)
    tree$edge.length=tree$edge.length*scale
    tree$tip.label=names(dat)
    if (model==3) tree=unroot(tree)
  } else {
    n=Ntip(initTree)
    tree=initTree
  }
  mypml=phangorn::pml(tree=tree,data=dat,model='JC')
  r=phylo2contphylo(tree,model=model)
  p=posterior(mypml,r,scale)
  res=vector('list',niter/thin)
  if (showProgress) pb <- txtProgressBar(min=0,max=niter,style = 3)
  adapt_target=0.234
  adapt_scale=1.0/adapt_target/(1.0-adapt_target)
  sd=rep(1,2)
  adapt_nstart=rep(5.0*adapt_scale,2)
  adapt_lambda=rep(1.0,2)
  naccepted=c(0,0)
  for (i in 1:niter) {
    if (showProgress) setTxtProgressBar(pb, i)
    for (move in 1:(1+(model>1))) {
      r2=r
      w=1+sample.int(n-1-(model==3 && move==1),1)
      r2$coordinates[w,move]=r2$coordinates[w,move]+rnorm(1,0,sd[move])
      p2=posterior(mypml,r2,scale)
      mh=p2-p
      if (log(runif(1))<mh) {
        p=p2
        r=r2
        naccepted[move]=naccepted[move]+1
      }
      prob_accept=min(1,exp(mh))
      sd[move]=sd[move]*exp(0.5*adapt_scale/(i+adapt_nstart[move])*(prob_accept-adapt_target))
      lambda=sd[move]*sd[move]
      if (abs(log(lambda)-log(adapt_lambda[move]))>log(3.0)) {
        #Restart Robbins-Munro
        adapt_lambda[move]=lambda
        adapt_nstart[move]=5.0*adapt_scale-i
      }
    }
    if (i%%thin==0) {
      res[[i/thin]]=contphylo2phylo(r)
      res[[i/thin]]$posterior=p
      res[[i/thin]]$sd=sd
      res[[i/thin]]$paccepted=naccepted/i
    }
  }
  if (showProgress) close(pb)
  return(res)
}

#' MCMC algorithm using package adaptMCMC
#' @param dat Object of class phyDat containing the genetic data
#' @param model Model of tree: 1 for ultrametric, 2 for rooted, 3 for unrooted
#' @param initTree Initial tree to use as start of MCMC
#' @param niter Number of MCMC iterations
#' @param thin Thining interval
#' @param scale Scale of the branch lengths
#' @param showProgress Whether to show a progress bar
#' @returns MCMC result
#' @export
mcmc2 = function(dat=phyDatEmpty(5),model=2,initTree,niter=1e3,thin=max(1,round(niter/1e3)),scale=1,showProgress=T) {
  if (missing(initTree)) {
    n=length(dat)
    if (model==1) tree=rcoal(n) else tree=rtree(n)
    tree$edge.length=tree$edge.length*scale
    tree$tip.label=names(dat)
    if (model==3) tree=unroot(tree)
  } else {
    n=Ntip(initTree)
    tree=initTree
  }
  mypml=phangorn::pml(tree=tree,data=dat,model='JC')
  r=phylo2contphylo(tree,model=model)
  if (model==1) init=c(r$coordinates[-1,1])
  if (model==2) init=c(r$coordinates[-1,1],r$coordinates[-1,2])
  if (model==3) init=c(r$coordinates[c(-1,-n),1],r$coordinates[-1,2])
  target = function(v) {
    r=list()
    if (model==1) r$coordinates=matrix(c(0,v),ncol=1)
    if (model==2) r$coordinates=rbind(c(0,0),matrix(v,ncol=2,byrow=F))
    if (model==3) r$coordinates=rbind(c(0,0),matrix(c(v[1:(n-2)],0,v[(n-1):(2*n-3)]),ncol=2,byrow=F))
    r$tip.label=tree$tip.label
    r$model=model
    class(r)<-'contphylo'
    p=posterior(mypml,r,scale)
    return(p)
  }
  r=adaptMCMC::MCMC(p=target,n=niter,init=init,adapt=T,acc.rate=0.234,showProgressBar=showProgress)
  res=vector('list',niter/thin)
  l=list()
  l$tip.label=tree$tip.label
  l$model=model
  class(l)<-'contphylo'
  for (i in seq(thin,niter,thin)) {
    v=r$samples[i,]
    if (model==1) l$coordinates=matrix(c(0,v),ncol=1)
    if (model==2) l$coordinates=rbind(c(0,0),matrix(v,ncol=2,byrow=F))
    if (model==3) l$coordinates=rbind(c(0,0),matrix(c(v[1:(n-2)],0,v[(n-1):(2*n-3)]),ncol=2,byrow=F))
    res[[i/thin]]=contphylo2phylo(l)
    if (model==3) res[[i/thin]]=unroot(res[[i/thin]])
    res[[i/thin]]$posterior=r$log.p[i]
  }
  res[[1]]$cov.jump=r$cov.jump#interesting
  return(res)
}

