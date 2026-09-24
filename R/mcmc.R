#' MCMC algorithm
#' @param n Number of tips
#' @param niter Number of MCMC iterations
#' @param thin Thining interval
#' @param shape Shape of the gamma distribution for branch lengths
#' @param scale Scale of the gamma distribution for branch lengths
#' @param showProgress Whether to show a progress bar
#' @returns MCMC result
#' @export
mcmc = function(n=10,niter=1e3,thin=max(1,round(niter/1e3)),shape=1,scale=1,showProgress=T) {
  tree=rtree(n)
  r=phylo2contphylo(tree)
  p=prior(r,shape,scale)
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
    for (move in 1:2) {
      r2=r
      new=r$coordinates[-1,]
      w=sample.int(nrow(new),1)
      new[w,move]=new[w,move]+rnorm(1,0,sd[move])
      new=rbind(c(0,0),new)
      r2$coordinates=new
      p2=prior(r2,shape,scale)
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
      res[[i/thin]]$sd=sd
      res[[i/thin]]$paccepted=naccepted/i
    }
  }
  if (showProgress) close(pb)
  return(res)
}

#' MCMC algorithm using package adaptMCMC
#' @param n Number of tips
#' @param niter Number of MCMC iterations
#' @param thin Thining interval
#' @param shape Shape of the gamma distribution for branch lengths
#' @param scale Scale of the gamma distribution for branch lengths
#' @param showProgress Whether to show a progress bar
#' @returns MCMC result
#' @export
mcmc2 = function(n=10,niter=1e3,thin=max(1,round(niter/1e3)),shape=1,scale=1,showProgress=T) {
  target = function(v) {
    r=list()
    r$coordinates=rbind(c(0,0),matrix(v,ncol=2,byrow=F))
    r$tip.label=paste0('t',1:n)
    class(r)<-'contphylo'
    p=prior(r,shape,scale)
    return(p)
  }
  init=runif(n*2-2)
  r=adaptMCMC::MCMC(p=target,n=niter,init=init,adapt=T,acc.rate=0.234,showProgressBar=showProgress)
  res=vector('list',niter/thin)
  l=list()
  l$tip.label=paste0('t',1:n)
  class(l)<-'contphylo'
  for (i in seq(thin,niter,thin)) {
    l$coordinates=rbind(c(0,0),matrix(r$samples[i,],ncol=2,byrow=F))
    res[[i/thin]]=contphylo2phylo(l)
  }
  res[[1]]$cov.jump=r$cov.jump#interesting
  return(res)
}

