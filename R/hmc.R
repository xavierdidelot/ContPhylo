#' Compute prior using stan
#' @param r Object of class ContPhylo
#' @param scale Scale of the exponential distribution for branch lengths
#' @returns Log of the prior
#' @export
prior_stan <- function(r,scale=1) {
  a=r$coordinates[-1,1]
  b=r$coordinates[-1,2]
  n <- length(a) + 1L
  invisible(capture.output(suppressMessages({
  mod <- cmdstan_model(system.file('stan','prior.stan',package='ContPhylo',mustWork = T))
  fit <- mod$sample(
    data=list(n=n,scale=scale),
    chains=1,parallel_chains=1,iter_warmup=0,iter_sampling=1,
    fixed_param=TRUE,init=list(list(a=a,b=b)),seed=1L,refresh=0)
  })))
  log_prior=as.numeric(fit$draws(variables='log_prior',format='matrix'))
  log_prior=log_prior+(n-1)*log(2)+lfactorial(n-1)-lfactorial(2*n-2)
  return(log_prior)
}

#' HMC algorithm using stan
#' @param dat Object of class phyDat containing the genetic data
#' @param initTree Initial tree to use as start of MCMC
#' @param niter Number of iterations
#' @param thin Thining interval
#' @param scale Scale of the exponential distribution for branch lengths
#' @param showProgress Whether to show a progress bar
#' @returns HMC result
#' @export
hmc = function(dat=phyDatEmpty(10),initTree,niter=1e3,thin=max(1,round(niter/1e3)),scale=1,showProgress=T) {
  if (missing(initTree)) {
    n=length(dat)
    tree=rtree(n)
    tree$tip.label=names(dat)
  } else {
    n=Ntip(initTree)
    tree=initTree
  }
  r=phylo2contphylo(tree)
  a=r$coordinates[-1,1]
  b=r$coordinates[-1,2]
#  invisible(capture.output(suppressMessages({
    mod <- cmdstan_model(system.file('stan','prior.stan',package='ContPhylo',mustWork = T))
    fit <- mod$sample(
      data=list(n=n,scale=scale),
      chains=1,parallel_chains=1,iter_warmup=niter/2,iter_sampling=niter/2,
      init=list(list(a=a,b=b)),seed=1L,refresh=ifelse(showProgress,thin,0))
#  })))
  mat=fit$draws(format = 'matrix')
  lp=fit$lp()
  res=vector('list',niter/2/thin)
  l=list()
  l$tip.label=tree$tip.label
  class(l)<-'contphylo'
  for (i in seq(thin,niter/2,thin)) {
    l$coordinates=rbind(c(0,0),matrix(mat[i,2:(2*n-1)],ncol=2,byrow=F))
    res[[i/thin]]=contphylo2phylo(l)
    res[[i/thin]]$posterior=lp[i]
  }
  return(res)
}
