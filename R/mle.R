#' MLE algorithm
#' @param dat Object of class phyDat containing the genetic data
#' @param model Model of tree: 1 for ultrametric, 2 for rooted, 3 for unrooted
#' @param initTree Initial tree to use as start of MCMC
#' @returns MLE result
#' @export
mle = function(dat=phyDatEmpty(5),model=2,initTree) {
  if (missing(initTree)) {
    n=length(dat)
    tree=randtree(n,model,0.01)
    tree$tip.label=names(dat)
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
    tree=contphylo2phylo(r)
    mypml=stats::update(mypml,tree=tree)
    return(mypml$logLik)
  }
  o=optim(init,target,method='BFGS',control=list(fnscale=-1))
  #o=nlm(target,init,fscale=-1,print.level=2);o$par=o$estimate;o$value=o$minimum
  v=o$par
  l=list()
  l$tip.label=tree$tip.label
  l$model=model
  class(l)<-'contphylo'
  if (model==1) l$coordinates=matrix(c(0,v),ncol=1)
  if (model==2) l$coordinates=rbind(c(0,0),matrix(v,ncol=2,byrow=F))
  if (model==3) l$coordinates=rbind(c(0,0),matrix(c(v[1:(n-2)],0,v[(n-1):(2*n-3)]),ncol=2,byrow=F))
  res=contphylo2phylo(l)
  if (model==3) res=unroot(res)
  res$ll=o$value
  return(res)
}

