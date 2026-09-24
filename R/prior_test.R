#' Test that a given set of trees (eg output from mcmc()) is distributed according to prior distribution
#' @param trees List of trees
#' @param scale Scale parameter of exponential distribution for branch lengths
#' @param analyse_topology Whether or not to analyse the topology
#' @export
prior_test=function(trees=mcmc(n=3),scale=1,analyse_topology=T) {
  n=Ntip(trees[[1]])
  clusters=rep(1,length(trees))
  if (analyse_topology) {
    reps=vector('list',length(trees))
    k=0
    for (i in seq_along(trees)) {
      assigned=FALSE
      if (k>0) for (j in 1:k)
        if (all.equal(trees[[i]], reps[[j]],use.edge.length = F)) {clusters[i]=j;assigned=TRUE;break}
      if (!assigned) {k=k+1;reps[[k]]=trees[[i]];clusters[i]=k}
    }
    counts=table(clusters)
    print(counts)
    expected=round(exp(lfactorial(2*n-2)-(n-1)*log(2)-lfactorial(n-1)))
    if (length(counts)<expected) print(sprintf('Found %d topologies but was expecting %d.',length(counts),expected))
    t=chisq.test(counts)
    print(t)
  }

  par(mfrow=c(3,2))
  sbl=sapply(trees,function(x) return(sum(x$edge.length)))
  plot(sbl,type='l',main='Sum branch lengths',xlab='',ylab='')
  hist(sbl,main='Sum branch lengths',xlab='',probability=T,breaks=20)
  xs=seq(min(sbl),max(sbl),length.out=1000)
  lines(xs,dgamma(xs,shape=(2*n-2),scale=scale))
  ess=coda::effectiveSize(coda::as.mcmc(sbl))
  print(sprintf('ESS of sum of branch lengths: %f',ess))

  b1=sapply(trees,function(x) {w=which(x$edge[,2]==1);return(x$edge.length[w])})
  plot(b1,type='l',main='Single branch length',xlab='',ylab='')
  hist(b1,main='Single branch length',xlab='',probability=T,breaks=20)
  xs=seq(min(b1),max(b1),length.out=1000)
  lines(xs,dexp(xs,rate=1/scale))
  ess=coda::effectiveSize(coda::as.mcmc(b1))
  print(sprintf('ESS of single branch length: %f',ess))

  #Trace of topology
  ma=max(clusters)
  map=sample(1:ma,ma)
  plot(map[clusters],xlab='',ylab='')

  #Checking that consensus tree looks completely star-like
  m=trees
  class(m)<-'multiPhylo'
  cons_topo=consensus(m, p = 0.5, rooted = T)
  cons=phangorn::add_edge_length(cons_topo, m, rooted = T,fun=mean)
  plot(cons,show.tip.label = T)
  axisPhylo(1)
}
