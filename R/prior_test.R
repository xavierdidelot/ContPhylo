#' Test that a given set of trees (eg output from mcmc()) is distributed according to prior distribution
#' @param trees List of trees
#' @param model Model of tree: 1 for ultrametric, 2 for rooted, 3 for unrooted
#' @param scale Scale parameter of exponential distribution for branch lengths
#' @param analyse_topology Whether or not to analyse the topology
#' @export
prior_test=function(trees,model=2,scale=1,analyse_topology=T) {
  n=Ntip(trees[[1]])
  par(mfrow=c(3,2))
  sbl=sapply(trees,function(x) return(sum(x$edge.length)))
  plot(sbl,type='l',main='Sum branch lengths',xlab='',ylab='')
  hist(sbl,main='Sum branch lengths',xlab='',probability=T,breaks=20)
  xs=seq(min(sbl),max(sbl),length.out=1000)
  dcoal=function(x) {(n-1)/(2*scale)*exp(-x/(2*scale))*(1-exp(-x/(2*scale)))^(n-2)}
  par(xpd=NA)
  if (model==1) lines(xs,dcoal(xs))
  if (model==2) lines(xs,dgamma(xs,shape=(2*n-2),scale=scale))
  if (model==3) lines(xs,dgamma(xs,shape=(2*n-3),scale=scale))
  par(xpd=F)
  ess=coda::effectiveSize(coda::as.mcmc(sbl))
  print(sprintf('ESS of sum of branch lengths: %f',ess))

  b1=sapply(trees,function(x) {wi=which(x$tip.label==trees[[1]]$tip.label[1]);w=which(x$edge[,2]==wi);return(x$edge.length[w])})
  #b1=as.vector(sapply(trees,function(x) {w=which(x$edge[,2]<=n);return(x$edge.length[w])}))
  plot(b1,type='l',main='Single branch length',xlab='',ylab='')
  hist(b1,main='Single branch length',xlab='',probability=T,breaks=20)
  xs=seq(min(b1),max(b1),length.out=1000)
  if (model==1) lines(xs,sapply(xs,function (t) p_t(t,n,scale)))
  if (model>1)  lines(xs,dexp(xs,rate=1/scale))
  ess=coda::effectiveSize(coda::as.mcmc(b1))
  print(sprintf('ESS of single branch length: %f',ess))

  clusters=rep(1,length(trees))
  if (analyse_topology) {
    clusters=topologies(trees)
    counts=table(clusters)
    print(counts)
    expected=round(exp(lfactorial(2*n-2)-(n-1)*log(2)-lfactorial(n-1)))
    if (length(counts)<expected) print(sprintf('Found %d topologies but was expecting %d.',length(counts),expected))
    t=chisq.test(counts)
    print(t)

    #Trace and barplot of topology
    plot(clusters,xlab='',ylab='')
    names(counts)<-1:length(counts)
    barplot(counts)
  }

  if (F) {
  #Checking that consensus tree looks completely star-like
  m=trees
  class(m)<-'multiPhylo'
  cons_topo=consensus(m, p = 0.5, rooted = T)
  cons=phangorn::add_edge_length(cons_topo, m, rooted = T,fun=mean)
  plot(cons,show.tip.label = T)
  axisPhylo(1)
  }
}

#This is the distribution of a terminal branch length under the coalescent
p_t <- function(t, n, theta) {
  stopifnot(n >= 2, theta > 0)
  lambda <- 1/theta * choose(1:n, 2)          # lambda[j] for j = 1..n (lambda[1] = 0, unused)
  # inner term for given k: vector over t
  dens_k <- function(k) {
    js <- k:n
    out <- numeric(length(t))
    for (j in js) {
      rs <- setdiff(k:n, j)
      coef <- lambda[j] * prod(lambda[rs] / (lambda[rs] - lambda[j]))
      out <- out + coef * exp(-lambda[j] * t)
    }
    out
  }
  total <- numeric(length(t))
  for (k in 2:n) {
    total <- total + (k - 1) * dens_k(k)
  }
  2 / (n * (n - 1)) * total
}
