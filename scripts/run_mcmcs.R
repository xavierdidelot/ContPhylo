rm(list=ls())
library(ape)
library(ContPhylo)
source('scripts/fig_label.R')

set.seed(0)

# Canonical integer key for a ranked labeled topology.
# `ref` is the common tip-label ordering shared by all trees.
ranked_key <- function(tree, ref = sort(tree$tip.label), tol = 1e-8) {
  n    <- Ntip(tree)
  tree <- reorder(tree, "postorder")          # node ids unchanged
  lab  <- c(match(tree$tip.label, ref), rep(Inf, n - 1L))  # min ref-id below each node

  p  <- tree$edge[, 1]
  ch <- tree$edge[, 2]
  for (i in seq_along(p)) {                    # postorder: children before parents
    if (lab[ch[i]] < lab[p[i]]) lab[p[i]] <- lab[ch[i]]
  }

  # the two child minima for each internal node
  o   <- order(p)
  cm  <- lab[ch][o]
  par <- p[o][c(TRUE, FALSE)]                  # internal nodes n+1, ..., 2n-1
  a   <- cm[c(TRUE, FALSE)]
  b   <- cm[c(FALSE, TRUE)]

  age <- branching.times(tree)[par - n]        # ordered by node id n+1, ...
  ord <- order(-age)                           # oldest (root) first
  if (any(diff(age[ord]) > -tol))
    warning("Tied (or nearly tied) coalescence times: ranking is ambiguous.")

  as.vector(rbind(pmin(a, b), pmax(a, b))[, ord])
}

key_string <- function(key) paste(key, collapse = " ")

n=5#number of leaves in tree
model=1#model of tree
scale=0.01#scaling factor on branch lengths

pdf('mcmcs.pdf',7,7)
let=1
par(mfrow=c(3,3))
for (model in 1:3) {
  trees=mcmc2(phyDatEmpty(n),model=model,scale=scale,niter=2e5,thin=10)
  trees=trees[(length(trees)/2+1):length(trees)]

  clusters=rep(1,length(trees))
  reps=vector('list',length(trees))
  k=0
  if (model==1) {
    ref   <- sort(trees[[1]]$tip.label)
    keys  <- vapply(trees, function(t) key_string(ranked_key(t, ref)), character(1))
  }
  for (i in seq_along(trees)) {
    assigned=FALSE
    if (k>0) for (j in 1:k) {
      if (model>1  && all.equal(trees[[i]], trees[[reps[[j]]]],use.edge.length = F)) {clusters[i]=j;assigned=TRUE;break}
      if (model==1 && keys[[i]]==keys[[reps[[j]]]]) {clusters[i]=j;assigned=TRUE;break}
    }
    if (!assigned) {k=k+1;reps[[k]]=i;clusters[i]=k}
  }
  counts=table(clusters)
  expected=round(exp(lfactorial(2*n-2)-(n-1)*log(2)-lfactorial(n-1)))
  ma=max(clusters)
  map=sample(1:ma,ma)
  b=counts[map]
  names(b)<-1:length(b)
  barplot(b)
  lines(c(0,2*length(b)),c(1,1)*sum(b)/length(b),col='red')
  fig_label(LETTERS[let],cex=2);let=let+1

  b1=sapply(trees,function(x) {wi=which(x$tip.label==trees[[1]]$tip.label[1]);w=which(x$edge[,2]==wi);return(x$edge.length[w])})
  #b1=as.vector(sapply(trees,function(x) {w=which(x$edge[,2]<=n);return(x$edge.length[w])}))
  hist(b1[b1<=0.05],main='',xlab='',probability=T,breaks=seq(0,0.05,length.out=20))
  xs=seq(0,0.05,length.out=1000)
  if (model==1) lines(xs,sapply(xs,function (t) ContPhylo:::p_t(t,n,scale)),col='red')
  if (model>1)  lines(xs,dexp(xs,rate=1/scale),col='red')
  fig_label(LETTERS[let],cex=2);let=let+1

  sbl=sapply(trees,function(x) return(sum(x$edge.length)))
  hist(sbl[sbl<=0.2],main='',xlab='',probability=T,breaks=seq(0,0.2,length.out=20))
  xs=seq(0,0.2,length.out=1000)
  dcoal=function(x) {(n-1)/(2*scale)*exp(-x/(2*scale))*(1-exp(-x/(2*scale)))^(n-2)}
  par(xpd=NA)
  if (model==1) lines(xs,dcoal(xs),col='red')
  if (model==2) lines(xs,dgamma(xs,shape=(2*n-2),scale=scale),col='red')
  if (model==3) lines(xs,dgamma(xs,shape=(2*n-3),scale=scale),col='red')
  par(xpd=F)
  fig_label(LETTERS[let],cex=2);let=let+1
}
dev.off()
browseURL('mcmcs.pdf')
