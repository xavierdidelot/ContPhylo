#' Simulate a random tree
#' @param n Number of leaves
#' @param model Model of tree: 1 for ultrametric, 2 for rooted, 3 for unrooted
#' @param scale Scale parameter of exponential distribution for branch lengths
#' @returns Random tree
#' @export
randtree=function(n,model,scale) {
  if (model==1) {y=rcoal(n);y$edge.length=y$edge.length*scale;return(y)}
  if (model==2) return(rtopology(n,rooted=T,br=function (n) rexp(n,rate=1/scale)))
  if (model==3) return(rtopology(n,rooted=F,br=function (n) rexp(n,rate=1/scale)))
}
