#' Topology identification
#' @param trees List of phylogenetic tree
#' @param ranked Whether or not to identify ranked topologies (aka histories)
#' @returns Index of topologies
#' @export
topologies=function(trees,ranked=F) {
  #trees <- ape::.compressTipLabel(trees)
  if (ranked) keys  <- vapply(trees, ranked_key, "")
  else keys  <- vapply(trees, unranked_key, "")
  clusters <- match(keys, unique(keys))
  #ma=max(clusters)
  #remap=sample(1:ma,ma)
  #clusters=remap[clusters]
  return(clusters)
}

unranked_key <- function(tr) {
  rooted=is.rooted(tr)
  pp <- ape::prop.part(tr)            # list of tip-index vectors, one per internal node
  cl <- lapply(pp, function(x) {
    if (!rooted && 1L %in% x) x <- setdiff(seq_len(n), x)  # orient split away from tip 1
    sort(x)
  })
  len <- lengths(cl)
  cl <- cl[len >= 2L & (rooted | len <= n - 2L)]           # drop trivial splits
  s <- unique(vapply(cl, paste, "", collapse = ","))
  paste(sort(s), collapse = "|")
}

# Canonical integer key for a ranked labeled topology.
ranked_key <- function(tree) {
  n    <- Ntip(tree)
  tree <- reorder(tree, "postorder")          # node ids unchanged
  lab  <- c(match(tree$tip.label, tree$tip.label), rep(Inf, n - 1L))  # min ref-id below each node

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
  #if (any(diff(age[ord]) > -1e-8))
  #  warning("Tied (or nearly tied) coalescence times: ranking is ambiguous.")

  key=as.vector(rbind(pmin(a, b), pmax(a, b))[, ord])
  return(paste(key, collapse = " "))
}

