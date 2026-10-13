#' Convert a phylo into a contphylo
#'
#' @param x Object of class `phylo`.
#' @param model Model of tree: 1 for ultrametric, 2 for rooted, 3 for unrooted.
#' @return Object of class `contphylo`.
#' @export
phylo2contphylo <- function(x,model=2) {
  if ((model==1 || model==2) && !is.rooted(x)) stop('Tree needs to be rooted.')
  if (model==1 && !is.ultrametric(x)) stop('Tree needs to be ultrametric.')
  if (model==3 && is.rooted(x)) stop('Tree needs to be unrooted.')
  if (model==3) x=root(x,tail(x$tip.label,1),resolve.root=T)
  edge <- x$edge
  n <- length(x$tip.label)
  parent <- length_to_parent <- rep(NA, max(edge))
  children <- vector("list", max(edge))

  for (k in seq_len(nrow(edge))) {
    p <- edge[k, 1]
    child <- edge[k, 2]
    parent[child] <- p
    length_to_parent[child] <- x$edge.length[k]
    children[[p]] <- c(children[[p]], child)
  }

  root <- setdiff(edge[, 1], edge[, 2])

  tip_key <- function(node) {
    if (node <= n) return(as.numeric(node))
    min(vapply(children[[node]], tip_key, numeric(1)))
  }

  canonicalise <- function(node) {
    if (node <= n) return(invisible())
    keys <- vapply(children[[node]], tip_key, numeric(1))
    children[[node]] <<- children[[node]][order(keys)]
    for (child in children[[node]]) canonicalise(child)
  }
  canonicalise(root)

  tree_length <- sum(x$edge.length)
  coordinates <- matrix(NA_real_, n, 2)
  coordinates[1, ] <- c(0,0)

  contour_offset <- function(target) {
    distance <- 0

    visit <- function(node) {
      for (child in children[[node]]) {
        if (child == target) return(TRUE)
        distance <<- distance + length_to_parent[child]
        if (visit(child)) return(TRUE)
        distance <<- distance + length_to_parent[child]
      }
      FALSE
    }

    visit(root)
    distance
  }

  for (i in n:2) {
    p <- parent[i]
    sibling <- children[[p]][children[[p]] != i]
    pendant_length <- length_to_parent[i]
    graft_distance <- length_to_parent[sibling]
    grandparent <- parent[p]

    subtree_length <- function(node) {
      s <- 0
      for (child in children[[node]]) s <- s + length_to_parent[child] + subtree_length(child)
      s
    }

    if (is.na(grandparent)) {
      tree_length <- tree_length - pendant_length - graft_distance
      coordinates[i, ] <- c(2 * tree_length + graft_distance,
                            pendant_length)
      if (model == 3 && i == n) {
        coordinates[i, ] <- c(NA, pendant_length + graft_distance)
      }
      root <- sibling
      parent[sibling] <- NA
      children[[p]] <- integer()
    } else {
      upper_length <- length_to_parent[p]
      children[[grandparent]][children[[grandparent]] == p] <- sibling
      parent[sibling] <- grandparent
      length_to_parent[sibling] <- graft_distance + upper_length
      children[[p]] <- integer()
      coordinates[i, ] <- c(contour_offset(sibling) + length_to_parent[sibling] +
                              2 * subtree_length(sibling) + graft_distance,
                            pendant_length)
      tree_length <- tree_length - pendant_length
    }
  }
  if (model==1) coordinates=coordinates[,1,drop=F] else coordinates[-1,2] <- log(coordinates[-1,2])
  out=list()
  out$coordinates=coordinates
  out$tip.label=x$tip.label
  out$model=model
  class(out)<-'contphylo'
  return(out)
}
