#' Convert a phylo into a contphylo
#'
#' @param x Object of class `phylo`.
#' @return Object of class `contphylo`.
#' @export
phylo2contphylo <- function(x) {
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
  tree_length <- sum(x$edge.length)
  coordinates <- matrix(NA_real_, n, 2)
  coordinates[1, ] <- c(1,1)

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

    if (is.na(grandparent)) {
      tree_length <- tree_length - pendant_length - graft_distance
      coordinates[i, ] <- c(2 * tree_length + graft_distance,
                            pendant_length)
      root <- sibling
      parent[sibling] <- NA
      children[[p]] <- integer()
    } else {
      upper_length <- length_to_parent[p]
      children[[grandparent]][children[[grandparent]] == p] <- sibling
      parent[sibling] <- grandparent
      length_to_parent[sibling] <- graft_distance + upper_length
      children[[p]] <- integer()
      coordinates[i, ] <- c(contour_offset(sibling) + upper_length,
                            pendant_length)
      tree_length <- tree_length - pendant_length
    }
  }

  coordinates[,2] <- log(coordinates[,2])
  rownames(coordinates) <- x$tip.label
  colnames(coordinates) <- c("contour", "pendant")

  out <- structure(list(coordinates = coordinates), class = "contphylo")
  out$tip.label <- x$tip.label
  out
}
