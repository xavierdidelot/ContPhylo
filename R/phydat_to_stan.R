phyDat_to_stan <- function(dat) {
  stopifnot(inherits(dat, "phyDat"))

  X <- unclass(dat)
  weights <- attr(dat, "weight")

  # phyDat generally stores taxa as list elements.
  taxa <- names(X)
  n <- length(X)
  n_patterns <- length(weights)

  tip_state <- matrix(
    0L,
    nrow = n,
    ncol = n_patterns,
    dimnames = list(taxa, NULL)
  )

  for (i in seq_len(n)) {
    xi <- X[[i]]

    # For standard DNA phyDat objects, unambiguous nucleotides use
    # bit-coded states: A=1, C=2, G=4, T=8.
    tip_state[i, ] <- match(xi, c(1L, 2L, 4L, 8L))
    tip_state[i, is.na(tip_state[i, ])] <- 0L
  }

  list(
    n = n,
    n_patterns = n_patterns,
    tip_state = tip_state,
    pattern_count = as.integer(weights)
  )
}
