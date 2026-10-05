# phyDat_to_stan <- function(dat, compress_jc69 = TRUE) {
#   stopifnot(inherits(dat, "phyDat"))
#
#   X <- unclass(dat)
#   weights <- as.integer(attr(dat, "weight"))
#   taxa <- names(X)
#   n <- length(X)
#   n_patterns <- length(weights)
#
#   # phyDat entries are indices into allLevels.
#   all_levels <- tolower(attr(dat, "allLevels"))
#   state_lookup <- match(all_levels, c("a", "c", "g", "t"))
#
#   tip_state <- matrix(
#     0L,
#     nrow = n,
#     ncol = n_patterns,
#     dimnames = list(taxa, NULL)
#   )
#
#   for (i in seq_len(n)) {
#     states <- state_lookup[X[[i]]]
#     states[is.na(states)] <- 0L
#     tip_state[i, ] <- states
#   }
#
#   if (compress_jc69 && n_patterns > 0L) {
#     # Canonicalize each column under a global permutation of A/C/G/T.
#     # The first observed state becomes 1, the next new state 2, etc.;
#     # missing/fully ambiguous state 0 remains 0.
#     for (s in seq_len(n_patterns)) {
#       map <- integer(4L)
#       next_state <- 0L
#
#       for (i in seq_len(n)) {
#         old_state <- tip_state[i, s]
#
#         if (old_state != 0L) {
#           if (map[old_state] == 0L) {
#             next_state <- next_state + 1L
#             map[old_state] <- next_state
#           }
#           tip_state[i, s] <- map[old_state]
#         }
#       }
#     }
#
#     keys <- apply(tip_state, 2L, paste0, collapse = "")
#     keep <- !duplicated(keys)
#     group <- match(keys, keys[keep])
#     n_groups <- sum(keep)
#
#     weights <- as.integer(
#       rowsum(matrix(weights, ncol = 1L), group, reorder = FALSE)[, 1L]
#     )
#     tip_state <- tip_state[, keep, drop = FALSE]
#     n_patterns <- n_groups
#   }
#
#   list(
#     n = n,
#     n_patterns = n_patterns,
#     tip_state = tip_state,
#     pattern_count = weights
#   )
# }
