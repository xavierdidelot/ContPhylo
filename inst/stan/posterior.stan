functions {
#include decode.stan
//this imports the function build_contphylo which from a and b generates the tree represente as a tuple:
//(edge parent IDs, edge child IDs, edge lengths), matching ape::phylo numbering.
#include pruning.stan
//this imports the function pruning which does the Felsenstein pruning for JC69
}

data {
  int<lower=2> n;
  real<lower=0> scale;
  int<lower=0> n_patterns;
  array[n, n_patterns] int<lower=0, upper=4> tip_state;
  array[n_patterns] int<lower=0> pattern_count;
  // Encoding: 1=A, 2=C, 3=G, 4=T, 0=missing/fully ambiguous.
}

parameters {
  vector[n - 1] a;
  vector[n - 1] b;
}

model {
  tuple(array[2 * n - 2] int,
        array[2 * n - 2] int,
        vector[2 * n - 2]) tree = decode(a, b);

  // Exponential prior on every branch length.
  target += exponential_lpdf(tree.3 | 1/scale);

  // Change-of-variables Jacobian.
  target += sum(b);

  target += (n - 1) * log(2.0) + lgamma(n) - lgamma(2 * n - 1);

  // Sequence likelihood under JC69.
  target += pruning(n, tip_state, pattern_count, tree.1, tree.2, tree.3);
}

