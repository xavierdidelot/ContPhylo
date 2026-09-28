#include decode.stan
//this imports build_contphylo which from a and b generates the tree represente as a tuple:
//(edge parent IDs, edge child IDs, edge lengths), matching ape::phylo numbering.

data {
  int<lower=2> n;
  real<lower=0> scale;
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
}

generated quantities {
  tuple(array[2 * n - 2] int,
        array[2 * n - 2] int,
        vector[2 * n - 2]) tree = decode(a, b);

  real log_prior;
  log_prior = exponential_lpdf(tree.3 | 1/scale)+sum(b);
}

