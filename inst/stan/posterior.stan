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
  //real<lower=0> root_distance;
  //real<lower=0> root_sd;
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

array[2 * n - 1] int parent = rep_array(0, 2 * n - 1);
vector[2 * n - 1] length_to_parent =
    rep_vector(0.0, 2 * n - 1);

for (e in 1:(2 * n - 2)) {
  parent[tree.2[e]] = tree.1[e];
  length_to_parent[tree.2[e]] = tree.3[e];
}

/*
//This block of code encourages root position
real root_to_leaf1 = 0;
int node = 1;
for (step in 1:(2 * n - 2)) {
  if (node != n + 1) {
    root_to_leaf1 += length_to_parent[node];
    node = parent[node];
  }
}
target += normal_lpdf(root_to_leaf1 | root_distance, root_sd);
*/

}

