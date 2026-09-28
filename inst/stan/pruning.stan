functions {
  real pruning(
      int n,
      array[,] int tip_state,
      array[] int pattern_count,
      array[] int edge_parent,
      array[] int edge_child,
      vector edge_length) {
    int n_patterns = size(pattern_count);
    int n_nodes = 2 * n - 1;
    int n_edges = 2 * n - 2;

    array[n_nodes, 2] int child_node;
    array[n_nodes, 2] int child_edge;
    array[n_nodes] int n_children;
    array[n_nodes] int scheduled;
    array[n_nodes] int is_child;
    array[n - 1] int postorder;

    int root = 0;
    int n_roots = 0;
    real log_lik = 0;

    if (size(edge_parent) != n_edges ||
        size(edge_child) != n_edges ||
        num_elements(edge_length) != n_edges)
      reject("Tree must contain exactly 2*n-2 edges");

    if (size(tip_state) != n)
      reject("The first dimension of tip_state must equal n");

    for (v in 1:n_nodes) {
      n_children[v] = 0;
      scheduled[v] = v <= n;
      is_child[v] = 0;

      for (k in 1:2) {
        child_node[v, k] = 0;
        child_edge[v, k] = 0;
      }
    }

    // Record the two descendants and their edge indices for each parent.
    for (e in 1:n_edges) {
      int p = edge_parent[e];
      int c = edge_child[e];

      if (p < 1 || p > n_nodes || c < 1 || c > n_nodes)
        reject("Tree contains an invalid node ID");

      n_children[p] += 1;
      if (n_children[p] > 2)
        reject("The tree is not binary at node ", p);

      child_node[p, n_children[p]] = c;
      child_edge[p, n_children[p]] = e;
      is_child[c] = 1;
    }

    // Find the unique internal root.
    for (v in (n + 1):n_nodes) {
      if (is_child[v] == 0) {
        root = v;
        n_roots += 1;
      }
    }

    if (n_roots != 1)
      reject("Expected exactly one internal root, found ", n_roots);

    // Construct a postorder independently of internal-node numbering.
    for (step in 1:(n - 1)) {
      int found = 0;

      for (v in (n + 1):n_nodes) {
        if (found == 0 &&
            scheduled[v] == 0 &&
            n_children[v] == 2 &&
            scheduled[child_node[v, 1]] == 1 &&
            scheduled[child_node[v, 2]] == 1) {
          postorder[step] = v;
          scheduled[v] = 1;
          found = 1;
        }
      }

      if (found == 0)
        reject("Edges do not describe a rooted binary tree");
    }

    // Scaled Felsenstein pruning over compressed site patterns.
    for (s in 1:n_patterns) {
      array[n_nodes] vector[4] partial;
      real log_scale = 0;

      for (tip in 1:n) {
        int state = tip_state[tip, s];

        if (state == 0) {
          partial[tip] = rep_vector(1.0, 4);
        } else {
          partial[tip] = rep_vector(0.0, 4);
          partial[tip][state] = 1.0;
        }
      }

      for (step in 1:(n - 1)) {
        int v = postorder[step];
        int c1 = child_node[v, 1];
        int c2 = child_node[v, 2];
        int e1 = child_edge[v, 1];
        int e2 = child_edge[v, 2];

        real decay1 = exp(-4.0 * edge_length[e1] / 3.0);
        real decay2 = exp(-4.0 * edge_length[e2] / 3.0);
        real off1 = 0.25 * (1.0 - decay1);
        real off2 = 0.25 * (1.0 - decay2);

        vector[4] message1 =
          rep_vector(off1 * sum(partial[c1]), 4) + decay1 * partial[c1];
        vector[4] message2 =
          rep_vector(off2 * sum(partial[c2]), 4) + decay2 * partial[c2];

        real normalizer;
        partial[v] = message1 .* message2;
        normalizer = sum(partial[v]);

        if (normalizer <= 0)
          reject("Non-positive pruning normalizer");

        partial[v] /= normalizer;
        log_scale += log(normalizer);
      }

      // JC69 stationary root frequencies are all 1/4.
      log_lik += pattern_count[s] *
                 (log_scale + log(sum(partial[root])) - log(4.0));
    }

    return log_lik;
  }
}

