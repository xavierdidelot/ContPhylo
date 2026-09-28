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

    vector[n_edges] decay;
    vector[n_edges] off;

    int root = 0;
    int n_roots = 0;
    real log_lik = 0;

    if (size(edge_parent) != n_edges ||
        size(edge_child) != n_edges ||
        num_elements(edge_length) != n_edges)
      reject("Tree must contain exactly 2*n-2 edges");

    if (size(tip_state) != n)
      reject("The first dimension of tip_state must equal n");

    // Compute each branch's JC69 quantities only once per log-density
    // evaluation, rather than once per site pattern.
    decay = exp((-4.0 / 3.0) * edge_length);
    off = 0.25 * (rep_vector(1.0, n_edges) - decay);

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
    if (n_patterns > 0) {
      for (s in 1:n_patterns) {
        array[n_nodes] vector[4] partial;
        real log_scale = 0;

        for (step in 1:(n - 1)) {
          int v = postorder[step];
          int c1 = child_node[v, 1];
          int c2 = child_node[v, 2];
          int e1 = child_edge[v, 1];
          int e2 = child_edge[v, 2];

          vector[4] message1;
          vector[4] message2;
          real normalizer;

          // For a tip, form its branch message directly without allocating
          // and multiplying a separate tip partial vector.
          if (c1 <= n) {
            int state1 = tip_state[c1, s];

            if (state1 == 0) {
              message1 = rep_vector(1.0, 4);
            } else {
              message1 = rep_vector(off[e1], 4);
              message1[state1] += decay[e1];
            }
          } else {
            // Internal partials are normalized, so their elements sum to 1.
            message1 = rep_vector(off[e1], 4) + decay[e1] * partial[c1];
          }

          if (c2 <= n) {
            int state2 = tip_state[c2, s];

            if (state2 == 0) {
              message2 = rep_vector(1.0, 4);
            } else {
              message2 = rep_vector(off[e2], 4);
              message2[state2] += decay[e2];
            }
          } else {
            message2 = rep_vector(off[e2], 4) + decay[e2] * partial[c2];
          }

          partial[v] = message1 .* message2;
          normalizer = sum(partial[v]);

          if (normalizer <= 0)
            reject("Non-positive pruning normalizer");

          partial[v] /= normalizer;
          log_scale += log(normalizer);
        }

        // The normalized root partial sums to one. JC69 root frequencies
        // are all 1/4, so the root contribution is simply -log(4).
        log_lik += pattern_count[s] * (log_scale - log(4.0));
      }
    }

    return log_lik;
  }

