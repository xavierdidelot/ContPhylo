functions {
  // Return (graft node, distance from graft node to the new attachment point).
  tuple(int, real) contour_traversal(array[] int child1,
                                     array[] int child2,
                                     array[] int parent,
                                     vector dist_to_parent,
                                     real d) {
    int N = size(parent);
    int root = 0;
    int from;
    int cur;
    real remainder = d;

    for (v in 1:N) {
      if (child1[v] != 0 && parent[v] == 0) {
        root = v;
        break;
      }
    }
    if (root == 0)
      reject("contour_traversal: could not find root");

    if (d < 0)
      return (root, -d);

    from = root;
    cur = child1[root];

    // A complete contour uses O(N) state changes; the guard catches bugs.
    for (step in 1:(20 * N + 20)) {
      if (cur == root && from == child2[root])
        return (cur, remainder);

      if (cur == root && from == child1[root]) {
        from = root;
        cur = child2[root];
        continue;
      }

      if (remainder < dist_to_parent[cur] && from == parent[cur]) {
        remainder = dist_to_parent[cur] - remainder;
        return (cur, remainder);
      }

      if (child1[cur] == 0 && remainder < 2 * dist_to_parent[cur]) {
        remainder -= dist_to_parent[cur];
        return (cur, remainder);
      }

      if (child1[cur] == 0) {
        remainder -= 2 * dist_to_parent[cur];
        from = cur;
        cur = parent[cur];
        continue;
      }

      if (from == parent[cur]) {
        remainder -= dist_to_parent[cur];
        from = cur;
        cur = child1[cur];
        continue;
      }

      if (from == child1[cur]) {
        from = cur;
        cur = child2[cur];
        continue;
      }

      if (from == child2[cur]) {
        if (remainder < dist_to_parent[cur])
          return (cur, remainder);
        remainder -= dist_to_parent[cur];
        from = cur;
        cur = parent[cur];
        continue;
      }

      reject("contour_traversal: invalid traversal state");
    }

    reject("contour_traversal: iteration guard exceeded");
    //return (root, remainder); // Unreachable, but required by the type checker.
  }

  // Implements the final R relabelling that swaps the actual root with n + 1.
  int relabel_node(int node, int root, int conventional_root) {
    if (node == root)
      return conventional_root;
    if (node == conventional_root)
      return root;
    return node;
  }

  // Input convention:
  //   a[i - 1] = coords[i, 1], b[i - 1] = coords[i, 2], for R leaf i = 2,...,n.
  // Output tuple:
  //   (edge parent IDs, edge child IDs, edge lengths), matching phylo numbering.
  tuple(array[] int, array[] int, vector) decode(vector a,vector b) {
    int n = num_elements(a) + 1;
    int N = 2 * n - 1;
    int E = 2 * n - 2;
    int conventional_root = n + 1;
    int root = 0;
    int edge_index = 0;
    array[N] int child1 = rep_array(0, N);
    array[N] int child2 = rep_array(0, N);
    array[N] int parent = rep_array(0, N);
    vector[N] dist_to_parent = rep_vector(0.0, N);
    array[E] int edge_parent;
    array[E] int edge_child;
    vector[E] edge_length;

    if (num_elements(b) != n - 1)
      reject("build_contphylo: a and b must have equal lengths");

    for (i in 2:n) {
      int graft;
      int old_parent;
      int new_node = n + i - 1;
      real graft_dist;

      if (i == 2) {
        graft = 1;
        graft_dist = abs(a[i - 1]);
      } else {
        tuple(int, real) location = contour_traversal(
          child1, child2, parent, dist_to_parent, a[i - 1]
        );
        graft = location.1;
        graft_dist = location.2;
      }

      old_parent = parent[graft];

      parent[i] = new_node;
      dist_to_parent[i] = exp(b[i - 1]);

      child1[new_node] = graft;
      child2[new_node] = i;
      parent[new_node] = old_parent;
      dist_to_parent[new_node] = dist_to_parent[graft] - graft_dist;

      if (old_parent != 0) {
        if (child1[old_parent] == graft)
          child1[old_parent] = new_node;
        if (child2[old_parent] == graft)
          child2[old_parent] = new_node;
      }

      parent[graft] = new_node;
      dist_to_parent[graft] = graft_dist;
    }

    for (v in 1:N) {
      if (child1[v] != 0 && parent[v] == 0) {
        root = v;
        break;
      }
    }
    if (root == 0)
      reject("build_contphylo: could not find final root");

    // Emit rows in the same child-label order as the R phylo object after its
    // root swap: 1,...,n,n+2,...,2n-1. The root n+1 is omitted.
    for (new_child in 1:N) {
      if (new_child != conventional_root) {
        int old_child = relabel_node(new_child, root, conventional_root);
        int old_parent = parent[old_child];

        if (old_parent == 0)
          reject("build_contphylo: non-root node has no parent");

        edge_index += 1;
        edge_child[edge_index] = new_child;
        edge_parent[edge_index] = relabel_node(
          old_parent, root, conventional_root
        );
        edge_length[edge_index] = dist_to_parent[old_child];
      }
    }

    return (edge_parent, edge_child, edge_length);
  }
}

