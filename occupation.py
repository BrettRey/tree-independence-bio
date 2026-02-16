"""Per-vertex occupation probabilities P(v) for trees.

Uses a two-pass tree DP (O(n) with big-integer arithmetic):

  Bottom-up (post-order): compute subtree independent-set counts at x=1.
    T0[v] = I(subtree(v); 1) with v excluded
    T1[v] = I(subtree(v); 1) with v included

  Top-down (pre-order): compute rest-of-tree counts.
    up0[v] = # IS of T \\ subtree(v) with parent(v) NOT in the set
    up1[v] = # IS of T \\ subtree(v) with parent(v) IN the set

  Then: P(v) = T1[v] * up0[v] / (T0[root] + T1[root])

Usage:
    from occupation import occupation_probabilities
    probs = occupation_probabilities(n, adj)
"""

from __future__ import annotations

import sys
from collections import deque
from pathlib import Path


def occupation_probabilities(n: int, adj: list[list[int]], root: int = 0) -> list[float]:
    """Compute P(v) for every vertex in a tree.

    Parameters
    ----------
    n : int
        Number of vertices.
    adj : list[list[int]]
        Adjacency list (undirected tree).
    root : int
        Root vertex (default 0).

    Returns
    -------
    list[float]
        P(v) for each vertex, where P(v) = fraction of independent sets
        containing v.
    """
    if n == 0:
        return []
    if n == 1:
        return [0.5]  # 2 IS: {} and {0}; P(0) = 1/2

    # Build rooted tree structure via BFS
    parent = [-1] * n
    children: list[list[int]] = [[] for _ in range(n)]
    visited = [False] * n
    visited[root] = True
    queue = deque([root])
    bfs_order = []
    while queue:
        v = queue.popleft()
        bfs_order.append(v)
        for u in adj[v]:
            if not visited[u]:
                visited[u] = True
                parent[u] = v
                children[v].append(u)
                queue.append(u)

    # Bottom-up pass (post-order): compute T0, T1 using big integers
    T0 = [0] * n
    T1 = [0] * n

    # Process in reverse BFS order (leaves first)
    for v in reversed(bfs_order):
        if not children[v]:
            # Leaf
            T0[v] = 1
            T1[v] = 1
        else:
            prod_sum = 1   # product of (T0[c] + T1[c])
            prod_excl = 1  # product of T0[c]
            for c in children[v]:
                prod_sum *= (T0[c] + T1[c])
                prod_excl *= T0[c]
            T0[v] = prod_sum
            T1[v] = prod_excl

    total = T0[root] + T1[root]  # I(T; 1)

    # Top-down pass (pre-order): compute up0, up1
    up0 = [0] * n
    up1 = [0] * n
    up0[root] = 1
    up1[root] = 0

    for v in bfs_order:
        ch = children[v]
        if not ch:
            continue

        # Precompute prefix and suffix products for sibling exclusion
        k = len(ch)

        # Products of (T0[c] + T1[c]) for siblings
        prefix_sum = [1] * (k + 1)
        suffix_sum = [1] * (k + 1)
        for i in range(k):
            prefix_sum[i + 1] = prefix_sum[i] * (T0[ch[i]] + T1[ch[i]])
        for i in range(k - 1, -1, -1):
            suffix_sum[i] = suffix_sum[i + 1] * (T0[ch[i]] + T1[ch[i]])

        # Products of T0[c] for siblings
        prefix_excl = [1] * (k + 1)
        suffix_excl = [1] * (k + 1)
        for i in range(k):
            prefix_excl[i + 1] = prefix_excl[i] * T0[ch[i]]
        for i in range(k - 1, -1, -1):
            suffix_excl[i] = suffix_excl[i + 1] * T0[ch[i]]

        for i, c in enumerate(ch):
            sib_sum = prefix_sum[i] * suffix_sum[i + 1]
            sib_excl = prefix_excl[i] * suffix_excl[i + 1]
            up0[c] = (up0[v] + up1[v]) * sib_sum
            up1[c] = up0[v] * sib_excl

    # Compute P(v) = T1[v] * up0[v] / total
    probs = [0.0] * n
    for v in range(n):
        probs[v] = (T1[v] * up0[v]) / total

    return probs


def _verify():
    """Spot-check against hand-computed values."""
    print("Verifying occupation probabilities...")

    # P3: path 0-1-2, rooted at 1
    # IS: {}, {0}, {1}, {2}, {0,2} => I(T;1) = 5
    # P(0) = 2/5 ({0}, {0,2}), P(1) = 1/5 ({1}), P(2) = 2/5 ({2}, {0,2})
    adj_p3 = [[1], [0, 2], [1]]
    probs = occupation_probabilities(3, adj_p3, root=1)
    assert abs(probs[0] - 2 / 5) < 1e-12, f"P3: P(0) = {probs[0]}, expected 0.4"
    assert abs(probs[1] - 1 / 5) < 1e-12, f"P3: P(1) = {probs[1]}, expected 0.2"
    assert abs(probs[2] - 2 / 5) < 1e-12, f"P3: P(2) = {probs[2]}, expected 0.4"
    print(f"  P3 (root=1): P = [{probs[0]:.4f}, {probs[1]:.4f}, {probs[2]:.4f}]  OK")

    # K_{1,5}: star with hub 0, leaves 1..5
    # IS: {}, {0}, {1}, {2}, {3}, {4}, {5}, and all subsets of {1..5}
    # Total = 2^5 + 1 = 33
    # P(0) = 1/33 (only {0} includes hub)
    # P(leaf) = 2^4 / 33 = 16/33 (any subset of other 4 leaves, plus this leaf)
    adj_star = [[1, 2, 3, 4, 5], [0], [0], [0], [0], [0]]
    probs = occupation_probabilities(6, adj_star, root=0)
    assert abs(probs[0] - 1 / 33) < 1e-12, f"K_1,5: P(hub) = {probs[0]}, expected {1/33}"
    for i in range(1, 6):
        assert abs(probs[i] - 16 / 33) < 1e-12, f"K_1,5: P({i}) = {probs[i]}, expected {16/33}"
    print(f"  K_{{1,5}}: P(hub) = {probs[0]:.4f}, P(leaf) = {probs[1]:.4f}  OK")

    # P5: path 0-1-2-3-4
    # IS of P5: {}, {0},{1},{2},{3},{4}, {0,2},{0,3},{0,4},{1,3},{1,4},{2,4},
    #           {0,2,4},{0,3,4} wait...
    # Actually I(P5;1) = F(7) = 13 (Fibonacci)
    # Let me just check sum of P(v) = mean IS size
    adj_p5 = [[1], [0, 2], [1, 3], [2, 4], [3]]
    probs = occupation_probabilities(5, adj_p5, root=0)
    # Sum of P(v) should equal (sum of sizes of all IS) / (number of IS)
    # = (0+1+1+1+1+1+2+2+2+2+2+2+3+3) / 13 ... let me just verify sum
    # IS: {}, {0},{1},{2},{3},{4},
    #     {0,2},{0,3},{0,4},{1,3},{1,4},{2,4},
    #     {0,2,4},{1,4} wait that's wrong. {1,4} has distance 3, so yes independent.
    #     3-element: {0,2,4}, {0,3,... 3 adj to 2 and 4, so {0,3} need 3 not adj to 0.
    #     3 is not adj to 0, so {0,3} is IS. {0,3,?} ? must not be adj to 0 or 3.
    #     Not adj to 0: {2,3,4}; not adj to 3: {0,1}. Intersection empty besides already.
    #     So only 3-elem IS is {0,2,4}.
    # Total: 1 + 5 + 6 + 1 = 13. Correct!
    # Sum of sizes: 0 + 5*1 + 6*2 + 1*3 = 20
    # Mean IS size = 20/13
    # Sum of P(v) should = 20/13
    psum = sum(probs)
    assert abs(psum - 20 / 13) < 1e-12, f"P5: sum(P) = {psum}, expected {20/13}"
    print(f"  P5: sum(P(v)) = {psum:.4f}, expected {20/13:.4f}  OK")

    # Check that P(v) sums correctly for all trees above
    # For K_{1,5}: mean IS size = (0 + 1 + 5*16) / 33 = 81/33 = 27/11
    psum_star = sum(probs_s for probs_s in occupation_probabilities(6, adj_star, root=0))
    # Total size sum: 0*(1) + 1*(1 hub set) + 1*5(singleton leaves) + 2*C(5,2) + 3*C(5,3) + 4*C(5,4) + 5*C(5,5)
    # = 0 + 1 + 5 + 20 + 30 + 20 + 5 = 81
    # Divided by 33 = 81/33
    assert abs(psum_star - 81 / 33) < 1e-12, f"K_1,5: sum(P) = {psum_star}, expected {81/33}"
    print(f"  K_{{1,5}}: sum(P(v)) = {psum_star:.4f}, expected {81/33:.4f}  OK")

    print("All verifications passed.")


if __name__ == "__main__":
    _verify()

    # If bio_trees is available, compute P(v) for a sample tree
    try:
        PROJECT = Path(__file__).parent
        sys.path.insert(0, str(PROJECT))
        from bio_trees import read_swc, read_newick

        # Small example: Homininae
        path = PROJECT / "data" / "phylogenies" / "homininae.nwk"
        if path.exists():
            n, adj, meta = read_newick(path)
            probs = occupation_probabilities(n, adj, root=meta.get("root", 0))
            degrees = [len(adj[i]) for i in range(n)]
            print(f"\nHomininae (n={n}):")
            print(f"  Mean P(v) = {sum(probs)/n:.4f}")
            print(f"  Max P(v)  = {max(probs):.4f} (degree {degrees[probs.index(max(probs))]})")
            print(f"  Min P(v)  = {min(probs):.4f} (degree {degrees[probs.index(min(probs))]})")
    except ImportError:
        pass
