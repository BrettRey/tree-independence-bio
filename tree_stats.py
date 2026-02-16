"""Tree-shape statistics: Sackin index, Colless-like index, subtree sizes.

Standard phylogenetic balance measures computed on rooted trees.
All functions take (n, adj, root) and run in O(n).

Usage:
    from tree_stats import sackin_index, colless_like_index, subtree_sizes
"""

from __future__ import annotations

import sys
from collections import deque
from pathlib import Path


def _rooted_structure(n: int, adj: list[list[int]], root: int = 0):
    """BFS to get parent, children, and BFS order for a rooted tree."""
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
    return parent, children, bfs_order


def subtree_sizes(n: int, adj: list[list[int]], root: int = 0) -> list[int]:
    """Compute the subtree size of every vertex.

    Returns a list where sizes[v] = number of vertices in the subtree rooted at v.
    """
    _, children, bfs_order = _rooted_structure(n, adj, root)
    sizes = [1] * n
    for v in reversed(bfs_order):
        for c in children[v]:
            sizes[v] += sizes[c]
    return sizes


def sackin_index(n: int, adj: list[list[int]], root: int = 0) -> int:
    """Compute the Sackin index: sum of leaf depths.

    Measures tree imbalance. Higher values = more imbalanced.
    """
    _, children, bfs_order = _rooted_structure(n, adj, root)
    depth = [0] * n
    for v in bfs_order:
        for c in children[v]:
            depth[c] = depth[v] + 1

    return sum(depth[v] for v in range(n) if not children[v])


def colless_like_index(n: int, adj: list[list[int]], root: int = 0) -> int:
    """Compute the Colless-like index for possibly multifurcating trees.

    For each internal node, adds (max child subtree size - min child subtree size).
    This generalises the standard Colless index (which uses |n_left - n_right|
    for bifurcating trees) to multifurcating trees.
    """
    _, children, bfs_order = _rooted_structure(n, adj, root)
    sizes = [1] * n
    for v in reversed(bfs_order):
        for c in children[v]:
            sizes[v] += sizes[c]

    total = 0
    for v in range(n):
        if children[v]:
            child_sizes = [sizes[c] for c in children[v]]
            total += max(child_sizes) - min(child_sizes)
    return total


def compute_all_stats(n: int, adj: list[list[int]], root: int = 0) -> dict:
    """Compute all tree-shape statistics at once (single traversal)."""
    _, children, bfs_order = _rooted_structure(n, adj, root)

    # Depths (for Sackin)
    depth = [0] * n
    for v in bfs_order:
        for c in children[v]:
            depth[c] = depth[v] + 1

    # Subtree sizes (for Colless-like)
    sizes = [1] * n
    for v in reversed(bfs_order):
        for c in children[v]:
            sizes[v] += sizes[c]

    sackin = sum(depth[v] for v in range(n) if not children[v])
    n_leaves = sum(1 for v in range(n) if not children[v])

    colless = 0
    for v in range(n):
        if children[v]:
            child_sizes = [sizes[c] for c in children[v]]
            colless += max(child_sizes) - min(child_sizes)

    max_depth = max(depth) if n > 0 else 0

    return {
        "sackin": sackin,
        "colless_like": colless,
        "n_leaves": n_leaves,
        "max_depth": max_depth,
        "mean_leaf_depth": sackin / n_leaves if n_leaves > 0 else 0,
    }


def _verify():
    """Spot-check on small trees."""
    print("Verifying tree-shape statistics...")

    # Path P4: 0-1-2-3, rooted at 0
    # Depths: 0,1,2,3; leaves = {3}; Sackin = 3
    # Wait, P4 has internal nodes 0,1,2 each with 1 child, leaf = {3}
    # Actually in a path, only the endpoint without children is a leaf.
    # With root=0: children[0]=[1], children[1]=[2], children[2]=[3], children[3]=[]
    # Leaf = {3}, Sackin = 3
    adj_p4 = [[1], [0, 2], [1, 3], [2]]
    s = sackin_index(4, adj_p4, root=0)
    assert s == 3, f"P4: Sackin = {s}, expected 3"
    print(f"  P4 (root=0): Sackin = {s}  OK")

    # Star K_{1,4}: hub 0, leaves 1,2,3,4
    # All leaves at depth 1, Sackin = 4
    adj_star = [[1, 2, 3, 4], [0], [0], [0], [0]]
    s = sackin_index(5, adj_star, root=0)
    assert s == 4, f"K_1,4: Sackin = {s}, expected 4"
    print(f"  K_{{1,4}}: Sackin = {s}  OK")

    # Colless-like for star: root has 4 children all with subtree size 1
    # max - min = 1 - 1 = 0
    c = colless_like_index(5, adj_star, root=0)
    assert c == 0, f"K_1,4: Colless = {c}, expected 0"
    print(f"  K_{{1,4}}: Colless-like = {c}  OK")

    # Caterpillar: 0-1-2, with 1 having extra leaf child 3
    # adj: 0:[1], 1:[0,2,3], 2:[1], 3:[1]
    # Root at 0: children[0]=[1], children[1]=[2,3], children[2]=[], children[3]=[]
    # Subtree sizes: 0->4, 1->3, 2->1, 3->1
    # Leaves: {2,3}, depths: 2,2. Sackin = 4
    # Colless: node 0 has child 1 with size 3, max-min = 3-3=0 (only 1 child)
    #   node 1 has children 2,3 with sizes 1,1: max-min=0
    # Total Colless = 0
    adj_cat = [[1], [0, 2, 3], [1], [1]]
    s = sackin_index(4, adj_cat, root=0)
    assert s == 4, f"Caterpillar: Sackin = {s}, expected 4"
    c = colless_like_index(4, adj_cat, root=0)
    assert c == 0, f"Caterpillar: Colless = {c}, expected 0"
    print(f"  Caterpillar: Sackin = {s}, Colless-like = {c}  OK")

    # Unbalanced tree for non-zero Colless:
    #      0
    #     / \
    #    1   4
    #   / \
    #  2   3
    # Root 0: children [1,4]
    # Subtree sizes: 0->5, 1->3, 2->1, 3->1, 4->1
    # Node 0: children have sizes 3,1 -> 3-1=2
    # Node 1: children have sizes 1,1 -> 0
    # Colless = 2
    adj_unbal = [[1, 4], [0, 2, 3], [1], [1], [0]]
    c = colless_like_index(5, adj_unbal, root=0)
    assert c == 2, f"Unbalanced: Colless = {c}, expected 2"
    s = sackin_index(5, adj_unbal, root=0)
    # Leaves: 2 (depth 2), 3 (depth 2), 4 (depth 1). Sackin = 5
    assert s == 5, f"Unbalanced: Sackin = {s}, expected 5"
    print(f"  Unbalanced tree: Sackin = {s}, Colless-like = {c}  OK")

    print("All verifications passed.")


if __name__ == "__main__":
    _verify()

    # Compute stats for all biological trees if available
    try:
        PROJECT = Path(__file__).parent
        sys.path.insert(0, str(PROJECT))
        from bio_trees import read_swc, read_newick
        from generate_figures import NEURONS, PHYLOGENIES

        print("\nTree-shape statistics for biological trees:")
        print(f"{'Tree':<35} {'n':>6} {'Sackin':>8} {'Colless':>8} {'Leaves':>7} {'MaxD':>5}")
        print("-" * 75)

        for name, species, n_val, alpha, mode, nm, fname in NEURONS:
            path = PROJECT / "data" / "neuromorpho" / fname
            if path.exists():
                n, adj, meta = read_swc(path)
                stats = compute_all_stats(n, adj, root=meta.get("root", 0))
                print(f"{name:<35} {n:>6} {stats['sackin']:>8} "
                      f"{stats['colless_like']:>8} {stats['n_leaves']:>7} "
                      f"{stats['max_depth']:>5}")

        for name, n_val, alpha, mode, nm, fname in PHYLOGENIES:
            path = PROJECT / "data" / "phylogenies" / fname
            if path.exists():
                n, adj, meta = read_newick(path)
                stats = compute_all_stats(n, adj, root=meta.get("root", 0))
                print(f"{name:<35} {n:>6} {stats['sackin']:>8} "
                      f"{stats['colless_like']:>8} {stats['n_leaves']:>7} "
                      f"{stats['max_depth']:>5}")
    except (ImportError, Exception) as e:
        print(f"Could not compute biological tree stats: {e}")
