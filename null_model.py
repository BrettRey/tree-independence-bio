"""Null-model comparison: random trees matched to biological tree sizes.

Generates random trees (Prüfer-sequence uniform and Yule birth-only),
computes their independence polynomials, and compares nm to biological trees.

Sample sizes (per biological tree, per model):
  n < 500:   1000 random trees
  n < 5000:  200 random trees
  n < 10000: 50 random trees
  n >= 10000: 10 random trees

Usage:
    python null_model.py                  # run all (may take hours)
    python null_model.py --sizes 46 113   # run specific sizes
    python null_model.py --quick          # 10 samples per size (testing)
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

PROJECT = Path(__file__).parent
sys.path.insert(0, str(PROJECT))

from indpoly import independence_poly, near_miss_ratio


# ---------------------------------------------------------------------------
# Random tree generators
# ---------------------------------------------------------------------------

def prufer_random_tree(n: int, rng: np.random.Generator) -> list[list[int]]:
    """Generate a uniform random labelled tree via Prüfer sequence.

    Returns adjacency list for n vertices (0-indexed).
    """
    if n == 1:
        return [[]]
    if n == 2:
        return [[1], [0]]

    seq = rng.integers(0, n, size=n - 2)

    degree = [1] * n
    for s in seq:
        degree[s] += 1

    adj: list[list[int]] = [[] for _ in range(n)]

    ptr = 0
    while degree[ptr] != 1:
        ptr += 1
    leaf = ptr

    for s in seq:
        adj[leaf].append(s)
        adj[s].append(leaf)
        degree[s] -= 1
        if degree[s] == 1 and s < ptr:
            leaf = s
        else:
            ptr += 1
            while ptr < n and degree[ptr] != 1:
                ptr += 1
            leaf = ptr

    # Last edge: two vertices with degree 1
    remaining = [i for i in range(n) if degree[i] == 1]
    if len(remaining) == 2:
        adj[remaining[0]].append(remaining[1])
        adj[remaining[1]].append(remaining[0])

    return adj


def yule_random_tree(n: int, rng: np.random.Generator) -> list[list[int]]:
    """Generate a Yule (birth-only) random tree with n vertices.

    Starts with a root (vertex 0). At each step, pick a random leaf and
    give it two children (replacing the leaf with an internal node).
    Continues until we have at least n vertices.

    The resulting tree has (n-1)//2 internal nodes and (n+1)//2 leaves
    when n is odd (we stop at n). For even n, we stop at n+1 and
    truncate.
    """
    if n <= 2:
        adj: list[list[int]] = [[] for _ in range(n)]
        if n == 2:
            adj[0].append(1)
            adj[1].append(0)
        return adj

    # Start with single root
    adj = [[]]  # vertex 0
    leaves = [0]
    next_id = 1

    while len(adj) < n:
        # Pick a random leaf to split
        idx = rng.integers(0, len(leaves))
        parent = leaves[idx]

        # Need at least 2 more vertices
        if len(adj) + 2 > n:
            # Add just one child to reach n exactly
            child = next_id
            next_id += 1
            adj.append([])
            adj[parent].append(child)
            adj[child].append(parent)
            # parent stays a leaf only if it had no other children
            # Actually it now has a child, so it's internal
            leaves[idx] = child
            break

        # Add two children
        c1, c2 = next_id, next_id + 1
        next_id += 2
        adj.append([])
        adj.append([])
        adj[parent].append(c1)
        adj[c1].append(parent)
        adj[parent].append(c2)
        adj[c2].append(parent)

        # parent is no longer a leaf; replace with c1, append c2
        leaves[idx] = c1
        leaves.append(c2)

    return adj


# ---------------------------------------------------------------------------
# Null-model computation
# ---------------------------------------------------------------------------

def compute_null_nm(n: int, model: str, n_samples: int,
                    seed: int = 42) -> list[float]:
    """Generate n_samples random trees of size n and compute their nm values."""
    rng = np.random.default_rng(seed)
    gen = prufer_random_tree if model == "prufer" else yule_random_tree

    nms = []
    for i in range(n_samples):
        adj = gen(n, rng)
        poly = independence_poly(n, adj)
        nm_val, _ = near_miss_ratio(poly)
        nms.append(nm_val)

    return nms


def default_n_samples(n: int) -> int:
    """Default number of random trees for a given tree size."""
    if n < 500:
        return 1000
    elif n < 5000:
        return 200
    elif n < 10000:
        return 50
    else:
        return 10


def run_null_model(sizes: list[int] | None = None,
                   quick: bool = False,
                   models: list[str] | None = None) -> dict:
    """Run null-model comparison for all biological tree sizes.

    Returns dict mapping (n, model) -> list of nm values.
    """
    from generate_figures import NEURONS, PHYLOGENIES

    if models is None:
        models = ["prufer", "yule"]

    if sizes is None:
        bio_sizes = sorted(set(
            [row[2] for row in NEURONS] + [row[1] for row in PHYLOGENIES]
        ))
    else:
        bio_sizes = sorted(sizes)

    # Build lookup for biological nm
    bio_nm = {}
    for row in NEURONS:
        bio_nm[row[2]] = row[5]
    for row in PHYLOGENIES:
        bio_nm[row[1]] = row[4]

    results = {}
    total_trees = 0

    for n in bio_sizes:
        n_samp = 10 if quick else default_n_samples(n)
        total_trees += n_samp * len(models)

    print(f"Null model: {len(bio_sizes)} sizes, {len(models)} models, "
          f"~{total_trees} total random trees")
    print()

    for n in bio_sizes:
        n_samp = 10 if quick else default_n_samples(n)
        bio = bio_nm.get(n, None)

        for model in models:
            t0 = time.time()
            nms = compute_null_nm(n, model, n_samp, seed=42 + n)
            elapsed = time.time() - t0

            mean_nm = np.mean(nms)
            std_nm = np.std(nms)

            # How many standard deviations is the biological nm from the null?
            if bio is not None and std_nm > 0:
                z = (bio - mean_nm) / std_nm
                z_str = f"z = {z:+.2f}"
            else:
                z_str = "---"

            print(f"  n={n:>6}, {model:>6}: {n_samp:>4} trees, "
                  f"{elapsed:>7.1f}s, "
                  f"null nm = {mean_nm:.4f} +/- {std_nm:.4f}, "
                  f"bio nm = {bio:.4f}, {z_str}")

            results[(n, model)] = {
                "nms": nms,
                "mean": float(mean_nm),
                "std": float(std_nm),
                "bio_nm": bio,
                "n_samples": n_samp,
            }

    return results


def save_results(results: dict, path: Path):
    """Save null-model results to JSON."""
    # Convert tuple keys to strings for JSON
    out = {}
    for (n, model), data in results.items():
        key = f"{n}_{model}"
        out[key] = {
            "n": n, "model": model,
            "nms": data["nms"],
            "mean": data["mean"],
            "std": data["std"],
            "bio_nm": data["bio_nm"],
            "n_samples": data["n_samples"],
        }
    path.write_text(json.dumps(out, indent=2))
    print(f"\nResults saved to {path}")


def main():
    parser = argparse.ArgumentParser(description="Null-model nm comparison")
    parser.add_argument("--sizes", nargs="+", type=int,
                        help="Specific tree sizes to test")
    parser.add_argument("--quick", action="store_true",
                        help="Quick mode: 10 samples per size")
    parser.add_argument("--model", choices=["prufer", "yule", "both"],
                        default="both", help="Which random tree model(s)")
    parser.add_argument("--output", type=str, default=None,
                        help="Output JSON path")
    args = parser.parse_args()

    models = ["prufer", "yule"] if args.model == "both" else [args.model]
    results = run_null_model(sizes=args.sizes, quick=args.quick, models=models)

    out_path = Path(args.output) if args.output else PROJECT / "null_model_results.json"
    save_results(results, out_path)

    # Summary statistics
    print("\n" + "=" * 65)
    print("Summary: biological nm vs null-model nm")
    print("=" * 65)
    for model in models:
        entries = [(n, d) for (n, m), d in results.items() if m == model]
        entries.sort()
        bio_above = sum(1 for n, d in entries
                        if d["bio_nm"] is not None and d["bio_nm"] > d["mean"])
        bio_below = sum(1 for n, d in entries
                        if d["bio_nm"] is not None and d["bio_nm"] <= d["mean"])
        print(f"  {model}: bio nm > null mean in {bio_above}/{bio_above+bio_below} cases")


if __name__ == "__main__":
    main()
