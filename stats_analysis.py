"""Statistical analysis: Spearman correlations, confidence intervals, p-values.

Computes:
  1. Spearman correlation of nm vs log(n), pooled and by category
  2. Spearman correlation of nm vs tree-shape statistics (Sackin, Colless-like)
  3. Bootstrap 95% CIs for all correlations

Usage:
    python stats_analysis.py
"""

from __future__ import annotations

import sys
from pathlib import Path
import math

import numpy as np
from scipy.stats import spearmanr

PROJECT = Path(__file__).parent
sys.path.insert(0, str(PROJECT))

from generate_figures import NEURONS, PHYLOGENIES
from bio_trees import read_swc, read_newick
from tree_stats import compute_all_stats


def _bootstrap_spearman_ci(x, y, n_boot=10_000, alpha=0.05, rng=None):
    """Bootstrap 95% CI for Spearman's rho."""
    if rng is None:
        rng = np.random.default_rng(42)
    n = len(x)
    rhos = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, size=n)
        rhos[i] = spearmanr(x[idx], y[idx]).statistic
    lo = np.percentile(rhos, 100 * alpha / 2)
    hi = np.percentile(rhos, 100 * (1 - alpha / 2))
    return lo, hi


def _load_tree_stats():
    """Load all trees and compute shape statistics."""
    results = []

    for name, species, n_val, alpha, mode, nm, fname in NEURONS:
        path = PROJECT / "data" / "neuromorpho" / fname
        if not path.exists():
            continue
        n, adj, meta = read_swc(path)
        stats = compute_all_stats(n, adj, root=meta.get("root", 0))
        results.append({
            "name": name, "category": "neuron", "n": n,
            "alpha": alpha, "mode": mode, "nm": nm,
            **stats,
        })

    for name, n_val, alpha, mode, nm, fname in PHYLOGENIES:
        path = PROJECT / "data" / "phylogenies" / fname
        if not path.exists():
            continue
        n, adj, meta = read_newick(path)
        stats = compute_all_stats(n, adj, root=meta.get("root", 0))
        results.append({
            "name": name, "category": "phylo", "n": n,
            "alpha": alpha, "mode": mode, "nm": nm,
            **stats,
        })

    return results


def _print_correlation(label, x, y, var_names=("x", "y")):
    """Compute and print Spearman correlation with bootstrap CI."""
    rho, p = spearmanr(x, y)
    lo, hi = _bootstrap_spearman_ci(np.array(x), np.array(y))
    sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
    print(f"  {label}: rho = {rho:.4f} [{lo:.4f}, {hi:.4f}], "
          f"p = {p:.2e} {sig}, n = {len(x)}")
    return rho, p, lo, hi


def main():
    print("Loading trees and computing shape statistics...")
    data = _load_tree_stats()
    print(f"  Loaded {len(data)} trees\n")

    neurons = [d for d in data if d["category"] == "neuron"]
    phylos = [d for d in data if d["category"] == "phylo"]

    # --- nm vs log(n) ---
    print("=" * 65)
    print("Spearman correlation: nm vs log(n)")
    print("=" * 65)

    all_logn = [math.log(d["n"]) for d in data]
    all_nm = [d["nm"] for d in data]
    _print_correlation("Pooled (all 27)", all_logn, all_nm)

    neuro_logn = [math.log(d["n"]) for d in neurons]
    neuro_nm = [d["nm"] for d in neurons]
    _print_correlation("Neurons (15)", neuro_logn, neuro_nm)

    phylo_logn = [math.log(d["n"]) for d in phylos]
    phylo_nm = [d["nm"] for d in phylos]
    _print_correlation("Phylogenies (12)", phylo_logn, phylo_nm)

    # --- nm vs tree-shape statistics ---
    print()
    print("=" * 65)
    print("Spearman correlation: nm vs tree-shape statistics")
    print("=" * 65)

    all_sackin = [d["sackin"] for d in data]
    all_colless = [d["colless_like"] for d in data]
    all_leaves = [d["n_leaves"] for d in data]
    all_maxd = [d["max_depth"] for d in data]

    _print_correlation("nm vs Sackin", all_nm, all_sackin)
    _print_correlation("nm vs Colless-like", all_nm, all_colless)
    _print_correlation("nm vs n_leaves", all_nm, all_leaves)
    _print_correlation("nm vs max_depth", all_nm, all_maxd)

    # --- alpha/n ratio ---
    print()
    print("=" * 65)
    print("Independence number ratio alpha/n")
    print("=" * 65)
    for d in data:
        ratio = d["alpha"] / d["n"]
        d["alpha_ratio"] = ratio

    all_ratio = [d["alpha_ratio"] for d in data]
    print(f"  Mean alpha/n: {np.mean(all_ratio):.4f}")
    print(f"  Range: [{min(all_ratio):.4f}, {max(all_ratio):.4f}]")
    neuro_ratio = [d["alpha_ratio"] for d in neurons]
    phylo_ratio = [d["alpha_ratio"] for d in phylos]
    print(f"  Neurons mean: {np.mean(neuro_ratio):.4f}")
    print(f"  Phylogenies mean: {np.mean(phylo_ratio):.4f}")

    # --- mode/alpha ratio ---
    print()
    print("=" * 65)
    print("Mode position ratio mode/alpha")
    print("=" * 65)
    for d in data:
        d["mode_ratio"] = d["mode"] / d["alpha"]
    all_mode_ratio = [d["mode_ratio"] for d in data]
    print(f"  Mean mode/alpha: {np.mean(all_mode_ratio):.4f}")
    print(f"  Range: [{min(all_mode_ratio):.4f}, {max(all_mode_ratio):.4f}]")

    # --- Summary for LaTeX insertion ---
    print()
    print("=" * 65)
    print("LaTeX-ready summary")
    print("=" * 65)
    rho_all, p_all, lo_all, hi_all = spearmanr(all_logn, all_nm)[0], \
        spearmanr(all_logn, all_nm)[1], \
        *_bootstrap_spearman_ci(np.array(all_logn), np.array(all_nm)),
    rho_all = spearmanr(all_logn, all_nm).statistic
    p_all = spearmanr(all_logn, all_nm).pvalue
    lo_all, hi_all = _bootstrap_spearman_ci(np.array(all_logn), np.array(all_nm))

    if p_all < 0.001:
        p_str = "p < 0.001"
    else:
        p_str = f"p = {p_all:.3f}"

    print(f"  Spearman's $\\rho = {rho_all:.3f}$ "
          f"(95\\% CI $[{lo_all:.3f},\\, {hi_all:.3f}]$, "
          f"${p_str}$, $n = {len(data)}$)")


if __name__ == "__main__":
    main()
