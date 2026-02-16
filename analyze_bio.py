"""Compute independence polynomial properties of biological trees.

Usage:
    python analyze_bio.py
"""

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).parent.parent / "Erdos_Problem_993"))

from bio_trees import read_swc, read_newick
from indpoly import independence_poly, is_unimodal, is_log_concave, near_miss_ratio, log_concavity_ratio


def analyze_tree(n, adj, name, timeout=120):
    """Compute and print independence polynomial properties."""
    degrees = [len(adj[i]) for i in range(n)]
    leaves = sum(1 for d in degrees if d <= 1)
    max_deg = max(degrees) if degrees else 0
    d_leaf_counts = []
    for v in range(n):
        leaf_children = sum(1 for u in adj[v] if len(adj[u]) == 1)
        d_leaf_counts.append(leaf_children)
    max_d_leaf = max(d_leaf_counts) if d_leaf_counts else 0
    hubs = sum(1 for c in d_leaf_counts if c >= 2)

    print(f"\n{'='*60}")
    print(f"  {name}")
    print(f"  n={n}, leaves={leaves}, max_degree={max_deg}")
    print(f"  vertices with d_leaf >= 2 (hubs): {hubs}")
    print(f"  max d_leaf: {max_d_leaf}")

    t0 = time.time()
    try:
        poly = independence_poly(n, adj)
    except Exception as e:
        print(f"  FAILED: {e}")
        return None
    elapsed = time.time() - t0

    alpha = len(poly) - 1
    mode_val = max(poly)
    mode_k = poly.index(mode_val)
    uni = is_unimodal(poly)
    lc = is_log_concave(poly)
    nm, nm_pos = near_miss_ratio(poly)
    lcr, lcr_pos = log_concavity_ratio(poly)

    print(f"  independence number alpha = {alpha}")
    print(f"  mode at k = {mode_k} (i_{mode_k} = {mode_val})")
    print(f"  unimodal: {uni}")
    print(f"  log-concave: {lc}")
    if nm_pos >= 0:
        print(f"  near-miss ratio: {nm:.6f} at k={nm_pos}")
    else:
        print(f"  near-miss ratio: n/a (monotone)")
    if lcr_pos >= 0:
        print(f"  log-concavity ratio: {lcr:.6f} at k={lcr_pos}")
    print(f"  time: {elapsed:.2f}s")

    # Print first/last few coefficients for large polynomials
    if alpha <= 20:
        print(f"  polynomial: {poly}")
    else:
        print(f"  polynomial: [{', '.join(str(c) for c in poly[:5])}, ..., {', '.join(str(c) for c in poly[-5:])}]")

    return {
        "name": name, "n": n, "leaves": leaves, "max_deg": max_deg,
        "hubs": hubs, "alpha": alpha, "mode": mode_k,
        "unimodal": uni, "log_concave": lc,
        "nm": nm, "nm_pos": nm_pos, "lcr": lcr, "time": elapsed,
    }


def main():
    base = Path(__file__).parent
    results = []

    # --- Neuromorpho SWC files ---
    swc_dir = base / "data" / "neuromorpho"
    if swc_dir.exists():
        print("\n" + "#" * 60)
        print("# NEUROMORPHO.ORG DENDRITIC ARBORS")
        print("#" * 60)

        swc_files = sorted(swc_dir.glob("*.swc"), key=lambda p: p.stat().st_size)
        for f in swc_files:
            try:
                n, adj, meta = read_swc(f)
            except Exception as e:
                print(f"\n  {f.name}: parse error: {e}")
                continue

            # Skip very large trees (would take too long)
            if n > 5000:
                print(f"\n  {f.name}: n={n}, skipping (too large for initial run)")
                continue

            r = analyze_tree(n, adj, f.name)
            if r:
                r["source"] = "neuromorpho"
                results.append(r)

    # --- Phylogenetic Newick files ---
    nwk_dir = base / "data" / "phylogenies"
    if nwk_dir.exists():
        print("\n" + "#" * 60)
        print("# PHYLOGENETIC TREES (Open Tree of Life)")
        print("#" * 60)

        nwk_files = sorted(nwk_dir.glob("*.nwk"), key=lambda p: p.stat().st_size)
        for f in nwk_files:
            try:
                n, adj, meta = read_newick(f)
            except Exception as e:
                print(f"\n  {f.name}: parse error: {e}")
                continue

            if n > 5000:
                print(f"\n  {f.name}: n={n}, skipping (too large for initial run)")
                continue

            r = analyze_tree(n, adj, f.name)
            if r:
                r["source"] = "phylogeny"
                results.append(r)

    # --- Summary ---
    if results:
        print("\n" + "#" * 60)
        print("# SUMMARY")
        print("#" * 60)
        print(f"\n{'Name':<45} {'n':>6} {'alpha':>6} {'mode':>5} {'uni':>4} {'LC':>4} {'nm':>8} {'time':>7}")
        print("-" * 95)
        for r in results:
            print(f"{r['name']:<45} {r['n']:>6} {r['alpha']:>6} {r['mode']:>5} "
                  f"{'Y' if r['unimodal'] else 'N':>4} "
                  f"{'Y' if r['log_concave'] else 'N':>4} "
                  f"{r['nm']:>8.4f} {r['time']:>6.2f}s")

        all_uni = all(r["unimodal"] for r in results)
        all_lc = all(r["log_concave"] for r in results)
        print(f"\nAll unimodal: {all_uni}")
        print(f"All log-concave: {all_lc}")
        print(f"Max near-miss ratio: {max(r['nm'] for r in results):.6f}")
        print(f"Trees analyzed: {len(results)}")


if __name__ == "__main__":
    main()
