"""Modal app for parallelised null-model computation.

Dispatches independent random-tree nm computations to Modal cloud.
Use this for large trees (n > 5000) where local computation is slow.

Prerequisites:
    pip install modal
    modal token new  # authenticate

Usage:
    modal run modal_null.py           # run all large trees
    modal run modal_null.py --n 30310 # run specific size
"""

from __future__ import annotations

import modal

app = modal.App("tree-null-model")

image = modal.Image.debian_slim(python_version="3.12").pip_install("numpy")


@app.function(image=image, timeout=7200, retries=1)
def compute_null_nm(n: int, seed: int, model: str = "prufer") -> dict:
    """Generate one random tree of size n and compute its nm.

    Self-contained: includes the independence polynomial algorithm
    so no external files need to be mounted.
    """
    import numpy as np

    rng = np.random.default_rng(seed)

    # --- Generate random tree ---
    if model == "prufer":
        adj = _prufer_tree(n, rng)
    else:
        adj = _yule_tree(n, rng)

    # --- Compute independence polynomial ---
    poly = _independence_poly(n, adj)
    alpha = len(poly) - 1
    mode_idx = poly.index(max(poly))

    # --- Near-miss ratio ---
    nm_val = _near_miss_ratio(poly)

    return {
        "n": n, "seed": seed, "model": model,
        "alpha": alpha, "mode": mode_idx, "nm": nm_val,
    }


def _prufer_tree(n, rng):
    if n <= 2:
        adj = [[] for _ in range(n)]
        if n == 2:
            adj[0].append(1)
            adj[1].append(0)
        return adj

    seq = rng.integers(0, n, size=n - 2)
    degree = [1] * n
    for s in seq:
        degree[s] += 1

    adj = [[] for _ in range(n)]
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

    remaining = [i for i in range(n) if degree[i] == 1]
    if len(remaining) == 2:
        adj[remaining[0]].append(remaining[1])
        adj[remaining[1]].append(remaining[0])
    return adj


def _yule_tree(n, rng):
    if n <= 2:
        adj = [[] for _ in range(n)]
        if n == 2:
            adj[0].append(1)
            adj[1].append(0)
        return adj

    adj = [[]]
    leaves = [0]
    next_id = 1

    while len(adj) < n:
        idx = rng.integers(0, len(leaves))
        parent = leaves[idx]

        if len(adj) + 2 > n:
            child = next_id
            next_id += 1
            adj.append([])
            adj[parent].append(child)
            adj[child].append(parent)
            leaves[idx] = child
            break

        c1, c2 = next_id, next_id + 1
        next_id += 2
        adj.append([])
        adj.append([])
        adj[parent].append(c1)
        adj[c1].append(parent)
        adj[parent].append(c2)
        adj[c2].append(parent)
        leaves[idx] = c1
        leaves.append(c2)

    return adj


def _polyadd(a, b):
    la, lb = len(a), len(b)
    out = [0] * max(la, lb)
    for i in range(la):
        out[i] += a[i]
    for i in range(lb):
        out[i] += b[i]
    return out


def _polymul(a, b):
    if not a or not b:
        return []
    la, lb = len(a), len(b)
    out = [0] * (la + lb - 1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            out[i + j] += ai * bj
    return out


def _independence_poly(n, adj):
    if n == 0:
        return [1]
    if n == 1:
        return [1, 1]

    visited = [False] * n
    components_polys = []

    for start in range(n):
        if visited[start]:
            continue

        comp = []
        q = [start]
        visited[start] = True
        head = 0
        while head < len(q):
            u = q[head]
            head += 1
            comp.append(u)
            for v in adj[u]:
                if not visited[v]:
                    visited[v] = True
                    q.append(v)

        if len(comp) == 1:
            components_polys.append([1, 1])
            continue

        k = len(comp)
        mapping = {old: new for new, old in enumerate(comp)}
        cadj = [[] for _ in range(k)]
        for u in comp:
            for v in adj[u]:
                cadj[mapping[u]].append(mapping[v])

        parent = [-1] * k
        children = [[] for _ in range(k)]
        vis = [False] * k
        vis[0] = True
        bq = [0]
        bh = 0
        while bh < len(bq):
            v = bq[bh]
            bh += 1
            for u in cadj[v]:
                if not vis[u]:
                    vis[u] = True
                    parent[u] = v
                    children[v].append(u)
                    bq.append(u)

        stack = [(0, False)]
        order = []
        while stack:
            v, done = stack.pop()
            if done:
                order.append(v)
                continue
            stack.append((v, True))
            for c in children[v]:
                stack.append((c, False))

        dp0 = [[] for _ in range(k)]
        dp1 = [[] for _ in range(k)]

        for v in order:
            if not children[v]:
                dp0[v] = [1]
                dp1[v] = [0, 1]
            else:
                prod = [1]
                for c in children[v]:
                    prod = _polymul(prod, _polyadd(dp0[c], dp1[c]))
                dp0[v] = prod
                prod = [1]
                for c in children[v]:
                    prod = _polymul(prod, dp0[c])
                dp1[v] = [0] + prod

        components_polys.append(_polyadd(dp0[0], dp1[0]))

    result = components_polys[0]
    for p in components_polys[1:]:
        result = _polymul(result, p)
    return result


def _near_miss_ratio(seq):
    if len(seq) <= 2:
        return 0.0
    first_descent = -1
    for i in range(1, len(seq)):
        if seq[i] < seq[i - 1]:
            first_descent = i
            break
    if first_descent == -1:
        return 0.0
    worst = 0.0
    for j in range(first_descent, len(seq) - 1):
        if seq[j] == 0:
            continue
        ratio = seq[j + 1] / seq[j]
        if ratio > worst:
            worst = ratio
    return worst


@app.local_entrypoint()
def main(n: int = 0, n_samples: int = 10, model: str = "prufer"):
    """Run null-model computation for specified tree size."""
    import json

    if n == 0:
        # Default: run large biological tree sizes
        sizes = [6378, 7649, 10451, 17992, 30310]
        samples_per = {6378: 50, 7649: 50, 10451: 10, 17992: 10, 30310: 10}
    else:
        sizes = [n]
        samples_per = {n: n_samples}

    all_results = {}

    for sz in sizes:
        ns = samples_per.get(sz, n_samples)
        print(f"Dispatching {ns} {model} trees of size {sz}...")

        futures = []
        for i in range(ns):
            seed = 1000 * sz + i
            futures.append(compute_null_nm.spawn(sz, seed, model))

        results = [f.get() for f in futures]
        nms = [r["nm"] for r in results]

        import numpy as np
        mean_nm = np.mean(nms)
        std_nm = np.std(nms)
        print(f"  n={sz}: nm = {mean_nm:.4f} +/- {std_nm:.4f}")

        all_results[f"{sz}_{model}"] = {
            "n": sz, "model": model,
            "nms": nms,
            "mean": float(mean_nm),
            "std": float(std_nm),
            "n_samples": ns,
        }

    out_path = f"null_model_results_modal.json"
    with open(out_path, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\nResults saved to {out_path}")
