"""Convert biological tree formats to (n, adj) for indpoly.py.

Supported formats:
  - SWC (NeuroMorpho.org neuronal reconstructions)
  - Newick (phylogenetic trees from TreeHub, Open Tree of Life, etc.)

Usage:
    from bio_trees import read_swc, read_newick

    n, adj, meta = read_swc("neuron.swc")
    n, adj, meta = read_newick("phylo.nwk")

    # Then pass to indpoly:
    from indpoly import independence_poly, is_unimodal, near_miss_ratio
    poly = independence_poly(n, adj)
"""

from __future__ import annotations

import re
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# SWC reader
# ---------------------------------------------------------------------------

def read_swc(path: str | Path) -> tuple[int, list[list[int]], dict]:
    """Read an SWC file and return (n, adj, metadata).

    Parameters
    ----------
    path : str or Path
        Path to .swc file.

    Returns
    -------
    n : int
        Number of compartments (nodes).
    adj : list[list[int]]
        0-indexed adjacency list suitable for indpoly.independence_poly().
    metadata : dict
        'labels' maps new index -> original sample ID,
        'types' maps new index -> compartment type (1=soma, 2=axon, ...),
        'coords' maps new index -> (x, y, z, radius),
        'root' is the new index of the root node.
    """
    path = Path(path)
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 7:
                continue
            sample_id = int(parts[0])
            struct_type = int(parts[1])
            x, y, z, r = float(parts[2]), float(parts[3]), float(parts[4]), float(parts[5])
            parent_id = int(parts[6])
            rows.append((sample_id, struct_type, x, y, z, r, parent_id))

    if not rows:
        raise ValueError(f"No valid SWC data in {path}")

    # Map original sample IDs to 0-indexed
    id_to_idx = {}
    for i, (sample_id, *_) in enumerate(rows):
        id_to_idx[sample_id] = i

    n = len(rows)
    adj = [[] for _ in range(n)]
    labels = {}
    types = {}
    coords = {}
    root = -1

    for sample_id, struct_type, x, y, z, r, parent_id in rows:
        idx = id_to_idx[sample_id]
        labels[idx] = sample_id
        types[idx] = struct_type
        coords[idx] = (x, y, z, r)

        if parent_id == -1:
            root = idx
        else:
            if parent_id not in id_to_idx:
                raise ValueError(
                    f"Parent {parent_id} not found for sample {sample_id}"
                )
            pidx = id_to_idx[parent_id]
            adj[idx].append(pidx)
            adj[pidx].append(idx)

    meta = {"labels": labels, "types": types, "coords": coords, "root": root}
    return n, adj, meta


# ---------------------------------------------------------------------------
# Newick reader
# ---------------------------------------------------------------------------

def read_newick(source: str | Path) -> tuple[int, list[list[int]], dict]:
    """Read a Newick string or file and return (n, adj, metadata).

    Handles standard Newick features: labels, branch lengths, internal
    node labels, quoted labels. Ignores branch lengths (graph-theoretic
    tree only).

    Parameters
    ----------
    source : str or Path
        Either a file path (ending in .nwk, .newick, .nex, .tre, .tree,
        .treefile) or a raw Newick string.

    Returns
    -------
    n : int
        Number of nodes (tips + internal).
    adj : list[list[int]]
        0-indexed adjacency list suitable for indpoly.independence_poly().
    metadata : dict
        'labels' maps node index -> label (empty string if unlabelled),
        'tips' is the list of tip (leaf) indices,
        'internal' is the list of internal node indices,
        'root' is the index of the root node.
    """
    text = _get_newick_string(source)
    nodes, edges = _parse_newick(text)

    n = len(nodes)
    adj = [[] for _ in range(n)]
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)

    tips = [i for i in range(n) if len(adj[i]) <= 1]
    internal = [i for i in range(n) if len(adj[i]) > 1]

    meta = {
        "labels": {i: label for i, label in enumerate(nodes)},
        "tips": tips,
        "internal": internal,
        "root": 0,
    }
    return n, adj, meta


def read_newick_multi(source: str | Path) -> list[tuple[int, list[list[int]], dict]]:
    """Read a file containing multiple Newick trees (one per line).

    Returns a list of (n, adj, metadata) tuples.
    """
    path = Path(source)
    results = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "(" in line:
                results.append(read_newick(line))
    return results


def _get_newick_string(source: str | Path) -> str:
    """Extract a Newick string from a file or raw string."""
    source_str = str(source)
    tree_extensions = {".nwk", ".newick", ".nex", ".nexus", ".tre", ".tree", ".treefile"}

    path = Path(source_str)
    if path.suffix.lower() in tree_extensions and path.exists():
        text = path.read_text()
    elif path.exists() and "(" in path.read_text()[:200]:
        text = path.read_text()
    else:
        text = source_str

    # For NEXUS files, extract the tree block
    if "begin trees" in text.lower():
        text = _extract_nexus_tree(text)

    # Strip to just the Newick string
    text = text.strip()
    # Take first tree if multiple
    if ";" in text:
        text = text[: text.index(";") + 1]
    return text


def _extract_nexus_tree(text: str) -> str:
    """Extract the first tree string from a NEXUS file."""
    in_trees = False
    for line in text.splitlines():
        stripped = line.strip().lower()
        if stripped.startswith("begin trees"):
            in_trees = True
            continue
        if in_trees:
            if stripped.startswith("end"):
                break
            # Tree lines look like: tree NAME = [&R] ((...));
            if stripped.startswith("tree "):
                eq_pos = line.find("=")
                if eq_pos >= 0:
                    tree_str = line[eq_pos + 1 :].strip()
                    # Strip [&R] or [&U] annotations
                    tree_str = re.sub(r"\[&[RU]\]\s*", "", tree_str)
                    return tree_str
    raise ValueError("No tree found in NEXUS data")


def _parse_newick(text: str) -> tuple[list[str], list[tuple[int, int]]]:
    """Parse a Newick string into node labels and edge list.

    Returns (labels, edges) where labels[i] is the label for node i
    and edges is a list of (parent_idx, child_idx) pairs.
    """
    text = text.strip().rstrip(";").strip()
    if not text:
        raise ValueError("Empty Newick string")

    labels: list[str] = []
    edges: list[tuple[int, int]] = []

    def new_node(label: str = "") -> int:
        idx = len(labels)
        labels.append(label)
        return idx

    pos = 0

    def parse_subtree(parent: int | None) -> int:
        nonlocal pos

        if pos < len(text) and text[pos] == "(":
            # Internal node
            node = new_node()
            pos += 1  # skip '('

            while True:
                child = parse_subtree(node)
                edges.append((node, child))

                if pos < len(text) and text[pos] == ",":
                    pos += 1  # skip ','
                elif pos < len(text) and text[pos] == ")":
                    pos += 1  # skip ')'
                    break
                else:
                    break

            # Read optional internal node label
            label = _read_label()
            labels[node] = label

            # Skip branch length
            _skip_branch_length()

            return node
        else:
            # Leaf node
            label = _read_label()
            node = new_node(label)
            _skip_branch_length()
            return node

    def _read_label() -> str:
        nonlocal pos
        if pos >= len(text):
            return ""
        # Quoted label
        if text[pos] == "'":
            pos += 1
            start = pos
            while pos < len(text) and text[pos] != "'":
                pos += 1
            label = text[start:pos]
            if pos < len(text):
                pos += 1  # skip closing quote
            return label
        # Unquoted label
        start = pos
        while pos < len(text) and text[pos] not in ",:;()[]":
            pos += 1
        return text[start:pos].strip()

    def _skip_branch_length():
        nonlocal pos
        if pos < len(text) and text[pos] == ":":
            pos += 1
            while pos < len(text) and text[pos] not in ",;()[]":
                pos += 1

    root = parse_subtree(None)
    return labels, edges


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _print_summary(n, adj, meta, source_name):
    """Print a quick summary of a converted tree."""
    degrees = [len(adj[i]) for i in range(n)]
    leaves = sum(1 for d in degrees if d <= 1)
    max_deg = max(degrees) if degrees else 0
    print(f"  {source_name}: n={n}, leaves={leaves}, max_degree={max_deg}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python bio_trees.py <file.swc|file.nwk|file.newick|...>")
        print("Reads a biological tree file and prints summary statistics.")
        sys.exit(1)

    for fpath in sys.argv[1:]:
        path = Path(fpath)
        if not path.exists():
            print(f"  {fpath}: file not found", file=sys.stderr)
            continue

        try:
            if path.suffix.lower() == ".swc":
                n, adj, meta = read_swc(path)
                _print_summary(n, adj, meta, path.name)
            else:
                n, adj, meta = read_newick(path)
                _print_summary(n, adj, meta, path.name)
        except Exception as e:
            print(f"  {path.name}: error: {e}", file=sys.stderr)
