"""Generate results table and figures for the biological trees paper.

Produces:
  1. LaTeX table of all computational results  (stdout + figures/results_table.tex)
  2. Scatter plot of near-miss ratio vs n       (figures/nm_vs_n.pdf)
  3. Normalised polynomial shape overlay        (figures/poly_shapes.pdf)
  4. P(v) vs degree scatter plot                (figures/pv_vs_degree.pdf)
  5. Null-model nm comparison box plot          (figures/null_model_nm.pdf)

Usage:
    python generate_figures.py
"""

import json
import sys
from pathlib import Path

# Paths
PROJECT = Path(__file__).parent
ERDOS = PROJECT.parent / "Erdos_Problem_993"
FIGURES = PROJECT / "figures"

sys.path.insert(0, str(PROJECT))
sys.path.insert(0, str(ERDOS))

from bio_trees import read_swc, read_newick
from indpoly import independence_poly, near_miss_ratio
from occupation import occupation_probabilities

import matplotlib
matplotlib.use("Agg")  # non-interactive backend
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import numpy as np


# ---------------------------------------------------------------------------
# Matplotlib style: clean, publication-quality, serif font
# ---------------------------------------------------------------------------

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["STIX Two Text", "Times", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 12,
    "legend.fontsize": 9,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
})


# ---------------------------------------------------------------------------
# Data: hardcoded from the completed computational analysis
# ---------------------------------------------------------------------------

# Neuronal arbors (NeuroMorpho.org)
# Columns: display_name, species, n, alpha, mode, nm, swc_filename
NEURONS = [
    ("Drosophila ddaC (A5)", "Drosophila", 113, 58, 32, 0.8478,
     "20_062909e_02_3rd_ddaC_A5.CNG.swc"),
    ("Drosophila ddaC (A6)", "Drosophila", 166, 83, 46, 0.9317,
     "06_060609y_9-A6-ddaC.CNG.swc"),
    ("Monkey L3 pyramidal (041)", "Macaque", 1038, 525, 290, 0.9813,
     "cnic_041.CNG.swc"),
    ("Monkey L3 pyramidal (001)", "Macaque", 1274, 646, 356, 0.9873,
     "cnic_001.CNG.swc"),
    ("Human aspiny interneuron", "Human", 1498, 752, 416, 0.9865,
     "H17-03-010-11-13-01_656411100_m.CNG.swc"),
    ("Mouse Purkinje (P35)", "Mouse", 1716, 891, 491, 0.9896,
     "Purkinje-slice-ageP35-1.CNG.swc"),
    ("Rat interneuron", "Rat", 1941, 978, 540, 0.9927,
     "SM080903A2.CNG.swc"),
    ("Mouse principal cell", "Mouse", 2212, 1108, 612, 0.9944,
     "TypeA-10.CNG.swc"),
    ("Rat interneuron (IDC)", "Rat", 2774, 1392, 770, 0.9930,
     "SM080902A1-5_IDC.CNG.swc"),
    ("Rat pyramidal (n419)", "Rat", 4229, 2122, 1174, 0.9958,
     "n419.CNG.swc"),
    ("Mouse Purkinje (P43)", "Mouse", 4156, 2126, 1172, 0.9954,
     "Purkinje-slice-ageP43-6.CNG.swc"),
    ("Rat basket interneuron", "Rat", 6378, 3199, 1767, 0.9970,
     "5-27-11cell1-nonFS-BC.CNG.swc"),
    ("Rat interneuron (A3)", "Rat", 7649, 3843, 2123, 0.9981,
     "SM080902A3-1.CNG.swc"),
    ("Human pyramidal (06)", "Human", 10451, 5236, 2895, 0.9980,
     "H17-03-010-11-13-06_651089035_m.CNG.swc"),
    ("Human pyramidal (05)", "Human", 17992, 9011, 4981, 0.9990,
     "H17-03-011-11-04-05_650978964_m.CNG.swc"),
]

# Phylogenetic trees (Open Tree of Life)
# Columns: display_name, n, alpha, mode, nm, nwk_filename
PHYLOGENIES = [
    ("Homininae", 46, 28, 15, 0.7700, "homininae.nwk"),
    ("Delphinidae", 107, 71, 36, 0.9184, "delphinidae.nwk"),
    ("Felidae", 144, 100, 52, 0.9479, "felidae.nwk"),
    ("Salamandridae", 157, 105, 54, 0.9423, "salamandridae.nwk"),
    ("Mustelidae", 188, 122, 64, 0.9488, "mustelidae.nwk"),
    ("Canidae", 367, 241, 125, 0.9812, "canidae.nwk"),
    ("Equidae", 372, 223, 117, 0.9612, "equidae.nwk"),
    ("Accipitridae", 884, 601, 311, 0.9857, "accipitridae.nwk"),
    ("Papilionidae", 1170, 929, 472, 0.9917, "papilionidae.nwk"),
    ("Primates", 1333, 872, 452, 0.9941, "primates.nwk"),
    ("Cetacea", 1387, 971, 498, 0.9911, "cetacea.nwk"),
    ("Aves", 30310, 20357, 10532, 0.9997, "aves.nwk"),
]


# ---------------------------------------------------------------------------
# 1. LaTeX table
# ---------------------------------------------------------------------------

def generate_latex_table():
    """Generate a LaTeX results table and write to file."""
    lines = []
    lines.append(r"\begin{table}[htbp]")
    lines.append(r"\centering")
    lines.append(r"\caption{Independence polynomial properties of biological"
                 r" trees. $n$~= number of nodes, $\alpha$~= independence"
                 r" number, mode~= index of largest coefficient,"
                 r" $\mathrm{nm}$~= near-miss ratio.}")
    lines.append(r"\label{tab:results}")
    lines.append(r"\small")
    lines.append(r"\begin{tabular}{@{}llrrrr@{}}")
    lines.append(r"\toprule")
    lines.append(r"Tree & Species & $n$ & $\alpha$ & mode & $\mathrm{nm}$ \\")
    lines.append(r"\midrule")
    lines.append(r"\multicolumn{6}{@{}l}{\textit{Neuronal arbors"
                 r" (NeuroMorpho.org)}} \\")
    lines.append(r"\addlinespace[2pt]")

    for name, species, n, alpha, mode, nm, _ in NEURONS:
        lines.append(
            f"  {name} & {species} & {n:,} & {alpha:,} & {mode:,}"
            f" & {nm:.4f} \\\\"
        )

    lines.append(r"\addlinespace[6pt]")
    lines.append(r"\multicolumn{6}{@{}l}{\textit{Phylogenetic trees"
                 r" (Open Tree of Life)}} \\")
    lines.append(r"\addlinespace[2pt]")

    for name, n, alpha, mode, nm, _ in PHYLOGENIES:
        lines.append(
            f"  {name} & --- & {n:,} & {alpha:,} & {mode:,}"
            f" & {nm:.4f} \\\\"
        )

    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")

    table_tex = "\n".join(lines)

    out_path = FIGURES / "results_table.tex"
    out_path.write_text(table_tex)
    print(f"LaTeX table written to {out_path}")
    print()
    print(table_tex)
    return table_tex


# ---------------------------------------------------------------------------
# 2. Near-miss ratio vs n (scatter plot)
# ---------------------------------------------------------------------------

def plot_nm_vs_n():
    """Scatter plot of near-miss ratio vs number of nodes."""
    # Extract data
    n_neuro = [row[2] for row in NEURONS]
    nm_neuro = [row[5] for row in NEURONS]
    n_phylo = [row[1] for row in PHYLOGENIES]
    nm_phylo = [row[4] for row in PHYLOGENIES]

    fig, ax = plt.subplots(figsize=(5.5, 3.8))

    # Neuronal arbors: filled circles
    ax.scatter(n_neuro, nm_neuro, marker="o", s=36, c="#2166ac",
               edgecolors="#2166ac", linewidths=0.6, zorder=3,
               label="Neuronal arbors")

    # Phylogenies: open triangles
    ax.scatter(n_phylo, nm_phylo, marker="^", s=42, facecolors="none",
               edgecolors="#b2182b", linewidths=1.0, zorder=3,
               label="Phylogenetic trees")

    # Reference curve: nm = 1 - C/n
    C = 6
    n_ref = np.linspace(30, 40000, 500)
    nm_ref = 1.0 - C / n_ref
    ax.plot(n_ref, nm_ref, "--", color="0.55", linewidth=0.9, zorder=2,
            label=f"$1 - {C}/n$")

    # Violation threshold
    ax.axhline(y=1.0, color="0.3", linestyle=":", linewidth=0.7, zorder=1)

    ax.set_xscale("log")
    ax.set_xlim(30, 50000)
    ax.set_ylim(0.70, 1.01)

    ax.set_xlabel("Number of nodes $n$")
    ax.set_ylabel("Near-miss ratio nm")

    # Clean up x-axis tick formatting (no scientific notation)
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.get_major_formatter().set_scientific(False)
    ax.xaxis.get_major_formatter().set_useOffset(False)
    ax.set_xticks([50, 100, 500, 1000, 5000, 10000, 30000])
    ax.set_xticklabels(["50", "100", "500", "1,000", "5,000", "10,000", "30,000"])

    # Minimal gridlines
    ax.yaxis.grid(True, linewidth=0.3, color="0.85", zorder=0)
    ax.xaxis.grid(False)

    ax.legend(loc="lower right", frameon=True, framealpha=0.95,
              edgecolor="0.8", fancybox=False)

    # Remove top and right spines for cleaner look
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(top=False, right=False)

    out_path = FIGURES / "nm_vs_n.pdf"
    fig.savefig(out_path)
    plt.close(fig)
    print(f"Scatter plot saved to {out_path}")


# ---------------------------------------------------------------------------
# 3. Normalised polynomial shapes (overlay)
# ---------------------------------------------------------------------------

def _load_tree(category, filename):
    """Load a biological tree from the data directory."""
    if category == "neuron":
        path = PROJECT / "data" / "neuromorpho" / filename
        n, adj, meta = read_swc(path)
    else:
        path = PROJECT / "data" / "phylogenies" / filename
        n, adj, meta = read_newick(path)
    return n, adj


def _load_tree_with_meta(category, filename):
    """Load a biological tree, returning metadata too."""
    if category == "neuron":
        path = PROJECT / "data" / "neuromorpho" / filename
        n, adj, meta = read_swc(path)
    else:
        path = PROJECT / "data" / "phylogenies" / filename
        n, adj, meta = read_newick(path)
    return n, adj, meta


def plot_poly_shapes():
    """Overlay normalised polynomial shapes for 4 representative trees."""

    # Representative trees (category, display_name, filename)
    representatives = [
        ("phylo", "Homininae ($n=46$)", "homininae.nwk"),
        ("neuron", "Drosophila ddaC A5 ($n=113$)",
         "20_062909e_02_3rd_ddaC_A5.CNG.swc"),
        ("phylo", "Accipitridae ($n=884$)", "accipitridae.nwk"),
        ("neuron", "Monkey L3 pyramidal ($n=1038$)", "cnic_041.CNG.swc"),
    ]

    # Line styles for visual distinction
    styles = [
        {"color": "#d95f02", "linestyle": "-", "linewidth": 1.4},
        {"color": "#1b9e77", "linestyle": "--", "linewidth": 1.4},
        {"color": "#7570b3", "linestyle": "-.", "linewidth": 1.4},
        {"color": "#e7298a", "linestyle": ":", "linewidth": 1.8},
    ]

    fig, ax = plt.subplots(figsize=(5.5, 3.8))

    for (cat, label, fname), style in zip(representatives, styles):
        print(f"  Computing polynomial for {label}...")
        n, adj = _load_tree(cat, fname)
        poly = independence_poly(n, adj)
        alpha = len(poly) - 1
        max_coeff = max(poly)

        # Normalise
        x = np.array([k / alpha for k in range(alpha + 1)])
        y = np.array([c / max_coeff for c in poly])

        ax.plot(x, y, label=label, **style)

    ax.set_xlabel(r"$k / \alpha$ (normalised index)")
    ax.set_ylabel(r"$i_k / \max(i_k)$ (normalised coefficient)")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.05)

    ax.legend(loc="upper right", frameon=True, framealpha=0.95,
              edgecolor="0.8", fancybox=False, fontsize=8)

    # Minimal gridlines
    ax.yaxis.grid(True, linewidth=0.3, color="0.85", zorder=0)
    ax.xaxis.grid(False)

    # Remove top and right spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(top=False, right=False)

    out_path = FIGURES / "poly_shapes.pdf"
    fig.savefig(out_path)
    plt.close(fig)
    print(f"Polynomial shapes saved to {out_path}")


# ---------------------------------------------------------------------------
# 4. P(v) vs degree scatter plot
# ---------------------------------------------------------------------------

def plot_pv_vs_degree():
    """Scatter plot of occupation probability P(v) vs vertex degree.

    Shows 3 representative trees (one small phylogeny, one neuronal arbor,
    one larger phylogeny) to illustrate Hub Exclusion quantitatively:
    high-degree vertices cluster at low P(v).
    """
    representatives = [
        ("phylo", "Homininae ($n=46$)", "homininae.nwk"),
        ("neuron", "Monkey L3 pyramidal ($n=1{,}038$)", "cnic_041.CNG.swc"),
        ("phylo", "Primates ($n=1{,}333$)", "primates.nwk"),
    ]

    markers = [
        {"marker": "^", "s": 18, "facecolors": "none",
         "edgecolors": "#d95f02", "linewidths": 0.7},
        {"marker": "o", "s": 10, "facecolors": "#2166ac",
         "edgecolors": "#2166ac", "linewidths": 0.4},
        {"marker": "s", "s": 12, "facecolors": "none",
         "edgecolors": "#7570b3", "linewidths": 0.7},
    ]

    fig, ax = plt.subplots(figsize=(5.5, 3.8))

    for (cat, label, fname), mstyle in zip(representatives, markers):
        print(f"  Computing P(v) for {label}...")
        n, adj, meta = _load_tree_with_meta(cat, fname)
        root = meta.get("root", 0)
        probs = occupation_probabilities(n, adj, root=root)
        degrees = [len(adj[i]) for i in range(n)]

        # Jitter degree slightly for visibility
        deg_jitter = np.array(degrees) + np.random.default_rng(42).uniform(
            -0.15, 0.15, size=n)

        ax.scatter(deg_jitter, probs, label=label, alpha=0.6, zorder=3,
                   **mstyle)

    # Reference line at P(v) = 1/3 (edge bound threshold)
    ax.axhline(y=1 / 3, color="0.5", linestyle=":", linewidth=0.7, zorder=1)
    ax.text(ax.get_xlim()[1] * 0.85, 1 / 3 + 0.01, "$1/3$",
            fontsize=8, color="0.5", ha="right")

    ax.set_xlabel("Vertex degree")
    ax.set_ylabel("Occupation probability $P(v)$")
    ax.set_xlim(0.5, None)
    ax.set_ylim(-0.01, 0.55)

    ax.legend(loc="upper right", frameon=True, framealpha=0.95,
              edgecolor="0.8", fancybox=False, fontsize=8)

    ax.yaxis.grid(True, linewidth=0.3, color="0.85", zorder=0)
    ax.xaxis.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(top=False, right=False)

    out_path = FIGURES / "pv_vs_degree.pdf"
    fig.savefig(out_path)
    plt.close(fig)
    print(f"P(v) vs degree saved to {out_path}")


# ---------------------------------------------------------------------------
# 5. Null-model nm comparison
# ---------------------------------------------------------------------------

def plot_null_model_nm():
    """Box plot comparing biological nm to null-model distributions.

    Reads precomputed null-model results from JSON. Falls back to a
    quick computation if no results file is found.
    """
    # Try to load precomputed results
    results_path = PROJECT / "null_model_results.json"
    if not results_path.exists():
        results_path = PROJECT / "null_model_results_quick.json"
    if not results_path.exists():
        results_path = PROJECT / "null_model_results_small.json"

    if not results_path.exists():
        print("  No null-model results found. Run null_model.py first.")
        return

    with open(results_path) as f:
        raw = json.load(f)

    # Build biological nm lookup
    bio_nm = {}
    bio_cat = {}
    for row in NEURONS:
        bio_nm[row[2]] = row[5]
        bio_cat[row[2]] = "neuron"
    for row in PHYLOGENIES:
        bio_nm[row[1]] = row[4]
        bio_cat[row[1]] = "phylo"

    # Collect Prüfer results (primary null model)
    entries = []
    for key, data in raw.items():
        if data["model"] != "prufer":
            continue
        n = data["n"]
        if n not in bio_nm:
            continue
        entries.append((n, data["nms"], bio_nm[n], bio_cat[n]))

    entries.sort(key=lambda x: x[0])

    if not entries:
        print("  No matching null-model entries found.")
        return

    fig, ax = plt.subplots(figsize=(6.5, 4.0))

    positions = list(range(len(entries)))
    labels = []

    for i, (n, null_nms, bio, cat) in enumerate(entries):
        # Box plot for null distribution
        bp = ax.boxplot([null_nms], positions=[i], widths=0.5,
                        patch_artist=True, showfliers=False,
                        boxprops=dict(facecolor="#d4e6f1", edgecolor="0.4",
                                      linewidth=0.6),
                        medianprops=dict(color="0.3", linewidth=0.8),
                        whiskerprops=dict(color="0.4", linewidth=0.6),
                        capprops=dict(color="0.4", linewidth=0.6))

        # Overlay biological nm as a coloured point
        colour = "#2166ac" if cat == "neuron" else "#b2182b"
        marker = "o" if cat == "neuron" else "^"
        ax.scatter([i], [bio], marker=marker, s=40, c=colour,
                   edgecolors=colour, linewidths=0.6, zorder=5)

        if n >= 1000:
            labels.append(f"{n // 1000}k")
        else:
            labels.append(str(n))

    ax.set_xticks(positions)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
    ax.set_xlabel("Tree size $n$")
    ax.set_ylabel("Near-miss ratio nm")

    # Violation threshold
    ax.axhline(y=1.0, color="0.3", linestyle=":", linewidth=0.7, zorder=1)

    ax.yaxis.grid(True, linewidth=0.3, color="0.85", zorder=0)
    ax.xaxis.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(top=False, right=False)

    # Legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#2166ac",
               markersize=6, label="Neuronal (observed)"),
        Line2D([0], [0], marker="^", color="w", markerfacecolor="#b2182b",
               markersize=6, label="Phylogeny (observed)"),
        plt.Rectangle((0, 0), 1, 1, fc="#d4e6f1", ec="0.4", linewidth=0.6,
                       label="Random tree null"),
    ]
    ax.legend(handles=legend_elements, loc="lower right", frameon=True,
              framealpha=0.95, edgecolor="0.8", fancybox=False, fontsize=8)

    out_path = FIGURES / "null_model_nm.pdf"
    fig.savefig(out_path)
    plt.close(fig)
    print(f"Null-model comparison saved to {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    FIGURES.mkdir(exist_ok=True)

    print("=" * 60)
    print("Generating results table")
    print("=" * 60)
    generate_latex_table()

    print()
    print("=" * 60)
    print("Generating Figure 1: nm vs n scatter plot")
    print("=" * 60)
    plot_nm_vs_n()

    print()
    print("=" * 60)
    print("Generating Figure 2: normalised polynomial shapes")
    print("=" * 60)
    plot_poly_shapes()

    print()
    print("=" * 60)
    print("Generating Figure 3: P(v) vs degree scatter plot")
    print("=" * 60)
    plot_pv_vs_degree()

    print()
    print("=" * 60)
    print("Generating Figure 4: null-model nm comparison")
    print("=" * 60)
    plot_null_model_nm()

    print()
    print("All outputs in:", FIGURES)


if __name__ == "__main__":
    main()
