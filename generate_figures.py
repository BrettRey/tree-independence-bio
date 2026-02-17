"""Generate results table and figures for the biological trees paper.

Produces:
  1. LaTeX table of all computational results  (stdout + figures/results_table.tex)
  2. Scatter plot of near-miss ratio vs n       (figures/nm_vs_n.pdf)
  3. Normalised polynomial shape overlay        (figures/poly_shapes.pdf)
  4. P(v) vs degree scatter plot                (figures/pv_vs_degree.pdf)
  5. nm regression with named residuals         (figures/nm_regression.pdf)
  6. Null-model nm comparison box plot          (figures/null_model_nm.pdf)

Usage:
    python generate_figures.py
"""

import json
import sys
from pathlib import Path

# Paths
PROJECT = Path(__file__).parent
HOUSE_STYLE = PROJECT.parent.parent / ".house-style"
FIGURES = PROJECT / "figures"

sys.path.insert(0, str(PROJECT))
sys.path.insert(0, str(HOUSE_STYLE))

from bio_trees import read_swc, read_newick
from indpoly import independence_poly, near_miss_ratio
from occupation import occupation_probabilities

import matplotlib
matplotlib.use("Agg")  # non-interactive backend
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import ScalarFormatter
import numpy as np

from plot_style import setup, COLORS, add_grid, save_figure

# Apply house style
setup()


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
    n_neuro = [row[2] for row in NEURONS]
    nm_neuro = [row[5] for row in NEURONS]
    n_phylo = [row[1] for row in PHYLOGENIES]
    nm_phylo = [row[4] for row in PHYLOGENIES]

    fig, ax = plt.subplots(figsize=(5.5, 3.8))

    ax.scatter(n_neuro, nm_neuro, marker="o", s=36, c=COLORS['primary'],
               edgecolors=COLORS['primary'], linewidths=0.6, zorder=3,
               label="Neuronal arbors")

    ax.scatter(n_phylo, nm_phylo, marker="^", s=42, facecolors="none",
               edgecolors=COLORS['secondary'], linewidths=1.0, zorder=3,
               label="Phylogenetic trees")

    # Reference curve: nm = 1 - C/n
    C = 6
    n_ref = np.linspace(30, 40000, 500)
    nm_ref = 1.0 - C / n_ref
    ax.plot(n_ref, nm_ref, "--", color=COLORS['dark'], linewidth=0.9,
            zorder=2, alpha=0.5, label=f"$1 - {C}/n$")

    ax.axhline(y=1.0, color=COLORS['dark'], linestyle=":", linewidth=0.7,
               zorder=1)

    ax.set_xscale("log")
    ax.set_xlim(30, 50000)
    ax.set_ylim(0.70, 1.01)
    ax.set_xlabel("Number of nodes $n$")
    ax.set_ylabel("Near-miss ratio nm")

    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.get_major_formatter().set_scientific(False)
    ax.xaxis.get_major_formatter().set_useOffset(False)
    ax.set_xticks([50, 100, 500, 1000, 5000, 10000, 30000])
    ax.set_xticklabels(["50", "100", "500", "1,000", "5,000",
                        "10,000", "30,000"])

    add_grid(ax)
    ax.legend(loc="lower right")

    save_figure(fig, FIGURES / "nm_vs_n", formats=("pdf",))
    plt.close(fig)


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

    representatives = [
        ("phylo", "Homininae ($n=46$)", "homininae.nwk"),
        ("neuron", "Drosophila ddaC A5 ($n=113$)",
         "20_062909e_02_3rd_ddaC_A5.CNG.swc"),
        ("phylo", "Accipitridae ($n=884$)", "accipitridae.nwk"),
        ("neuron", "Monkey L3 pyramidal ($n=1038$)", "cnic_041.CNG.swc"),
    ]

    styles = [
        {"color": COLORS['secondary'], "linestyle": "-"},
        {"color": COLORS['tertiary'], "linestyle": "--"},
        {"color": COLORS['quaternary'], "linestyle": "-."},
        {"color": COLORS['quinary'], "linestyle": ":"},
    ]

    fig, ax = plt.subplots(figsize=(5.5, 3.8))

    for (cat, label, fname), style in zip(representatives, styles):
        print(f"  Computing polynomial for {label}...")
        n, adj = _load_tree(cat, fname)
        poly = independence_poly(n, adj)
        alpha = len(poly) - 1
        max_coeff = max(poly)

        x = np.array([k / alpha for k in range(alpha + 1)])
        y = np.array([c / max_coeff for c in poly])

        ax.plot(x, y, label=label, **style)

    ax.set_xlabel(r"$k / \alpha$ (normalised index)")
    ax.set_ylabel(r"$i_k / \max(i_k)$ (normalised coefficient)")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.05)

    add_grid(ax)
    ax.legend(loc="upper right", fontsize=8)

    save_figure(fig, FIGURES / "poly_shapes", formats=("pdf",))
    plt.close(fig)


# ---------------------------------------------------------------------------
# 4. P(v) vs degree scatter plot
# ---------------------------------------------------------------------------

def plot_pv_vs_degree():
    """Scatter plot of occupation probability P(v) vs vertex degree."""
    representatives = [
        ("phylo", "Homininae ($n=46$)", "homininae.nwk"),
        ("neuron", "Monkey L3 pyramidal ($n=1{,}038$)", "cnic_041.CNG.swc"),
        ("phylo", "Primates ($n=1{,}333$)", "primates.nwk"),
    ]

    markers = [
        {"marker": "^", "s": 18, "facecolors": "none",
         "edgecolors": COLORS['secondary'], "linewidths": 0.7},
        {"marker": "o", "s": 10, "facecolors": COLORS['primary'],
         "edgecolors": COLORS['primary'], "linewidths": 0.4},
        {"marker": "s", "s": 12, "facecolors": "none",
         "edgecolors": COLORS['quaternary'], "linewidths": 0.7},
    ]

    fig, ax = plt.subplots(figsize=(5.5, 3.8))

    for (cat, label, fname), mstyle in zip(representatives, markers):
        print(f"  Computing P(v) for {label}...")
        n, adj, meta = _load_tree_with_meta(cat, fname)
        root = meta.get("root", 0)
        probs = occupation_probabilities(n, adj, root=root)
        degrees = [len(adj[i]) for i in range(n)]

        deg_jitter = np.array(degrees) + np.random.default_rng(42).uniform(
            -0.15, 0.15, size=n)

        ax.scatter(deg_jitter, probs, label=label, alpha=0.6, zorder=3,
                   **mstyle)

    # Reference line at P(v) = 1/3
    ax.axhline(y=1 / 3, color=COLORS['dark'], linestyle=":", linewidth=0.7,
               zorder=1, alpha=0.5)
    ax.text(ax.get_xlim()[1] * 0.85, 1 / 3 + 0.01, "$1/3$",
            fontsize=8, color=COLORS['dark'], alpha=0.6, ha="right")

    # Theoretical lower envelope: P(hub) = 1/(2^k + 1) for star K_{1,k}
    k_vals = np.arange(1, 15)
    p_star = 1.0 / (2.0 ** k_vals + 1.0)
    ax.plot(k_vals, p_star, "-", color=COLORS['dark'], linewidth=1.0,
            zorder=2, alpha=0.4, marker=".", markersize=3,
            label=r"$K_{1,k}$ hub: $1/(2^k+1)$")

    ax.set_xlabel("Vertex degree")
    ax.set_ylabel("Occupation probability $P(v)$")
    ax.set_xlim(0.5, None)
    ax.set_ylim(-0.01, 0.55)

    add_grid(ax)
    ax.legend(loc="upper right", fontsize=8)

    save_figure(fig, FIGURES / "pv_vs_degree", formats=("pdf",))
    plt.close(fig)


# ---------------------------------------------------------------------------
# 5. nm regression with named residuals
# ---------------------------------------------------------------------------

def plot_nm_regression():
    """Two-panel figure: nm vs n with fitted 1 - C/n curve, and named
    residuals."""
    from scipy.optimize import curve_fit

    names, ns, nms, cats = [], [], [], []
    for row in NEURONS:
        names.append(row[0])
        ns.append(row[2])
        nms.append(row[5])
        cats.append("neuron")
    for row in PHYLOGENIES:
        names.append(row[0])
        ns.append(row[1])
        nms.append(row[4])
        cats.append("phylo")

    ns = np.array(ns, dtype=float)
    nms = np.array(nms, dtype=float)

    def model(n, C):
        return 1.0 - C / n

    popt, pcov = curve_fit(model, ns, nms, p0=[6.0])
    C_fit = popt[0]
    C_se = np.sqrt(pcov[0, 0])
    nm_pred = model(ns, C_fit)
    residuals = nms - nm_pred

    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((nms - np.mean(nms)) ** 2)
    r_sq = 1.0 - ss_res / ss_tot

    print(f"  Fit: nm = 1 - {C_fit:.2f}/n  (SE = {C_se:.2f}, "
          f"R\u00b2 = {r_sq:.4f})")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6.0, 5.5),
                                    height_ratios=[2, 1],
                                    gridspec_kw={"hspace": 0.35})

    # --- Top panel: nm vs n with fit ---
    n_curve = np.linspace(30, 40000, 500)
    nm_curve = model(n_curve, C_fit)

    ax1.plot(n_curve, nm_curve, "-", color=COLORS['dark'], linewidth=1.0,
             zorder=2, alpha=0.5,
             label=f"$1 - {C_fit:.1f}/n$  ($R^2 = {r_sq:.3f}$)")

    for i, cat in enumerate(cats):
        colour = COLORS['primary'] if cat == "neuron" else COLORS['secondary']
        marker = "o" if cat == "neuron" else "^"
        ax1.scatter(ns[i], nms[i], marker=marker, s=30, c=colour,
                    edgecolors=colour, linewidths=0.5, zorder=3)

    ax1.axhline(y=1.0, color=COLORS['dark'], linestyle=":", linewidth=0.7,
                zorder=1)
    ax1.set_xscale("log")
    ax1.set_xlim(30, 50000)
    ax1.set_ylim(0.70, 1.01)
    ax1.set_ylabel("Near-miss ratio nm")
    ax1.xaxis.set_major_formatter(ScalarFormatter())
    ax1.xaxis.get_major_formatter().set_scientific(False)
    ax1.set_xticks([50, 100, 500, 1000, 5000, 10000, 30000])
    ax1.set_xticklabels(["50", "100", "500", "1k", "5k", "10k", "30k"])
    add_grid(ax1)

    legend_elements = [
        Line2D([0], [0], marker="o", color="w",
               markerfacecolor=COLORS['primary'], markersize=5,
               label="Neuronal arbors"),
        Line2D([0], [0], marker="^", color="w",
               markerfacecolor=COLORS['secondary'], markersize=5,
               label="Phylogenies"),
        Line2D([0], [0], color=COLORS['dark'], linewidth=1.0, alpha=0.5,
               label=f"$1 - {C_fit:.1f}/n$  ($R^2 = {r_sq:.3f}$)"),
    ]
    ax1.legend(handles=legend_elements, loc="lower right", fontsize=8)

    # --- Bottom panel: named residuals ---
    ax2.axhline(y=0, color=COLORS['dark'], linestyle="-", linewidth=0.5,
                zorder=1, alpha=0.4)

    for i, cat in enumerate(cats):
        colour = COLORS['primary'] if cat == "neuron" else COLORS['secondary']
        marker = "o" if cat == "neuron" else "^"
        ax2.scatter(ns[i], residuals[i], marker=marker, s=24, c=colour,
                    edgecolors=colour, linewidths=0.5, zorder=3)

    short_names = []
    for nm in names:
        s = nm.replace("Drosophila ddaC ", "Dros. ")
        s = s.replace("Monkey L3 pyramidal ", "Monkey ")
        s = s.replace("Human aspiny interneuron", "Hum. aspiny")
        s = s.replace("Mouse Purkinje ", "Purk. ")
        s = s.replace("Rat interneuron (IDC)", "Rat IDC")
        s = s.replace("Rat interneuron (A3)", "Rat A3")
        s = s.replace("Rat interneuron", "Rat int.")
        s = s.replace("Rat pyramidal (n419)", "Rat pyr.")
        s = s.replace("Rat basket interneuron", "Rat basket")
        s = s.replace("Mouse principal cell", "Mouse princ.")
        s = s.replace("Human pyramidal ", "Hum. pyr. ")
        short_names.append(s)

    for i, sn in enumerate(short_names):
        y_off = 4 if residuals[i] >= 0 else -8
        ax2.annotate(sn, (ns[i], residuals[i]),
                     fontsize=5.0, color=COLORS['dark'], alpha=0.6,
                     ha="center", textcoords="offset points",
                     xytext=(0, y_off), rotation=45)

    ax2.set_xscale("log")
    ax2.set_xlim(30, 50000)
    ax2.set_xlabel("Number of nodes $n$")
    ax2.set_ylabel("Residual")
    ax2.xaxis.set_major_formatter(ScalarFormatter())
    ax2.xaxis.get_major_formatter().set_scientific(False)
    ax2.set_xticks([50, 100, 500, 1000, 5000, 10000, 30000])
    ax2.set_xticklabels(["50", "100", "500", "1k", "5k", "10k", "30k"])
    add_grid(ax2)

    save_figure(fig, FIGURES / "nm_regression", formats=("pdf",))
    plt.close(fig)


# ---------------------------------------------------------------------------
# 6. Null-model nm comparison
# ---------------------------------------------------------------------------

def plot_null_model_nm():
    """Box plot comparing biological nm to null-model distributions."""
    results_path = PROJECT / "null_model_results_small.json"
    if not results_path.exists():
        results_path = PROJECT / "null_model_results_quick.json"
    if not results_path.exists():
        results_path = PROJECT / "null_model_results.json"

    if not results_path.exists():
        print("  No null-model results found. Run null_model.py first.")
        return

    with open(results_path) as f:
        raw = json.load(f)

    bio_nm = {}
    bio_cat = {}
    for row in NEURONS:
        bio_nm[row[2]] = row[5]
        bio_cat[row[2]] = "neuron"
    for row in PHYLOGENIES:
        bio_nm[row[1]] = row[4]
        bio_cat[row[1]] = "phylo"

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
        bp = ax.boxplot([null_nms], positions=[i], widths=0.5,
                        patch_artist=True, showfliers=False,
                        boxprops=dict(facecolor=COLORS['light'],
                                      edgecolor=COLORS['dark'],
                                      linewidth=0.6),
                        medianprops=dict(color=COLORS['dark'],
                                         linewidth=0.8),
                        whiskerprops=dict(color=COLORS['dark'],
                                          linewidth=0.6),
                        capprops=dict(color=COLORS['dark'],
                                      linewidth=0.6))

        colour = COLORS['primary'] if cat == "neuron" else COLORS['secondary']
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

    ax.axhline(y=1.0, color=COLORS['dark'], linestyle=":", linewidth=0.7,
               zorder=1)

    add_grid(ax)

    legend_elements = [
        Line2D([0], [0], marker="o", color="w",
               markerfacecolor=COLORS['primary'], markersize=6,
               label="Neuronal (observed)"),
        Line2D([0], [0], marker="^", color="w",
               markerfacecolor=COLORS['secondary'], markersize=6,
               label="Phylogeny (observed)"),
        plt.Rectangle((0, 0), 1, 1, fc=COLORS['light'],
                       ec=COLORS['dark'], linewidth=0.6,
                       label="Random tree null"),
    ]
    ax.legend(handles=legend_elements, loc="lower right", fontsize=8)

    save_figure(fig, FIGURES / "null_model_nm", formats=("pdf",))
    plt.close(fig)


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
    print("Generating Figure 4: nm regression with named residuals")
    print("=" * 60)
    plot_nm_regression()

    print()
    print("=" * 60)
    print("Generating Figure 5: null-model nm comparison")
    print("=" * 60)
    plot_null_model_nm()

    print()
    print("All outputs in:", FIGURES)


if __name__ == "__main__":
    main()
