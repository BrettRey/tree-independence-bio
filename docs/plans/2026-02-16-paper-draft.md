# Paper Draft: Tree Independence Polynomials and Biological Network Motifs

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Draft a complete paper translating the Erdos 993 independence polynomial results into biological network language, grounded in computational results on real neuronal arbors and phylogenetic trees.

**Architecture:** Five-section paper (Intro, Background, Results, Discussion, Conclusion) targeting a computational biology audience. The math is established in the companion paper (cited via GitHub repo); this paper's contribution is the biological interpretation plus computational demonstration on 25+ real biological trees. Each section drafted in LaTeX following house style, with references added to `references-local.bib`.

**Tech Stack:** XeLaTeX, biber, house-style preamble (EB Garamond, APA citations), `amsmath` for formulas, `booktabs` for tables, `graphicx` for figures.

---

## Source Material

All mathematical results come from the companion paper:
- **Repo:** `papers/Erdos_Problem_993/` (cite as GitHub repo until arXiv preprint available)
- **Main tex:** `papers/Erdos_Problem_993/paper/main_v2.tex`

Computational results from this project:
- **Analysis script:** `analyze_bio.py`
- **Converters:** `bio_trees.py`
- **Data:** `data/neuromorpho/` (15 SWC files), `data/phylogenies/` (12 Newick files)

Key reviewer input shaping the paper:
- Simulated Uri Alon review (network motifs, SIMs, robustness, dilution)
- Simulated Noga Alon review (original motivation for biological connections)

## Computational Results Summary

25 biological trees analysed, all unimodal and log-concave:
- 15 neuronal arbors (113--17,992 nodes), 5 species, 5 cell types
- 11 phylogenies (46--1,387 nodes; Aves 30,310 still computing)
- Near-miss ratio ranges from 0.770 (Homininae, n=46) to 0.999 (human pyramidal, n=17,992)
- Hub Exclusion Lemma applies: phylogenies have up to 222-degree polytomies with many hubs

---

### Task 1: Set Up Document Structure and Math Macros

**Files:**
- Modify: `main.tex`

**Step 1: Add math macros and restructure sections**

Add to `main.tex` after `\input{.house-style/preamble.tex}`:

```latex
% Project-specific macros
\newcommand{\indpoly}{I(T;\,x)}
\newcommand{\nm}{\mathrm{nm}}
\newcommand{\occ}{P}
\DeclareMathOperator{\mode}{mode}
```

Restructure sections to:

```
\section{Introduction}
\section{Mathematical background}
\subsection{Independence polynomials and unimodality}
\subsection{The hard-core model}
\subsection{Structural reduction: the 1-Private framework}
\section{Biological trees as independence structures}
\subsection{Neuronal dendritic arbors}
\subsection{Phylogenetic trees}
\subsection{Transcription regulatory motifs}
\section{Computational results}
\subsection{Data and methods}
\subsection{Unimodality and log-concavity}
\subsection{Near-miss ratio and robustness}
\subsection{Hub exclusion in biological networks}
\section{Discussion}
\subsection{The independence polynomial as a network descriptor}
\subsection{Universality: physics, chemistry, biology}
\subsection{Limitations and open questions}
\section{Conclusion}
```

**Step 2: Build to verify no LaTeX errors**

Run: `make quick`
Expected: PDF compiles with section headings, no content yet.

**Step 3: Commit**

```bash
git add main.tex
git commit -m "restructure paper sections and add math macros"
```

---

### Task 2: Write `references-local.bib`

**Files:**
- Create: `references-local.bib`

**Step 1: Add verified references**

These references are needed across the paper. Every entry MUST be verified against an authoritative source (DOI, publisher page). Do not generate from training data.

Required references (verify each before adding):
- Alavi, Malde, Schwenk, Erdos (1987) -- the original conjecture
- Chudnovsky & Seymour (2007) -- real-rootedness of tree independence polynomials
- Reynolds companion paper -- cite as GitHub repo (no arXiv yet)
- Alon, U. (2007) -- network motifs review (Nature Reviews Genetics)
- Ascoli et al. (2007) -- NeuroMorpho.org database paper
- Hinchliff et al. (2015) -- Open Tree of Life
- Hosoya (1971) -- topological index (Bull. Chem. Soc. Jpn.)
- Merrifield & Simmons (1989) -- topological methods in chemistry (book)
- Scott & Sokal (2005) or similar -- hard-core model / lattice gas reference
- Barabasi & Albert (1999) -- scale-free networks (if needed for PPI discussion)

**Step 2: Build to verify bib parses**

Run: `make`
Expected: Full build succeeds, bibliography section appears (empty until citations added to text).

**Step 3: Commit**

```bash
git add references-local.bib
git commit -m "add initial bibliography entries"
```

---

### Task 3: Draft Introduction

**Files:**
- Modify: `main.tex` (Introduction section)

**Step 1: Write introduction (~400 words)**

Structure:
1. **Opening hook** (1 para): Trees are ubiquitous in biology. Give three concrete examples: phylogenies, dendritic arbors, regulatory network backbones. Note that the combinatorial structure of these trees encodes functional properties.
2. **The independence polynomial** (1 para): Define $\indpoly$ in accessible terms. A vertex is "occupied" or "free"; an independent set is a collection of occupied vertices with no two adjacent. The polynomial counts these by size. Explain what unimodality means (single-peaked distribution) and state the conjecture (Alavi et al. 1987, Erdos Problem 993).
3. **Why biologists should care** (1 para): In a protein interaction network, an independent set is a set of non-interacting proteins. In a dendritic arbor, it's a set of non-adjacent compartments. The shape of the independence polynomial (mode, near-miss ratio) characterises the combinatorial robustness of the tree's topology.
4. **This paper's contribution** (1 para): We apply recent results from [companion paper] to real biological trees. We compute independence polynomials for 25+ neuronal reconstructions and phylogenetic trees, confirming unimodality and log-concavity in all cases. We interpret the Hub Exclusion Lemma, the hard-core edge bound, and the near-miss ratio in biological terms.

House style reminders:
- ~60 word paragraphs, max 100
- Contractions OK
- No throat-clearers
- Concrete before abstract
- `\textcite{}` for narrative citations, `\citep{}` for parenthetical
- En-dashes with spaces for parenthetical asides (`~-- text~--`)

**Step 2: Build**

Run: `make quick`
Expected: Compiles, introduction renders.

**Step 3: Commit**

```bash
git add main.tex
git commit -m "draft introduction"
```

---

### Task 4: Draft Mathematical Background

**Files:**
- Modify: `main.tex` (Section 2)

**Step 1: Write Section 2 (~800 words across 3 subsections)**

This section presents the math for a biology audience. Keep notation minimal. Define everything. Avoid proofs (point to companion paper).

**2.1 Independence polynomials and unimodality (~250 words)**
- Define: graph, tree, independent set, independence number $\alpha(T)$
- Define: $I(T;\,x) = \sum_{k=0}^{\alpha} i_k x^k$
- State unimodality conjecture (Alavi et al. 1987)
- Note: verified computationally for all trees on $\le 26$ vertices (447.7 million trees)
- Note: Chudnovsky & Seymour (2007) proved all roots are real, which gives log-concavity "for free" once unimodality is established

**2.2 The hard-core model (~250 words)**
- Define: uniform distribution over all independent sets (fugacity $\lambda = 1$)
- Define: occupation probability $\occ(v)$ = probability vertex $v$ appears in a random independent set
- State: $\mu = \sum_v \occ(v) = I'(1)/I(1)$ is the expected independent set size
- State edge bound: $\occ(u) + \occ(v) < 2/3$ for any adjacent $u, v$ (Theorem from companion paper)
- Interpret: $\{v : \occ(v) > 1/3\}$ is automatically an independent set

**2.3 Structural reduction: the 1-Private framework (~300 words)**
- Define: private neighbour, 1-Private maximal independent set
- State Hub Exclusion Lemma: if $d_{\mathrm{leaf}}(v) \ge 2$, then $v \notin S$ for any 1-Private MIS $S$, and all leaf-children of $v$ are in $S$
- State Transfer Lemma: pruning hub-leaf clusters preserves 1-Private structure
- Define near-miss ratio $\nm(T)$: measures how close $T$ is to violating unimodality
- State asymptotic: $\nm(s) = 1 - C/s + O(1/s^2)$

**Step 2: Build**

Run: `make quick`
Expected: Math renders correctly.

**Step 3: Commit**

```bash
git add main.tex
git commit -m "draft mathematical background section"
```

---

### Task 5: Draft Section 3 (Biological Interpretation)

**Files:**
- Modify: `main.tex` (Section 3)

**Step 1: Write Section 3 (~600 words across 3 subsections)**

This is the paper's main conceptual contribution: translating the math into biology.

**3.1 Neuronal dendritic arbors (~200 words)**
- Dendritic trees are literal trees (SWC format encodes parent-child relationships)
- Each compartment is a vertex; adjacency = physical connection
- An independent set = a collection of non-adjacent compartments (could model non-interacting functional units)
- Hub Exclusion predicts: high-branching points (soma, primary branch points) are excluded from maximal independent configurations; terminal dendrites are included
- This matches the biological intuition that terminal arbors are the "active" sites (synaptic input) while branch points are structural

**3.2 Phylogenetic trees (~200 words)**
- Phylogenies are trees by definition
- Tip nodes = extant species; internal nodes = ancestral divergence events
- Independence polynomial encodes the combinatorial structure of evolutionary branching
- Polytomies (unresolved branching, high degree) are exactly the hubs where Hub Exclusion applies
- The near-miss ratio characterises how "robust" the phylogeny's combinatorial structure is to perturbation

**3.3 Transcription regulatory motifs (~200 words)**
- Single-Input Modules (SIMs) in transcription networks are star-shaped trees (Alon 2007)
- One transcription factor regulates $N$ target genes
- The independence polynomial of a SIM is easily computed
- Hub Exclusion says: the regulator is excluded from any 1-Private MIS; the targets are included
- This mirrors the biological logic: targets can be independently active, but the regulator's state is determined by its targets
- The edge bound $\occ(\text{TF}) + \occ(\text{target}) < 2/3$ constrains co-expression

**Step 2: Build**

Run: `make quick`

**Step 3: Commit**

```bash
git add main.tex
git commit -m "draft biological interpretation section"
```

---

### Task 6: Generate Results Table and Figure

**Files:**
- Create: `figures/` directory
- Create: `generate_figures.py` (or extend `analyze_bio.py`)

**Step 1: Generate LaTeX table of computational results**

Write a script that outputs a LaTeX `booktabs` table with columns:
- Tree name (abbreviated), Source (NM/OToL), Species/Clade, n, $\alpha$, mode, nm

Include all 25+ trees. Group by source (neuromorpho first, then phylogenies).

**Step 2: Generate near-miss ratio vs. n scatter plot**

Plot nm on y-axis vs. n on x-axis. Two series: neuronal arbors (filled circles) and phylogenies (open triangles). Log scale on x-axis. Include reference line showing the 1 - C/s asymptotic envelope.

Save as PDF: `figures/nm_vs_n.pdf`

**Step 3: Generate independence polynomial shape figure**

Pick 3-4 representative trees (one small phylogeny, one small neuron, one large neuron, one large phylogeny). Plot the normalised coefficient sequence $i_k / \max(i_k)$ vs. $k/\alpha$ for each, overlaid. This shows the universal unimodal shape.

Save as PDF: `figures/poly_shapes.pdf`

**Step 4: Build and verify figures render**

Run: `make`
Expected: Figures appear in PDF.

**Step 5: Commit**

```bash
git add figures/ generate_figures.py main.tex
git commit -m "add results table and figures"
```

---

### Task 7: Draft Section 4 (Computational Results)

**Files:**
- Modify: `main.tex` (Section 4)

**Step 1: Write Section 4 (~600 words across 4 subsections)**

**4.1 Data and methods (~150 words)**
- Data sources: NeuroMorpho.org (15 SWC files, 5 species, 5 cell types, 113--17,992 nodes), Open Tree of Life (12 Newick files, 46--30,310 nodes)
- Software: `indpoly.py` for polynomial computation, `bio_trees.py` for format conversion
- Computation: all polynomials computed exactly (arbitrary-precision integer arithmetic)
- Runtime: sub-second for n < 1000; ~19 minutes for n = 17,992
- Cite companion paper GitHub repo for the algorithm

**4.2 Unimodality and log-concavity (~100 words)**
- All 25+ trees are unimodal AND log-concave
- This extends the exhaustive verification (n $\le$ 26, abstract trees) to real biological trees of much larger size
- Reference the results table

**4.3 Near-miss ratio and robustness (~200 words)**
- nm ranges from 0.770 (Homininae, n=46) to 0.999 (human pyramidal, n=17,992)
- Clear positive correlation between n and nm (reference scatter plot)
- Consistent with the $1 - C/s$ asymptotic from the companion paper
- Interpretation: larger biological trees sit closer to the combinatorial phase transition, but never cross it
- The "dilution effect": adding nodes smooths the independence distribution

**4.4 Hub exclusion in biological networks (~150 words)**
- Phylogenies have many hub vertices (polytomies): up to max $d_{\mathrm{leaf}} = 176$ (Papilionidae)
- Hub Exclusion Lemma is heavily active in phylogenies; less so in neurons (mostly degree-2 chains)
- Report hub counts for representative trees
- Biological meaning: in phylogenies, unresolved radiations (polytomies) are structurally excluded from maximal independent configurations

**Step 2: Build**

Run: `make`

**Step 3: Commit**

```bash
git add main.tex
git commit -m "draft computational results section"
```

---

### Task 8: Draft Discussion

**Files:**
- Modify: `main.tex` (Section 5)

**Step 1: Write Discussion (~500 words across 3 subsections)**

**5.1 The independence polynomial as a network descriptor (~200 words)**
- The independence polynomial captures combinatorial structure that scalar measures (degree distribution, clustering coefficient) miss
- The mode position, near-miss ratio, and occupation probabilities provide a richer characterisation
- For biological networks with tree-like backbone, these are exactly computable (polynomial time)
- Contrast with general graphs where computing the independence polynomial is #P-hard

**5.2 Universality: physics, chemistry, biology (~150 words)**
- The same polynomial appears in statistical physics (hard-core lattice gas), chemical graph theory (Merrifield-Simmons index), and now biological network analysis
- In chemistry, the total count $I(T; 1)$ correlates with thermodynamic properties of acyclic hydrocarbons (Hosoya 1971, Merrifield & Simmons 1989). [Source-grounding flag: verify specific correlations before writing.]
- The full polynomial $I(T; x)$ carries richer information than the single number $I(T; 1)$
- The structural results (hub exclusion, edge bound, near-miss asymptotics) are new contributions that apply across all these domains

**5.3 Limitations and open questions (~150 words)**
- The unimodality conjecture remains unproved (though verified for 447.7M abstract trees and now 25+ biological trees)
- Real biological networks are rarely exact trees; extending to near-tree structures (bounded treewidth) is an open direction
- The occupation probabilities $\occ(v)$ deserve systematic study as centrality measures
- Larger-scale analysis (all of NeuroMorpho.org's 260K+ reconstructions) is computationally feasible and would strengthen the empirical base

**Step 2: Build**

Run: `make`

**Step 3: Commit**

```bash
git add main.tex
git commit -m "draft discussion section"
```

---

### Task 9: Draft Abstract and Conclusion

**Files:**
- Modify: `main.tex` (Abstract and Section 6)

**Step 1: Write Conclusion (~150 words)**

Summarise: we applied independence polynomial theory to biological trees. All 25+ trees are unimodal and log-concave. The Hub Exclusion Lemma, hard-core edge bound, and near-miss ratio have natural biological interpretations. The independence polynomial is a promising network descriptor for tree-structured biological systems.

**Step 2: Write Abstract (~200 words)**

The abstract should cover: (1) what independence polynomials are, (2) the biological motivation, (3) what we computed (25+ real trees, all unimodal), (4) the key interpretive results (hub exclusion as motif logic, near-miss as robustness, edge bound as co-occurrence constraint), (5) the universality claim (physics, chemistry, biology).

**Step 3: Build full paper**

Run: `make`
Expected: Complete paper compiles with all sections.

**Step 4: Commit**

```bash
git add main.tex
git commit -m "draft abstract and conclusion"
```

---

### Task 10: Polish and Proofread

**Files:**
- Modify: `main.tex`

**Step 1: Run house style check**

Run: `python .house-style/check-style.py main.tex` (if available) or manually audit for:
- Semantic macros used correctly (`\term{}`, `\mention{}`, `\enquote{}`)
- En-dashes with spaces, not em-dashes
- No `\paragraph{}` headings
- No throat-clearers
- Paragraph lengths within bounds
- `\newpage` before `\printbibliography`

**Step 2: Verify all citations resolve**

Run: `make` and check for undefined citation warnings.

**Step 3: Review source grounding**

Audit every factual claim:
- Are statistics from our actual computations? (Check against `analyze_bio.py` output.)
- Are citations verified? (Especially Merrifield-Simmons, Hosoya, Alon.)
- Any claims that need sources but don't have them?

**Step 4: Final build**

Run: `make clean && make`
Expected: Clean build, no warnings.

**Step 5: Commit**

```bash
git add main.tex references-local.bib
git commit -m "polish: house style, citations, source grounding"
```

---

## Dependencies Between Tasks

```
Task 1 (structure) ─┬─> Task 2 (bib) ─────> Task 3 (intro) ─> Task 4 (background)
                     │                                              │
                     └─> Task 6 (figures) ──────────────────────────┤
                                                                    v
                                                           Task 5 (bio interpretation)
                                                                    │
                                                                    v
                                                           Task 7 (results) ──> Task 8 (discussion)
                                                                                        │
                                                                                        v
                                                                               Task 9 (abstract/conclusion)
                                                                                        │
                                                                                        v
                                                                               Task 10 (polish)
```

Tasks 2 and 6 can run in parallel after Task 1. Tasks 3--5 are sequential (each builds on the previous). Task 6 (figures) can be done any time after Task 1 but must be done before Task 7.
