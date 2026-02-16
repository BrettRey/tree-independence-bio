# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Role: Researcher / Writer

## Project Overview

This paper applies the computational and theoretical machinery developed for Erdős Problem #993 (tree independence polynomial unimodality) to biological network analysis. The companion math paper lives at `papers/Erdos_Problem_993/`.

### Motivation

Simulated reviewer "Alon" (modelling Noga Alon's perspective) noted that tree independence polynomials connect naturally to biological networks:
- **Protein interaction networks**: tree-like local structure; independent sets correspond to non-interacting protein subsets
- **Phylogenetic trees**: independence polynomials encode combinatorial properties of evolutionary trees
- **Hard-core model**: the occupancy probabilities P(v) = I'_v(1)/I(1) map to network motif centrality measures
- **Near-miss ratio**: nm as a robustness/fragility measure for network motifs

### Key connections to develop

1. **Independent sets as non-interacting modules** in protein interaction networks (tree-like backbone)
2. **Hard-core edge bound** P(u)+P(v) < 2/3 as a constraint on co-occurrence of adjacent network motifs
3. **Hub Exclusion Lemma** as a structural prediction: high-degree hubs are excluded from maximal independent configurations
4. **Near-miss ratio** as a measure of how close a tree-structured network is to combinatorial phase transition
5. **Unimodality** as a structural regularity property that biological networks inherit from their tree-like topology

### Source material

All paths relative to the portfolio root (`LLM-CLI-projects/`):

- Math paper: `papers/Erdos_Problem_993/paper/main_v2.tex`
- Computational code: `papers/Erdos_Problem_993/indpoly.py`
- Hard-core model results: `papers/Erdos_Problem_993/notes/one_private_status.md`
- Subdivision identity: `papers/Erdos_Problem_993/notes/subdivision_new_findings.md`

### Target venue

Computational biology or network science journal (e.g., PLOS Computational Biology, Journal of Theoretical Biology, or a bioinformatics venue). The framing should be accessible to biologists who work with network models.

### Author note

Brett is a linguist, not a biologist or mathematician. This paper should be honest about the interdisciplinary bridging: the math is established in the companion paper, and the biology connection is the contribution here. The paper needs real biological data or at least concrete biological examples to be credible.

## Current State

**Early skeleton.** `main.tex` has section headings (Introduction through Conclusion) with TODO placeholders. No content, no references yet. The `references.bib` is empty.

**Dependency on companion paper:** The mathematical machinery (independence polynomials, hard-core model, Hub Exclusion Lemma, near-miss ratio, subdivision identity) is developed in `papers/Erdos_Problem_993/`. This paper's contribution is the biological application, not new math. Read the companion paper's results before drafting.

**Needs real biology:** The paper won't be credible without concrete biological data or examples (protein interaction networks, phylogenetic trees). Finding or generating these is a prerequisite for serious drafting.

## Build System

This LaTeX project requires **XeLaTeX** (not pdfLaTeX) due to the Charis SIL font requirement.

**Avoid LuaLaTeX** – it tends to run words together in the underlying PDF text layer, breaking copy-paste and accessibility.

### Compilation Commands

```bash
# Full build sequence
xelatex main.tex
biber main
xelatex main.tex
xelatex main.tex

# Or use automated build
make              # Full build
make quick        # Single pass
make clean        # Clean artifacts
```

The multiple runs are necessary to resolve all cross-references and citations.

## File Structure

```
Tree_Independence_Polynomials_and_Biological_Network_Motifs/
├── main.tex                  # Main document
├── references.bib            # Bibliography
├── .house-style/             # House style snapshot
│   ├── preamble.tex         # LaTeX preamble
│   └── style-rules.yaml     # Style conventions
├── Makefile                  # Build automation
├── CLAUDE.md                 # This file
├── AGENTS.md                 # Synced with this file
└── GEMINI.md                 # Synced with this file
```

## House Style

This project uses Brett Reynolds house style (see `.house-style/style-rules.yaml`).

### Key LaTeX Conventions

**Terms, Mentions, Quotations:**
- Use `\term{}` for terms/concepts (small caps): `\term{definiteness}`
- Use `\mention{}` for mentions/forms (italics): `\mention{the}`
- Use `\olang{}` for object language (italics): `\olang{der Hund}`
- Use `\enquote{}` for quotations: `\enquote{actual quote}`
- Never use raw quotes for mention

**Cross-linguistic Notation:**
- Cross-linguistic: `\textsc{subject}\crossmark`
- Language-specific: `\textsc{subject}\textsubscript{eng}`

**Dashes:**
- Parenthetical: `foo~-- bar~-- baz` (en dash with spaces)
- Ranges: `2001--2025` (en dash, no spaces)
- Compounds: `corpus-based` (hyphen)

**Citations:**
- Parenthetical: `\citep{key}`
- Textual: `\textcite{key}`

**Citations and BibTeX (LAW):**
- Citations and BibTeX entries must NEVER be placeholders
- Citations must NEVER be generated from training data
- LLMs MUST browse the web to find DOIs and verify bibliographic data
- Every citation must be confirmed against an authoritative source
- If you cannot verify a citation, say so. Do not guess. Do not fabricate.

### Writing Style

**Preferred:**
- Use contractions (don't, won't)
- Keep paragraphs ~60 words, max 100
- Direct verbs and short clauses
- Concrete before abstract

**Avoid:**
- Throat-clearers: "It is important to note that..."
- `\paragraph{}` headings (use topic sentences)
- Bold labels in prose
- Hackneyed adverbs: moreover, furthermore

**Document Structure:**
- Use `\section{}` and `\subsection{}` only
- Avoid bullet points for arguments (use prose)
- Use ordinal markers: "first," "second," "third"

**Examples (gb4e):**
```latex
\ea\label{ex:example}
\textit{Example sentence.}
\z
```

## Common Tasks

**Adding References:**
1. Add entries to `references-local.bib` (create if needed; preamble loads it via `\IfFileExists`)
2. Protect capitals: `title = {The {Cambridge} Grammar...}`
3. Use `\textcite{key}` or `\citep{key}`
4. The central bib (`.house-style/references.bib`) is protected; use `/push-bib` to merge local entries at polish time

**Building:**
```bash
make              # Full build
make quick        # Fast build
make clean        # Clean up
```

**Git Workflow:**
- Pre-commit hook keeps CLAUDE.md, AGENTS.md, GEMINI.md synced
- Commit often with meaningful messages
- Build before committing to ensure no LaTeX errors

## Multi-Agent Dispatch (MANDATORY)

**Before dispatching multiple agents, ALWAYS ask Brett:**

1. **Which model(s)?** Options: Claude, Codex, Gemini, Copilot
   - Codex is often best for Brett's work
   - Claude has the smallest context window
   - Different models have different strengths

2. **Redundant outputs?** Should multiple models tackle the same task?
   - Useful for judgment calls (e.g., "Should I add this figure?")
   - Not needed for factual tasks

### CLI Command Patterns

| CLI | Command | Notes |
|-----|---------|-------|
| **Codex** | `codex -p 'prompt' > output.txt &` | Include "Read [PATH] first" in prompt |
| **Gemini** | `cat file.tex \| gemini --yolo -o text 'prompt'` | Must pipe content (file reading broken in YOLO) |
| **Copilot** | `copilot -p 'prompt' > output.txt &` | Fast; add `--allow-all-tools` for file ops |

**Token limits:** Gemini > Codex > Claude (most constrained)

See portfolio `CLAUDE.md` or `HPC book/.agent/workflows/multi-agent-review.md` for full patterns.
