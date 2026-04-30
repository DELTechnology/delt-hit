# Revision Plan — DELT-Hit Manuscript
*Ordered trivial → major. Last updated: 2026-04-29.*

---

## Trivial — pure manuscript text edits (no new content, no code dependency)

**T2. Competing-tools tone** (para 63, and para 33–34)
- Para 63: "the only other open-source DEL framework introduced to date is DELi" → "another open-source DEL framework is DELi, which addresses several core aspects of DEL informatics…"
- Para 33–34: soften or remove any remaining claim that "no general, open-source software solution is available".

**T3. Figure 1 caption clarification** (para 39)
Add one sentence distinguishing data flow from processing steps: "Solid arrows indicate data flow between modules; blue rectangles represent processing steps; grey boxes indicate input/output file types."

**T4. Replicate and sequencing depth guidance** (Experimental Design > Quality control, ~para 88–92)
Add 2–3 sentences: "We recommend a minimum of three biological replicates per condition for edgeR-based analysis. For reliable hit identification, a sequencing depth of at least 10–20 reads per library member per selection is advisable."

**T5. Error-handling balance sentence** (demultiplexing section, ~para 84/279/491)
Add: "DELT-Hit balances sensitivity and specificity by allowing per-region configuration of maximum error rates and indel tolerance, so users can apply stricter matching to constant regions while permitting controlled mismatches in barcode regions."

---

## Small — new content blocks (self-contained, no code dependency)

**S1. Abstract additions** (para 9)
Three additions at the end of the abstract:
- Skill-level sentence: "The protocol is suitable for users with basic command-line experience; no prior expertise in cheminformatics or bioinformatics is required."
- Time estimate: "Depending on library size, the complete workflow takes 2–8 hours, with chemical enumeration being the most time-consuming step, scaling with library size."
- GitHub link: "DELT-Hit is freely available at https://github.com/DELTechnology/delt-hit."

**S2. Minimal YAML example inline** (Experimental Design > Input requirements, ~para 71–74)
Insert a short code block showing the minimal `analysis.yaml` structure (experiment name, fastq_path, save_dir, selections with group assignments).

**S3. Inputs/outputs table per module** (Experimental Design or Materials)
New table: Module | Key input(s) | Key output(s). One row per module: init, demultiplex, library enumerate, library properties, analyse, dashboard.

**S4. Quick-start overview** (after "Overview of the procedure", ~para 44–49)
A compact 5-row summary table: Stage | CLI command | Primary output. Maps each phase to the one command a user runs and what file comes out.

**S5. QC figure caption guidance** (all QC figure captions)
For each QC figure caption, add 1–2 interpretation sentences describing what a healthy output looks like vs. a warning sign.

**S6. Worked chemical example** (library chemistry section, ~para 82–84)
A small table or code block: one reaction (e.g. amide coupling), its SMIRKS string, one building block SMILES, and the resulting product SMILES, showing the enumeration logic concretely.

**S7. Low read retention troubleshooting expansion** (~para 533)
The existing "<50% retention" bullet is too brief. Expand to 4–5 concrete causes and fixes:
1. Primer / constant region sequence mismatch in the structure sheet.
2. Wrong barcode orientation (try reverse complement).
3. Max error rate set too low for the sequencing quality.
4. FASTQ quality issues (check Q-score distribution before demultiplexing).
5. Library structure sheet does not match the actual construct design.

---

## Medium — new content that requires describing new/updated CLI behaviour

**M1. Enumeration validation workflow** (~para 84, library chemistry CRITICAL STEP)

Replace the simple "validate on a small subset" note with a proper two-step workflow that reflects the current tooling:

> **CRITICAL STEP — Validating enumeration before full library runs**
>
> Before running enumeration on the full library, we strongly recommend the following validation workflow:
>
> **Step 1 — Visual inspection with `delt-hit visualize enumerate`**
> Run `delt-hit visualize enumerate` to generate reaction graphs, building-block grids, 2D reaction scheme panels, and example compound structures. This lets you visually confirm that building blocks, reaction SMIRKS, and the reaction graph are correctly configured before committing to a full enumeration run.
>
> **Step 2 — Debug enumeration run**
> Run the full enumeration with `--debug invalid --errors raise`. This stops enumeration at the first compound that cannot be successfully constructed and renders the reaction graph that was attempted, providing a precise entry point for diagnosing SMIRKS or connectivity issues.
>
> **On reaction product uniqueness:**
> DELT-Hit deliberately requires each reaction to produce exactly one well-defined product. Unlike some other frameworks, multiple distinct reaction outcomes are not supported. If enumeration reports multiple products for a given reactant combination, this indicates that the SMIRKS definition is too broad and needs to be tightened to a more specific pattern. This design choice reflects the expectation that DEL synthesis routes are well-defined chemical transformations, and that ambiguity in the product is a signal of an incorrect or overly general reaction definition rather than a feature to be silently handled.

*Note for rebuttal letter: explicitly contrast this design decision with DELi's `pick_fragment` / `ignore_multiple_products` approach and argue why deterministic, single-product reactions are the right constraint for DEL synthesis.*

---

## Major — require data, experiments, or significant new sections (not yet unblocked)

**X1. Re-analysis of NF2 dataset (R1.1)**
Re-run Nat Chem 2021 library with DELT-Hit; produce comparison table of counts (original tool vs. DELT-Hit). Add as a new subsection under "Applications of the method".
*Blocked on: Jörg providing FASTQ + original counts files.*

**X2. Quantitative DELi benchmark (R1.2 / R2.2)**
Run both DELT-Hit and DELi on the synthetic libraries using the new benchmarking infrastructure; collect timing, demultiplexing stats, and enrichment reliability numbers. Update / replace the "Comparison with other methods" section with quantitative data.
*Blocked on: Leandro setting scoring; benchmark runs collecting results.*

**X3. edgeR justification with positive-control retrieval (R1.3 / R3.3)**
Run edgeR and alternative scoring methods on NF library with triplicates; show positive-control compound recovery rates. Add as a new subsection or supplementary figure.
*Blocked on: Jörg providing NF library files + Johanna compiling positive-control list.*

**X4. Normalised Z-score enrichment mode (R4.4)**
Implement `--method zscore` in `delt-hit analyse enrichment`; add CLI documentation and a protocol step describing when to use it (single-replicate experiments). Full design spec in `issues/z-score.md`.
*Blocked on: Leandro implementing the feature.*

**X5. Dual-display case study (R3.2)**
Add GB library dual-display example to "Applications of the method" with enumeration and enrichment results.
*Blocked on: Alice producing the GB library pairing Excel; Adriano integrating dual-display into package.*

**X6. SMIRKS template library expansion (R4.1)**
Add a supplementary file with validated DEL-compatible reaction SMIRKS (with 2D chemical schemes). Reference from the library chemistry section.
*Blocked on: Alice compiling the file.*

**X7. Diverse selection conditions discussion (R4.3)**
Add a paragraph to Experimental Design explaining how the analysis modules handle competitive experiments, concentration gradients, and varied stringency washes.
*Blocked on: Adriano investigating and writing.*

---

## Rebuttal letter notes (not manuscript edits, but document here for drafting)

- **R3.1 / competing tools:** Explain that DELi exists but DELT-Hit is differentiated by its fully integrated, end-to-end workflow (single Excel config → enumeration → demultiplexing → enrichment → dashboard). DELi is stronger as a library-definition model; DELT-Hit is stronger as a complete analysis pipeline.
- **M1 / single-product design:** Contrast with DELi's `pick_fragment` / `ignore_multiple_products` / `fail_on_multiple_products` flags. Argue that for well-defined DEL synthesis, reaction ambiguity is always a SMIRKS specification problem, not a runtime choice. Single-product enforcement catches errors early and produces auditable enumeration outputs.
- **R4.5 / filtered enumeration:** Already implemented (`--counts_path`, `--top_n`). Show in rebuttal.
