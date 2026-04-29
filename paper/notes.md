# Paper Notes — DELT-Hit Manuscript & Revision

*Reference notes synthesised from `Martinelli_etal_manuscript.docx`, `260422_Nature Protocols_Revisions.docx`, and `CHANGELOG.md` (2026-04-24 to 2026-04-29). Last updated: 2026-04-29.*

---

## 1. Manuscript at a glance

**Title (current):** "DELT-Hit: An end-to-end computational framework for DNA-encoded chemical library analysis"
→ Editor asks to remove punctuation from the title.

**Target journal:** Nature Protocols (under major revision)

**Authors:** Adriano Martinelli†*, Alice Lessing†, Gary Hoppeler, Andreas Gloger, Louise Plais, Jörg Scheuermann* (†equal contribution, *corresponding)
→ Leandro and Johanna need to be added to the author list.

**One-sentence pitch:** DELT-Hit is an open-source Python CLI that converts raw DEL FASTQ reads into ML-ready chemical information through five interconnected modules: (1) adaptive demultiplexing, (2) chemical structure reconstruction, (3) molecular property/descriptor generation, (4) statistical hit ranking, and (5) QC & visualisation.

---

## 2. Revision action items by owner

### Jörg
- [ ] Take over Polybox and ensure all authors have access.
- [ ] Provide AM with FASTQ files, Excel files, and original counts files for the NF2 library used in the Nat Chem 2021 publication.
- [ ] Provide AM with files and info on positive controls (compound IDs etc.) for NF library edgeR validation.
- [ ] Once revision is ready: make the GitHub repository publicly available.

### Alice
- [ ] Improve Procedure formatting per Nature Protocols style (see editorial notes §2.3 below).
- [ ] Compile a pairing-library Excel file for the GB library (dual-display architectures case study).
- [ ] Compile a file with all current reaction SMIRKS (with chemical schemes) and add basic templates users can customise.
- [ ] Add reaction visualisation to the SMIRKS overview file.

### Leandro
- [ ] Set up scoring for DELi vs. DELT-Hit benchmark.
- [ ] Look into the normalised Z-score methodology (Faver et al. 2019, DOI: 10.1021/acscombsci.8b00116) as a third enrichment mode (see §4.4 below).

### Johanna
- [ ] Compile a list of positive control building blocks and their correspondence to proteins.

### Adriano
- [ ] Change compound indexing to be consistent throughout.
- [ ] Re-run NF2 library with DELT-Hit and compare counts to original results from Nat Chem 2021.
- [ ] Run both DELi and DELT-Hit on synthetic libraries; compare stats, compare speed; speed up enumeration.
- [ ] Compare edgeR vs. other scoring methods on NF library to show robust retrieval of positive controls.
- [ ] Integrate dual-display handling into the current package.
- [ ] Implement molecule visualisation (2D structures of top-N final compounds, reaction graph in SMIRKS overview).
- [ ] Look into handling of nuanced selection conditions (competitive experiments, concentration gradients, varied stringency washes).
- [ ] Look into "filtered enumeration" mode (ingest counts matrix, restrict SMILES/descriptor computation to compounds above a read threshold). **Status: implemented** — it is now possible to enumerate selected compounds (e.g. top 1000) out of raw counts or enrichment analysis files.

---

## 3. Editorial requests (Nature Protocols editor)

### 3.1 Administrative / code availability
- Obtain a DOI for the GitHub repo via Zenodo and cite it in the manuscript and Code Availability statement.
- Confirm availability of the web-based tool for the foreseeable future in Code Availability.
- All corresponding authors must link ORCIDs in the Manuscript Tracking System before acceptance.

### 3.2 Abstract
- Add a sentence indicating the required skill-set / experience level.
- Add a time estimate (benchmark across library sizes; note enumeration is the most time-consuming step and scales with library size).
- Include a link to the GitHub repository.

### 3.3 Technical Overview section
- This section does not conform to Nature Protocols format. Move content into "Experimental design" and/or the Materials > Software heading.

### 3.4 Materials
- Equipment section may have subheadings; software installation/setup can be a subheading under "Experimental Setup".
- Useful information can be formatted as `<CRITICAL>` callouts.

### 3.5 Procedure formatting
- Must be a continuously numbered sequence (1, 2, 3 … not 1.1, 1.2).
- Only one level of subheading within the Procedure.
- Options within a step use A, B, C; sub-steps within options use roman numerals.
- Optional procedures that break the linear flow go in boxes, referred to from the main steps.
- All explanatory text should be either inside a numbered step or inside a callout (CAUTION, CRITICAL STEP). Consider moving explanatory text to the Introduction / Experimental Design.

### 3.6 Timing
- Add corresponding Procedure step numbers to each entry.

### 3.7 Troubleshooting
- Add the word "Troubleshooting" after each relevant step in the Procedure.
- Add a "Step" column to the Troubleshooting table (the step where the problem becomes evident).

### 3.8 Figures
- All figures are new (not previously published).
- Screenshots: max 17 cm width; submit as TIFF, do not resize before saving; crop to relevant content.
- Follow NPG image manipulation guidelines.

### 3.9 Supplementary Information
- Use the Nature Protocols SI inventory template.
- Extended Data: ≤10 figures, each on one A4 page with legend; submitted as individual JPG/TIF/EPS.
- Supplementary figures and textual content in a single combined PDF, numbered "Supplementary Figure 1" etc.
- Source Data: full-length blots/gels as PDF; data underlying graphical figures in a single Excel file with one tab per figure.

---

## 4. Reviewer comments & required actions

### 4.1 Reviewer 1 (major revision)

**R1.1 — Validation on cited case studies**
The "Applications" section cites prior group publications (e.g. Favalli et al. Nat. Chem. 2021), but those analyses used different tools. Required: re-analyse at least one published dataset with DELT-Hit and provide a comparison table of key results (original method vs. DELT-Hit). Show that DELT-Hit enables additional analysis/evaluation steps.
→ *Action:* Adriano re-runs NF2 library; Jörg provides data files.

**R1.2 — Comprehensive DELi comparison**
The current "Comparison with other methods" section is not detailed enough. Required: expand to include quantitative performance comparisons on the same benchmark datasets.
→ *Action:* Run both tools on synthetic libraries; compare technical stats and user-friendliness. Leandro sets scoring; Adriano runs benchmarks and works on enumeration speed.

**R1.3 — Justification of edgeR for DEL data**
edgeR was developed for RNA-seq; DEL data can be sparse and may lack biological replicates. Required: additional discussion and/or supporting data justifying edgeR's applicability.
→ *Action:* Show on NF library that different scoring methods (especially edgeR) enable robust retrieval of positive controls. Ideally use a dataset with triplicates. Johanna compiles positive control info; Adriano runs the comparison.

**R1.4 — Language and formatting**
Careful proofreading for consistency in terminology and readability.

### 4.2 Reviewer 2 (minor revision)

**R2.1 — Timing per application example**
Table 14.1 gives timing estimates for different dataset sizes, but the "Applications" examples don't show actual runtimes. Annotate each case with approximate time and address whether runtime scales non-linearly with library size or architectural complexity.
→ *Action:* Covered by the DELi benchmark work (R1.2).

**R2.2 — Quantitative DELi comparison**
Same as R1.2: add numerical comparisons of enrichment reliability and replicate consistency on shared benchmarks.

### 4.3 Reviewer 3 (minor revision / clarifications)

**R3.1 — Tone re competing tools**
Replace strong "no existing open-source solution" claims with a nuanced statement that acknowledges DELi while emphasising DELT-Hit's integrated, end-to-end nature.

**R3.2 — Dual-display case study**
Explicitly state that case studies demonstrate robustness across diverse DEL architectures, including dual-display.
→ *Action:* Alice produces the GB library pairing Excel; Adriano integrates dual-display handling into the package.

**R3.3 — edgeR statistical methodology** — same as R1.3.

**R3.4 — Error handling**
Include a summary sentence explaining the balance between strict barcode matching and controlled error tolerance.

**R3.5 — Enumeration validation guidance**
Encourage users to validate reaction definitions on small subsets before full enumeration. Stress this point clearly in the manuscript.

**R3.6–3.13 — Minor additions** (each is a small manuscript edit):
- Quick-start overview summarising workflow steps.
- Table of inputs and outputs per module.
- Improve Figure 1 to distinguish data flow from processing steps.
- Expand QC figure captions with interpretation guidance.
- Provide a minimal YAML example inline.
- Clarify recommended number of replicates and sequencing depth.
- Include a small worked chemical example.
- Add more troubleshooting guidance for low read retention.

### 4.4 Reviewer 4 (major revision)

**R4.1 — Expand SMIRKS template library**
"Defining SMIRKS requires advanced cheminformatics expertise" is flagged as a barrier. Required: expand the template library with more validated, DEL-compatible chemistries.
→ *Action:* Alice compiles file with all current reactions (with chemical schemes) plus additional basic SMIRKS.

**R4.2 — 2D chemical visualisation**
Abstract chemical formats (SMILES strings, schematic graphs) are hard for chemists. Required: display 2D reaction schemes alongside SMIRKS strings; include 2D structural diagrams for building blocks and products from the tutorial dataset.
→ *Action:* Adriano implements molecule visualisation (top-N final compounds); Alice adds reaction visualisation to SMIRKS overview file.

**R4.3 — Diverse selection conditions**
The tool primarily demonstrates "protein vs. no-protein". Required: discuss how the analysis modules handle competitive experiments, concentration gradients, or varied stringency washes.
→ *Action:* Adriano investigates and adds to manuscript.

**R4.4 — Normalised Z-score enrichment mode**
"Simple Counts" lacks statistical rigour for single-replicate experiments. Required: implement the normalised Z-score method (Faver et al. 2019, ACS Comb. Sci., DOI: 10.1021/acscombsci.8b00116) as a third enrichment mode using a binomial distribution model.
→ *Action:* Leandro leads; full design documented in `issues/z-score.md`. Key formula:
  - Standard binomial z-score: `z = (p_o − p_i) / sqrt(p_i · q_i / n)`
  - Normalised (depth-insensitive): `z_n = (p_o − p_i) / sqrt(p_i · q_i)`
  - where `p_i = 1 / library_size`, `p_o = count / total_reads`, `n` = sequencing depth.

**R4.5 — Filtered enumeration mode**
Full in-silico enumeration of billion-member libraries is expensive. Required: optional mode that restricts SMILES construction and descriptor calculation to barcode combinations above a user-defined read threshold.
→ *Status: implemented.* Mention in manuscript.

---

## 5. Major revision checklist — status against CHANGELOG

The CHANGELOG covers work merged between 2026-04-24 and 2026-04-29. Each item below is one required revision point; status reflects what the changelog confirms as done vs. what still needs manuscript or experimental work.

| # | Requirement | Source | Code status | Manuscript / data still needed? |
|---|---|---|---|---|
| R1.1 | Re-analyse NF2 (Nat Chem 2021) dataset with DELT-Hit; provide comparison table | R1, R2 | CHANGELOG mentions "additional experiment templates and supporting-material runs for Favalli and related comparison workflows" — infrastructure is in place | Yes — actual counts comparison and results table must still be produced and written up |
| R1.2 / R2.2 | Quantitative DELi benchmark on synthetic libraries; compare demultiplexing stats, enrichment reliability, speed | R1, R2 | CHANGELOG: "Synthetic benchmarking utilities for demultiplexing and enrichment, including FASTQ generation, DELT/DELi input converters, runtime plotting, worker/chunk-size sweeps, and large-library benchmark matrices" — tooling done | Yes — results need to be collected, compared, and written into the "Comparison" section |
| R1.3 / R3.3 | Justify edgeR for DEL data; show positive-control retrieval on NF library with triplicates | R1, R3 | Not mentioned in CHANGELOG | Blocked — Jörg still needs to provide NF library files and positive-control info (Johanna) |
| R3.2 / R4.2-vis | Dual-display case study + 2D chemical visualisation (reaction graphs, BB grids, compound structure panels) | R3, R4 | CHANGELOG: "New `visualize enumerate` CLI workflow … including reaction graphs, reaction schemes, building-block grids, and compound structure panels" — feature implemented | Yes — GB library pairing Excel still needed from Alice; dual-display integration confirmation needed; screenshots/figures to be prepared |
| R4.1 | Expand SMIRKS template library with validated DEL-compatible chemistries + chemical schemes | R4 | Not explicitly mentioned in CHANGELOG | Alice still needs to compile the file |
| R4.3 | Discuss handling of nuanced selection conditions (competitive, gradient, stringency) | R4 | Not in CHANGELOG | Manuscript-only addition; Adriano to investigate and write |
| R4.4 | Implement normalised Z-score enrichment (Faver et al. 2019) as third enrichment mode | R4 | CHANGELOG mentions "z-score" under Documented — design documented but not yet implemented | Leandro to implement; CLI flag `--method zscore`; full spec in `issues/z-score.md` |
| R4.5 | Filtered enumeration from counts matrix | R4 | CHANGELOG: "Filtered library enumeration from observed `counts.txt` inputs with `--counts_path`, `--top_n`, and `--library_name`" — **implemented** ✓ | Mention in manuscript (one sentence in Procedure/Features) |
| R2.1 | Annotate each application example with actual runtime | R2 | CHANGELOG: benchmarking infrastructure and "configurable loops, estimated iteration controls" — tooling done | Yes — timing numbers for each case study must be collected and added to Table 14.1 / application examples |
| Editorial | Title punctuation removed | Editor | N/A | Manuscript edit only |
| Editorial | Abstract: skill-set sentence, time estimate, GitHub link | Editor | Time-estimate data may come from benchmarks (see R2.1) | Manuscript edit; GitHub link needs Zenodo DOI first |
| Editorial | Technical Overview section reformatted | Editor | N/A | Manuscript restructuring (Alice) |
| Editorial | Procedure renumbered continuously | Editor | N/A | Manuscript restructuring (Alice) |
| Editorial | Timing section: add step numbers | Editor | N/A | Manuscript edit |
| Editorial | Troubleshooting: add "Step" column | Editor | N/A | Manuscript edit |
| Editorial | Zenodo DOI for GitHub repo | Editor | N/A | Jörg / Adriano to create once repo is public |
| R3.1 | Soften claims about "no existing open-source solution" | R3 | N/A | Manuscript edit only |
| R3.4 | Add sentence on barcode error-tolerance balance | R3 | N/A | Manuscript edit only |
| R3.5 | Stress enumeration subset-validation guidance | R3 | N/A | Manuscript edit only |
| R3.6–3.13 | Minor additions (quick-start, I/O table, Figure 1 revision, YAML example, QC captions, replicate guidance, worked example, low-read-retention troubleshooting) | R3 | CHANGELOG: "reviewer-oriented examples and workflow notes" — partial | Most are manuscript edits; Figure 1 needs graphic work |

**Summary:** Of the ~20 required revision points, the CHANGELOG confirms that code-side work is complete or substantially done for **4 items** (filtered enumeration ✓, visualise-enumerate CLI ✓, synthetic benchmarking tooling ✓, Favalli experiment templates ✓). The remaining items are either blocked on data from Jörg/Johanna, still need Alice's file compilation, require Leandro's Z-score implementation, or are manuscript-only edits.

---

## 6. Open engineering issues (from `issues/`)

### 5.1 `get_counts` performance (`issues/get_counts_performance.md` + `…_plan.md`)
- Bottleneck: serial gzip decompression in the main process (not Python parsing).
- Current state: parallel path with `isal` + Python bytes ops achieves ~7× speedup on 5 M reads / 48-core node.
- Cython `parse_chunk` was implemented but gives no net gain because decompression is the true bottleneck.
- Next steps documented in `get_counts_performance_plan.md`:
  1. Keep isal + parallel Python as the production path.
  2. Unify the serial path (`num_workers=1`) to also use isal.
  3. Optionally explore Option C (pre-split FASTQ) or Option D (Polars/Arrow) for billion-read scale.

### 5.2 Chemical-properties API alignment (`issues/chemical-properties.md`)
- `library properties` should accept the same inputs as the enumeration API (`counts.txt` + library name).
- Output should be written to `properties/<name>.parquet`.

### 5.3 Enumeration UX (`issues/make_enumeration_more_intuitive.md`)
- Current model exposes the internal reaction graph; users must manually maintain `B*` sheet fields + `reaction_graph` sheet.
- Proposed: introduce a `route` column on `B*` sheets and an explicit `reaction_steps` sheet; compile these into the internal graph.
- Full roadmap in the issue file (8 phases from domain objects through stable compound IDs).

### 5.4 Normalised Z-score enrichment (`issues/z-score.md`)
- Full design spec is in the issue file.
- v1 scope: compound-level only; expected probability = `1 / library_size`; CLI flag `--method zscore`; outputs to `<save_dir>/<name>/zscore/`.
- Triggered by Reviewer 4 comment R4.4.

---

## 7. Comparison with DELi (`other_tools/report.md` summary)

| Dimension | DELT-Hit strength | DELi strength |
|---|---|---|
| Workflow integration | Single Excel workbook drives everything | Separate, reusable JSON/CSV library definitions |
| Enumeration model | Graph assembled implicitly from BB rows + `reaction_graph` | Explicit ordered reaction DAG; `bb_ids → structure` is direct |
| Ambiguity handling | Not configurable per reaction | Per-step `pick_fragment`, `fail_on_multiple_products`, alternative SMARTS lists |
| Demultiplexing | Cutadapt-based with Excel-driven structure | Built-in barcode decode with UMI support |
| Analysis / dashboard | Full pipeline: properties, enrichment, dashboard | Less central |
| Version control | Harder (binary Excel) | Good (text JSON/CSV) |

Immediate actions recommended by the report (relevant to revision):
1. Add `delt-hit library validate`.
2. Add `enumerate-one` for debugging decoded hits.
3. Extract RDKit reaction execution into a small tested module with ambiguity policies.
4. Add stable building-block IDs to enumeration outputs.
5. Add structured failure/warning outputs for enumeration.

---

## 8. Key external references to cite or address

- Favalli et al. *Nat. Chem.* 2021 — original NF2 dataset to re-analyse (R1.1).
- Faver et al. *ACS Comb. Sci.* 2019, DOI: 10.1021/acscombsci.8b00116 — normalised Z-score method (R4.4).
- Decurtins et al. *Nat. Protoc.* 2016 — group's prior Nature Protocols DEL selection paper (context for R3.1).
- GitHub repo DOI (Zenodo) — must be obtained and cited (editorial).
