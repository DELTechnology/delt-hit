    # NP-P250643A Revision Task List

This file turns the editor and reviewer comments for the Nature Protocols revision into a concrete checklist for the manuscript, code, figures, benchmarking, and supporting material.

## Editor

- [x] Confirm the author list and author order are correct.
- [x] Add a statement in Code Availability confirming foreseeable availability of any freely accessible web-based tool. Not applicable: DELT-Hit is not a freely accessible web-based tool.
- [ ] Make the GitHub repository citable via Zenodo and add the DOI citation to the manuscript, Code Availability statement, and reference list.
- [x] Change the title to avoid punctuation. Current title in `paper/main.tex`: `An end-to-end computational framework for DNA-encoded chemical library analysis`.
- [x] Add a sentence to the Abstract describing the required skill set or experience level. Added in `paper/main.tex`.
- [x] Add a sentence to the Abstract describing how long the protocol takes. Added in `paper/main.tex` with placeholder timing `[X-Y] h` still to be finalized.
- [x] Add a software link in the Abstract. Added in `paper/main.tex`.
- [x] Rework the Technical Overview section to match Nature Protocols format, moving content into Experimental Design and/or Materials > Software as needed.
- [x] Reformat the Procedure as one continuously numbered sequence.
- [x] Limit the Procedure hierarchy to one subheading level.
- [ ] Format options within steps as `A`, `B`, `C`, and sub-steps within options as roman numerals.
- [ ] Move optional branches that break linear flow into boxes and refer to them from the main Procedure.
- [x] Ensure explanatory text sits either inside numbered steps or in callouts such as `CRITICAL STEP` or `CAUTION`.
- [x] Add Procedure step numbers to the Timing section.
- [x] Add the word `Troubleshooting` after each relevant Procedure step.
- [x] Add a `Step` column to the Troubleshooting table.
- [ ] Check all figures against Nature Protocols image requirements and confirm that all figures are new.
- [ ] Prepare supplementary information in the required Nature Protocols format, including SI inventory, Extended Data, Supplementary Figures PDF, and Source Data files.

## Reviewer 1

- [ ] Re-analyse at least one previously published DEL dataset with DELT-Hit.
- [ ] Add a comparison table showing the original analysis versus DELT-Hit for that case study.
- [ ] Show what additional analysis or evaluation DELT-Hit enables on that dataset.
- [x] Expand the DELi comparison into a quantitative benchmark on shared datasets.
- [x] Compare technical performance, enrichment results, and practical usability against DELi. Addressed for technical performance and practical usability through the DELi comparison table plus shared runtime and memory benchmarking; quantitative enrichment reliability remains tracked separately below.
- [ ] Add discussion and/or data justifying the use of edgeR for DEL data.
- [ ] Show that edgeR and other scoring methods recover positive controls robustly.
- [ ] Preferably support the edgeR justification with a dataset containing replicates.
- [ ] **Confirm with Jörg** which selection indices in the Favalli dataset correspond to condition vs no-condition groups (currently `1_1/2_1/3_1` vs `4_1/5_1/6_1` as placeholder in `supporting_material/experiments/favalli/analysis.yaml`). Required before enrichment results for the Favalli re-analysis can be included in the manuscript.
- [ ] Proofread the manuscript for terminology, language consistency, and readability.

## Reviewer 2

- [ ] Add approximate runtimes for each application example, not just general timing tables.
- [x] Clarify how runtime scales with library size and/or architectural complexity.
- [ ] Add quantitative DELi comparisons for enrichment reliability and replicate consistency.

## Reviewer 3

- [ ] Soften claims that suggest no existing open-source solution exists; acknowledge DELi and emphasize DELT-Hit's end-to-end integration.
- [ ] Make it explicit that the case studies cover diverse DEL architectures, including dual-display libraries.
- [ ] Add or strengthen the statistical justification for edgeR on DEL data.
- [ ] Add a concise explanation of the balance between strict barcode matching and controlled error tolerance.
- [ ] Clearly recommend validating reaction definitions on small subsets before full enumeration.
- [ ] Add a quick-start overview summarizing the workflow steps.
- [ ] Add a table listing inputs and outputs for each module.
- [ ] Revise Figure 1 so data flow and processing steps are more clearly distinguished.
- [ ] Expand QC figure captions with interpretation guidance.
- [ ] Provide a minimal YAML example inline in the manuscript.
- [ ] Clarify the recommended number of replicates and the recommended sequencing depth.
- [ ] Add a small worked chemical example.
- [ ] Add more troubleshooting guidance for low read retention.

## Reviewer 4

- [x] Expand the SMIRKS template library with more validated DEL-compatible chemistries.
- [ ] Show 2D reaction schemes alongside the corresponding SMIRKS strings.
- [ ] Add 2D structural diagrams for building blocks from the tutorial dataset.
- [ ] Add 2D structural diagrams for final products from the tutorial dataset.
- [ ] Explain how the analysis handles more nuanced selection conditions such as competition, concentration gradients, and varied stringency washes.
- [ ] Implement a normalized Z-score enrichment mode for experiments without replicates, following Faver et al.
- [x] Add an optional filtered-enumeration mode that restricts SMILES construction and descriptor calculation to barcode combinations above a user-defined threshold.
- [ ] Mention the filtered-enumeration capability clearly in the revised manuscript.

## Cross-cutting Deliverables

- [ ] Update the manuscript text to reflect all requested clarifications and new features.
- [x] Generate any new benchmark data, case-study reruns, and comparison tables needed for the revision.
- [x] Prepare revised figures, captions, troubleshooting content, and supplementary files.
- [ ] Draft a point-by-point response letter mapping each change to the editor and reviewer comments.

## Reformatting Plan Against Nature Protocols Example

Reference files:
- Current manuscript: `paper/main.tex`
- Journal example: `paper/ExampleProtocol__1776416608_11.docx`

### 1. Technical Overview must be removed as a standalone section

Current state in `main.tex`:
- `\subsection{Technical Overview}` appears near the front of the manuscript.
- It contains a module list that overlaps with material better suited to `Experimental design` or `Materials`.

Required change:
- Remove `Technical Overview` as its own top-level section.
- Move conceptual workflow text into `Introduction` or `Experimental design`.
- Move software/module inventory into `Materials`, ideally under a software-focused subsection.

### 2. Procedure hierarchy must be flattened

Current state in `main.tex`:
- `\subsection{Procedure}`
- Six `\subsubsection{Phase ...}` blocks
- Within each phase, steps are introduced with `\paragraph{Step ...}`
- Nested content is further split with `\subparagraph{1.1 ...}`, `\subparagraph{1.2 ...}`, etc.

Why this conflicts with the revision request:
- Nature Protocols wants a continuously numbered Procedure.
- The editor explicitly asked for no `1.1`, `1.2`-style numbering and only one level of subheading within the Procedure.
- The example DOCX shows phase headings followed by plain numbered instructions, not a deep heading tree.

Required change:
- Keep `Procedure` as a single main section.
- Keep phase labels only as one level of organizational heading, if desired.
- Replace `\paragraph{Step 1 | ...}` with plain numbered step text.
- Remove all `\subparagraph{1.1 ...}`, `\subparagraph{1.2 ...}`, etc.
- Convert those nested subparts into prose, lettered options, roman substeps, or short lead-in labels inside the numbered step body.

### 3. Steps must become a single continuous numbering sequence

Current state in `main.tex`:
- The Procedure is grouped into phases and then subdivided into local substeps like `1.1`, `1.2`, `2.1`, etc.

Required change:
- Renumber the entire procedure as one continuous sequence from start to finish.
- Each actionable instruction should belong to exactly one numbered step.
- Any sheet-specific setup details currently in `1.1`, `1.2`, `1.3`, etc. should either:
  - become separate consecutive numbered steps, or
  - become compact sub-items embedded inside one numbered step when they are too fine-grained to stand alone.

### 4. Optional branches must be converted into boxed options or lettered alternatives

Current state in `main.tex`:
- Optionality is often expressed inline in prose.
- The current structure does not clearly distinguish mandatory flow from optional branches.

Required change:
- For genuine branches in workflow logic, use `Option A`, `Option B`, `Option C` formatting.
- If an optional procedure interrupts the main sequence, move it into a dedicated boxed section and refer to it from the main step.
- Reserve roman numerals for substeps within a lettered option.

### 5. Explanatory text must be pulled inside steps or converted to callouts

Current state in `main.tex`:
- Many sections in `Procedure`, `Materials`, and installation text contain explanatory paragraphs sitting outside a clearly defined step.
- Example: phase introductions and several note-like paragraphs are free-standing.

Required change:
- For the Procedure, attach explanatory content directly to the relevant numbered step.
- Convert high-value warnings or constraints into standardized callouts such as `CRITICAL STEP`, `CAUTION`, or `TROUBLESHOOTING`.
- Move broader background explanation out of `Procedure` and into `Introduction` or `Experimental design`.

### 6. Timing needs step references, not just phase-level timing blocks

Current state in `main.tex`:
- Timing is embedded in phase headings such as `Phase 1 ... TIMING 30-60 min`.
- There is also a separate `Timing` subsection with a summary table.

Required change:
- Keep overall timing guidance if useful, but add explicit Procedure step references in the `Timing` section/table.
- Make sure each timing entry maps to exact step numbers rather than only to phases.
- Review whether phase-level `TIMING` labels should remain in headers or be moved into the numbered steps / timing table for consistency with journal style.

### 7. Troubleshooting must point back to Procedure steps

Current state in `main.tex`:
- `Troubleshooting` exists as its own subsection with a table and extra guidance.
- The Procedure text contains at least one generic reference to the troubleshooting section, but not systematic tagging.

Required change:
- Add the word `Troubleshooting` after each relevant Procedure step where a problem may arise.
- Add a `Step` column to the troubleshooting table so each issue maps to a concrete numbered step.
- Check that every troubleshooting entry corresponds to a real point of failure in the Procedure.

### 8. Materials should be reorganized around Nature Protocols conventions

Current state in `main.tex`:
- `Materials` contains `Equipment`, `Software dependencies`, `Input file preparation`, and `Software installation`.
- Some content reads more like setup instructions than materials inventory.

Required change:
- Keep `Materials` focused on required equipment/software/resources.
- Consider moving hands-on installation instructions into an `Experimental setup` style subsection if that is the journal-preferred location.
- Keep software lists concise and move tutorial-style explanation into Procedure steps when possible.

### 9. Data and code availability should likely be split

Current state in `main.tex`:
- There is one combined `Data and code availability` subsection.

Required change:
- Check the latest journal guidance, but likely separate this into:
  - `Data availability`
  - `Code availability`
- Ensure the code statement includes the GitHub link and later the Zenodo DOI once minted.
- The “web-based tool availability” point is not applicable here because DELT-Hit is not a web tool.

### 10. A manuscript-wide structural pass is still needed after Procedure cleanup

Current state in `main.tex`:
- Several major sections already exist: `Introduction`, `Experimental design`, `Materials`, `Procedure`, `Troubleshooting`, `Timing`, `Anticipated results`, and `Data and code availability`.
- This is already closer to Nature Protocols structure than the earlier notes suggested.

Required change:
- The biggest remaining formatting gap is not the global manuscript skeleton; it is the internal formatting of the Procedure.
- After the Procedure is reworked, do a second pass for:
  - heading depth consistency
  - callout formatting
  - timing cross-references
  - troubleshooting cross-references
  - section naming consistency with the journal template
