# Draft Response to Referees

This draft groups related reviewer comments that are addressed together. Each item is marked with a status based on the current manuscript, supplementary files, and internal notes.

## Reviewer 1

### 1. Validation on published case studies and comparison to earlier in-house analysis
**STATUS: partial**

**Reviewer comment**

The reviewer asks for stronger validation on cited case studies, noting that some earlier publications from the group used predecessor analysis code rather than DELT-Hit. They suggest re-analysis of at least one published dataset with a comparison of key results between the original method and DELT-Hit.

**Draft response**

We agree that validation on previously published datasets strengthens the protocol. In the revised materials, we now include a dedicated supplementary section, “Analysis of published datasets”, describing re-analysis of the Favalli NF2 dataset and the Pure-DEL dataset with DELT-Hit. For the Favalli dataset, the supplementary text states that the DELT-Hit workflow reproduces the published count table after alignment of identifier conventions, and an example comparison table is provided in Supplementary Table 1. We therefore now document direct re-analysis of published DEL data with the current framework. A fuller side-by-side summary of key outputs across all cited applications remains limited to the datasets for which curated comparison material has already been prepared.

### 2. Comparison with DELi
**STATUS: completed**

**Reviewer comment**

The reviewer asks for a more comprehensive and informative comparison with DELi so that readers can better understand the relative scope and complementarity of the two frameworks.

**Draft response**

We have expanded the comparison substantially in both the main text and the supplementary information. The main manuscript now includes a structured comparison table covering workflow scope, supported DEL architectures, decoding and QC, enumeration flexibility, statistical analysis, and computational performance. In addition, Supplementary Note “DELT-Hit and DELi workflow comparison” provides the rationale behind each status assignment and includes a quantitative synthetic demultiplexing benchmark comparing runtime and peak memory usage between DELT-Hit and DELi across matched benchmark datasets.

### 3. Justification for edgeR in DEL data analysis
**STATUS: completed**

**Reviewer comment**

The reviewer asks for additional clarification and support for the use of edgeR, given that DEL count data can be sparse and may not always include standard biological replicates.

**Draft response**

We have added both methodological clarification and supporting analysis. In the Procedure, the statistical-analysis section now explicitly distinguishes three enrichment modes: simple count-based ranking, replicate-aware edgeR analysis, and a normalized z-score method for experiments without replicates. The edgeR subsection now states that edgeR is used as a significance-annotation layer for replicate-based comparisons and that it requires replicates. In the supplementary information, the “Statistical methods” note explains the rationale and limitations of edgeR for DEL hit picking, and Supplementary Note “Comparison of enrichment score rankings for the Favalli CA9 selection” shows a practical comparison of count-based, edgeR, and z-score rankings for recovery of a literature-supported positive-control motif.

### 4. Language, formatting, and readability
**STATUS: partial**

**Reviewer comment**

The reviewer recommends further proofreading and formatting improvements for consistency and readability.

**Draft response**

We have revised the manuscript extensively for readability and protocol formatting. Recent updates include conversion of list-like prose in the experimental-design section into continuous narrative text, restructuring of timing information in the Procedure, standardization of callout usage such as CAUTION and CRITICAL STEP, addition of a quick-start box, and revision of figure layouts and captions. We agree that language polishing remains an important part of revision, and further copyediting may still be carried out during final manuscript preparation.

## Reviewer 2

### 1. Application runtimes and scalability across library sizes
**STATUS: partial**

**Reviewer comment**

The reviewer asks whether the application examples reflect the timing trends summarized elsewhere and recommends annotating each application example with approximate runtime to help assess practicality and scalability.

**Draft response**

We now provide clearer timing guidance at the workflow level. The manuscript includes a dedicated Timing section with step-wise estimates, and the supplementary comparison material includes synthetic runtime and memory benchmarks across matched dataset sizes. In addition, the manuscript explicitly states that demultiplexing scales primarily with sequencing depth, whereas enumeration, property calculation, and representation generation scale with library size and architectural complexity. However, the application examples in the main text are not yet all annotated with measured end-to-end runtimes for each individual published case study; where those exact run logs are still being consolidated, we prefer not to overstate precision.

### 2. Quantitative comparison with DELi
**STATUS: completed**

**Reviewer comment**

The reviewer requests a more quantitative comparison with DELi, especially using shared benchmark datasets rather than only qualitative discussion.

**Draft response**

We have addressed this by adding both a structured comparison table in the main text and a quantitative benchmark in the supplementary information. Supplementary Note “DELT-Hit and DELi workflow comparison” now includes matched synthetic demultiplexing benchmarks reporting runtime and peak RSS memory across multiple library architectures and sequencing depths. These additions complement the qualitative workflow comparison with numerical evidence on scaling behavior under shared benchmark conditions.

## Reviewer 3

### 1. Positioning relative to existing open-source tools and robustness across DEL architectures
**STATUS: completed**

**Reviewer comment**

The reviewer asks us to replace overly strong claims with more nuanced wording acknowledging DELi, and to state more explicitly that the case studies demonstrate applicability across diverse DEL architectures.

**Draft response**

We have revised the manuscript to position DELT-Hit more carefully relative to existing tools. The Introduction now states that DELi is, to our knowledge, the only other open-source package that combines barcode decoding, enrichment analysis, and library enumeration, and the manuscript emphasizes DELT-Hit’s integrated end-to-end workflow rather than implying exclusivity. We have also clarified in the Applications section that DELT-Hit has been applied across several DEL campaigns spanning single-display, dual-display, and more complex library architectures, which directly supports the intended message of robustness across different DEL formats.

### 2. Statistical methodology, barcode-error handling, and reaction-definition validation
**STATUS: completed**

**Reviewer comment**

The reviewer asks for clearer explanation of the use of edgeR, a summary statement on the balance between strict barcode matching and controlled error tolerance, and stronger guidance to validate reaction definitions on small subsets before full enumeration.

**Draft response**

We have added all three clarifications. First, the statistical-analysis section now explains the role and scope of edgeR within DELT-Hit and complements it with normalized z-score scoring for non-replicate settings. Second, the Procedure and Troubleshooting material now explicitly recommend exact matching for selection barcodes while allowing limited error tolerance in constant regions to balance assignment accuracy and read recovery. Third, the chemical-enumeration sections in both the main text and the supplementary template-reaction note now explicitly advise validation of reaction definitions on representative subsets before full-library enumeration.

### 3. Quick-start guidance, worked examples, YAML illustration, troubleshooting, and interpretation aids
**STATUS: completed**

**Reviewer comment**

The reviewer asks for several usability-oriented additions, including a quick-start overview, clearer inputs and outputs, an inline YAML example, a worked chemical example, improved troubleshooting guidance, and clearer interpretation guidance around QC and outputs.

**Draft response**

We have incorporated these additions across the revised protocol. The Procedure now begins with a Quickstart guide, the Experimental Design and Procedure sections describe the required inputs and expected outputs for the major modules, and the statistical-comparison section includes an inline YAML example for defining replicate-aware analyses. We also added a worked SMILES/SMIRKS example with a corresponding structure image and reaction explanation, expanded troubleshooting guidance including low read-retention issues, and revised several figure captions and explanatory passages to better guide interpretation of QC and analysis outputs.

### 4. Figure 1 and other presentational refinements
**STATUS: partial**

**Reviewer comment**

The reviewer suggests improving Figure 1, for example by distinguishing data flow from processing steps more clearly.

**Draft response**

We have revised the manuscript presentation and several other figures, including workflow-adjacent explanations and chemical-output visualizations. Figure 1 now describes the alignment between wet-lab workflow stages and DELT-Hit analysis modules, but we acknowledge that some visual refinements of process-versus-data-flow emphasis remain largely a figure-design question rather than a change in scientific content. We can continue refining the graphical presentation during final figure preparation.

## Reviewer 4

### 1. Expanded SMIRKS guidance and chemistry templates
**STATUS: partial**

**Reviewer comment**

The reviewer notes that advanced cheminformatics expertise remains a barrier and asks for a broader set of validated, DEL-compatible reaction templates.

**Draft response**

We agree that chemistry specification is a major usability bottleneck. The revised manuscript now points users to a supplementary collection of commonly used DEL reaction templates and emphasizes that these templates should be validated and adapted to the chemistry of the specific library. The main text and supplementary notes now also stress subset validation before full enumeration. That said, the current revision mainly strengthens guidance and documentation; it does not claim to exhaustively cover all DEL-compatible chemistries.

### 2. 2D chemical visualization alongside SMILES and SMIRKS
**STATUS: completed**

**Reviewer comment**

The reviewer asks for more intuitive chemical visualization, including 2D structural diagrams and reaction context alongside abstract formats such as SMILES and SMIRKS.

**Draft response**

We have addressed this by expanding the worked chemistry examples in the protocol and by improving the figure set used in the enumeration section. The Procedure now includes a boxed SMILES/SMIRKS explanation with a concrete SMILES example, a corresponding structure image, and an example SMIRKS explanation at the reaction level. In addition, the enumeration figures now show representative reaction templates, building-block visualizations, and representative top-hit structure renderings, which together provide more direct chemical context for the reconstructed library outputs.

### 3. Support for more nuanced selection conditions
**STATUS: partial**

**Reviewer comment**

The reviewer asks for clarification of how DELT-Hit handles more nuanced experimental designs such as competitive selections, concentration gradients, or varied stringency conditions.

**Draft response**

We have clarified that DELT-Hit separates experimental design from downstream statistical comparison through the selection sheet and analysis configuration, which allows users to define arbitrary groupings and comparisons rather than only simple protein-versus-no-protein contrasts. The Materials and Procedure sections describe metadata and grouping fields that support this flexibility. At the same time, the revised protocol still presents condition-versus-control comparisons as the main worked example, so we do not claim that every advanced selection design is documented with a full dedicated example in the current manuscript.

### 4. Statistical support for experiments without replicates
**STATUS: completed**

**Reviewer comment**

The reviewer recommends adding a normalized z-score approach for non-replicate experiments because simple count-based ranking alone is not statistically rigorous.

**Draft response**

We have addressed this recommendation by adding the normalized z-score method to the statistical-analysis section of the protocol. The Procedure now presents three enrichment modes: count-based ranking, edgeR, and normalized z-score. The z-score method is described as suitable for individual selections without replicate-based variance estimation, and the supplementary statistical note provides the governing equations and rationale.

### 5. Filtered enumeration from observed counts
**STATUS: completed**

**Reviewer comment**

The reviewer recommends a filtered-enumeration mode that restricts structure reconstruction and descriptor calculation to observed barcode combinations above a chosen threshold, rather than requiring full enumeration for very large libraries.

**Draft response**

We have implemented and documented this workflow. In the Procedure, the selection-focused enumeration step now uses \texttt{--counts\_path}, \texttt{--top\_n}, and \texttt{--library\_name} to enumerate only the top-ranked compounds from observed count data rather than the full library. This provides exactly the intended restricted-enumeration mode for practical structure generation and property calculation on selected subsets.
