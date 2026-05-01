# DELT-Hit and DELi comparison

## Scope

This report compares DELT-Hit and DELi across the end-to-end analysis steps that are most relevant to the reviewer request: demultiplexing, enumeration, hit identification, and practical workflow usability. The earlier version of this note focused mostly on enumeration; this update broadens the comparison to include the latest DELT-Hit functionality documented in `CHANGELOG.md`, `README.md`, and the benchmark runbooks.

The analysis is based on the local repositories:

- DELT-Hit: `/Users/adrianomartinelli/projects/delt-hit`
- DELi: `/Users/adrianomartinelli/projects/DELi`

## Executive summary

DELT-Hit is strongest as an integrated end-to-end workflow. A single Excel workbook can define barcode structure, building blocks, reactions, demultiplexing inputs, enumeration settings, downstream properties, and analysis comparisons. Recent additions expand that workflow with explicit dual-display support, filtered enumeration from observed counts, reviewer-friendly visualization commands, and benchmark infrastructure for shared DELT-Hit versus DELi comparisons.

DELi is strongest as an explicit library-definition and decoding/enumeration model. It separates library JSON, building-block CSV files, barcode schemas, reaction definitions, and reusable data directories. Its chemistry model is more declarative and its decoding stack already exposes concepts such as barcode-schema sections, UMI handling, tool compounds, and stable building-block identifiers.

For the reviewer response, the most balanced framing is that DELT-Hit now offers the broader manuscript-facing end-to-end workflow, while DELi remains stronger in some internals of library-definition structure and decoding model expressiveness. The updated comparison should therefore emphasize complementary strengths instead of presenting either tool as universally superior.

## Reviewer-facing comparison

### Demultiplexing

| Area | DELT-Hit | DELi |
|---|---|---|
| Core approach | Uses a `cutadapt`-based demultiplexing pipeline prepared from the DELT-Hit config and then processed into `counts.txt` tables | Uses an internal decode pipeline driven by DELi barcode schemas and decode settings |
| Input definition | Barcode structure, constants, selections, and building-block codons are authored in the Excel workbook and compiled into `config.yaml` | Barcode regions are defined explicitly in `barcode_schema` plus decode settings YAML |
| Error handling | Region-level `max_error_rate` and `indels` settings are passed through the prepared demultiplex workflow | Section-specific error correction modes such as Hamming or Levenshtein are built into the decode model |
| Flexible input preparation | DELT-Hit supports `demultiplex prepare`, `run`, `process`, `report`, and `qc`, making it easy to inspect or rerun stages independently | DELi exposes decode settings directly and keeps decode behavior close to the barcode model |
| Dual-display and strand handling | Recent DELT-Hit updates add explicit dual-display-aware parsing plus reverse-complement/complement codon handling needed for B-sheet style inputs | DELi has an explicit barcode-schema model but the benchmark docs inspected here do not foreground comparable strand-aware dual-display validation |
| Benchmarking support | Repository includes synthetic FASTQ generation, DELi/DELT-Hit input converters, split-timing runs, worker/chunk-size sweeps, and runtime/memory plotting | Benchmarks can be run on the same synthetic datasets through converted DELi inputs |
| Benchmark-specific caveat | Current DELT-Hit benchmark docs report correct counts on the inspected synthetic error datasets | The local DELi benchmark README documents known count-loss modes on synthetic `err=1` datasets and advises against treating DELi as correctness-critical for those benchmark cases until resolved upstream |

### Enumeration

| Area | DELT-Hit | DELi |
|---|---|---|
| Primary input style | One Excel workbook converted to `config.yaml` | JSON library definitions plus CSV building-block files, usually under a DELi data directory |
| Library definition | Excel sheets: `reaction_graph`, `B*`, `compounds`, `reactions`, `structure`, `constant`, `selection` | Library JSON: `bb_sets`, `reactions`, optional `barcode_schema`, `scaffold`, `linker`, `doped` |
| Reaction definition | Flat catalog in the `reactions` sheet with `name`, `smirks`; graph edges decide how reactions connect | Reaction steps in JSON with `rxn_smarts`, `reactants`, optional `step_id`, `pick_fragment`, and multiple-product controls |
| Reaction topology | NetworkX graph assembled from building-block rows and `reaction_graph` | Explicit reaction DAG/tree built from dependencies via `product_<step_id>` references |
| Intermediate non-BB steps | Supported via `reaction_graph` if connected correctly to products | Supported naturally as reaction steps whose reactants are prior products or static reagents |
| Dual-display support | Recent DELT-Hit updates add explicit dual-display support across parsing, enumeration, and visualization | DELi supports flexible library definitions, but the local materials inspected here emphasize single-library JSON/CSV workflows rather than DELT-Hit-style dual-display reviewer outputs |
| Filtered enumeration | DELT-Hit can enumerate only observed barcode combinations from `counts.txt` using `--counts_path`, `--top_n`, and `--library_name` | DELi enumeration is explicit and reusable, but the inspected local docs do not present an equivalent counts-driven filtered enumeration workflow |
| Named filtered outputs | DELT-Hit named-library support lets downstream property generation target filtered enumerations directly | DELi enumerates compounds cleanly but the inspected workflow is less centered on named downstream artifact bundles |
| Single-compound lookup | Not exposed as a clean public primitive | First-class: `enumerate_by_bb_ids([...])` |
| Decode-to-structure mapping | Uses `code_*` columns and current 1-based indexing; convenient for DELT-Hit counts joins but still spreadsheet-order dependent | Decoded `bb_ids` are the same stable IDs used by building-block CSVs |
| Enumeration output | Parquet outputs plus visualization-oriented library folders | CSV/iterator outputs centered on explicit DEL and building-block identifiers |

### Hit identification

| Area | DELT-Hit | DELi |
|---|---|---|
| Main methods | Supports counts-based enrichment and `edgeR`-based enrichment from DELT-Hit-generated count tables | Local DELi materials describe analysis/enrichment support, but the inspected repo materials emphasize decoding and library utilities more strongly than a DELT-Hit-style end-to-end hit workflow |
| Replicate-aware comparisons | DELT-Hit config supports grouping selections into analysis comparisons and recommends `edgeR` when replicates are available | DELi includes analysis utilities, but this note avoids stronger claims about replicate modeling than what is clearly documented locally |
| QC framing | DELT-Hit includes demultiplex reports, QC plots, replicate correlation framing in the manuscript, and organized per-analysis output folders | DELi analysis utilities exist, but the inspected local materials do not foreground an equivalent manuscript-ready QC workflow |
| Downstream outputs | DELT-Hit writes per-analysis counts and `edgeR` artifacts, normalized tables, hit tables, and structured output folders that feed directly into reporting | DELi provides decoded and analyzed outputs, but the practical report-generation path appears less centralized in the inspected materials |
| Interactive review | DELT-Hit includes a dashboard for interactive inspection of selection count tables | No direct equivalent is highlighted as a core DELi feature in the inspected local repo materials |

### Workflow, outputs, and usability

| Area | DELT-Hit | DELi |
|---|---|---|
| Authoring model | Excel-driven authoring is approachable for experimental users and keeps experiment design, barcode setup, chemistry, and analysis metadata in one source artifact | JSON/CSV-based definitions are more modular, reusable, and version-control-friendly |
| Visualization outputs | New `visualize enumerate` workflow exports reaction graphs, reaction schemes, per-building-block structure panels, and configured compound structure panels | The inspected DELi materials do not present an equivalent reviewer-oriented chemistry-visualization workflow as a core feature |
| Library-member export | New `visualize library` workflow exports named-library members as one-structure-per-file PNGs with filenames derived from DEL code combinations | No comparable built-in per-member export workflow is highlighted in the inspected DELi repo materials |
| Downstream chemistry artifacts | DELT-Hit integrates enumeration, property generation, Morgan/BERT-style representations, and visualization-oriented outputs in one CLI family | DELi is strong on explicit library modeling; downstream property and representation generation are less central in the inspected workflow |
| Output organization | Recent DELT-Hit changes reorganize visualization and library outputs into cleaner per-library folders that map directly onto filtered or named enumerations | DELi keeps definitions and decode settings modular, which is advantageous for reuse and code review |
| Benchmark and supporting-material workflows | DELT-Hit now includes synthetic demultiplex and enrichment benchmarking helpers, converter scripts, runbooks, and reviewer-oriented example workflows | DELi participates in the shared benchmark setup through converted inputs and reusable decode settings |

## Take-home points for the manuscript table

- DELT-Hit now has enough breadth to support a four-section comparison table rather than an enumeration-only comparison.
- The cleanest DELT-Hit strengths to emphasize are end-to-end integration, dual-display support, counts-driven filtered enumeration, reviewer-friendly visualization outputs, and benchmark infrastructure.
- The cleanest DELi strengths to preserve are explicit JSON/CSV library definitions, stable building-block identifiers, barcode-schema-driven decoding, and a more declarative reaction-step model.
- Benchmark-specific limitations documented in the local DELi benchmark README should be described carefully as observations from the shared synthetic benchmark setup, not as global claims about all DELi use cases.

## DELT-Hit enumeration model

DELT-Hit begins with an Excel workbook. `delt-hit init` parses it into `config.yaml`, and `delt-hit library enumerate` consumes that config.

The chemistry-relevant sheets are:

- `B0`, `B1`, etc.: building-block rows with `codon`, `smiles`, `educt`, `reaction`, `product`.
- `compounds`: static named compounds with `name`, `smiles`.
- `reactions`: reaction catalog with `name`, `smirks`.
- `reaction_graph`: additional reaction edges with `educt_1`, optional `educt_2`, `reaction`, `product`.

The parser builds two edge sets:

- `bb_edges`: edges derived from each selected building-block row.
- `other_edges`: edges from `reaction_graph`.

Enumeration takes the Cartesian product of the building-block whitelists, builds a subgraph for each combination, adds applicable additional-reaction subgraphs, requires exactly one terminal sink, and repeatedly fills missing product SMILES with RDKit reactions.

Important current local behavior:

- Building-block whitelist indices are 1-based in `parser.py`.
- Enumeration writes `code_1`, `code_2`, etc., matching the current count-file convention more closely than the older 0-based form.
- Output is `library.parquet` with `code_*` columns and final `smiles`.

## DELi enumeration model

DELi loads a library definition into `CombinatorialLibrary` or `DELibrary`.

The chemistry-relevant files are:

- A library JSON file under `libraries/`.
- One building-block CSV file per cycle under `building_blocks/`.
- Optional named reaction `.rxn` files under `reactions/`, or inline reaction SMARTS in the library JSON.

The library JSON contains:

- `bb_sets`: ordered building-block sets by cycle.
- `reactions`: a dictionary of named reaction steps.
- Optional static chemistry fields: `scaffold`, `linker`, `truncated_linker`.
- Optional DNA fields: `barcode_schema`, `dna_barcode_on`.
- Optional controls: `doped`.

Building-block CSVs use case-sensitive columns:

```csv
id,tag,smiles,subset_id
```

For enumeration, `id` and `smiles` are essential. `tag` supports decoding. `subset_id` supports route-specific chemistry.

Reaction steps use a dictionary form:

```json
"reactions": {
  "bb1_coupling": {
    "step_id": "1",
    "rxn_smarts": "...",
    "reactants": ["scaffold", "BB1"]
  },
  "deprotect_after_bb1": {
    "step_id": "2",
    "rxn_smarts": "...",
    "reactants": ["product_1"]
  },
  "bb2_coupling": {
    "step_id": "3",
    "rxn_smarts": "...",
    "reactants": ["product_2", "BB2"]
  }
}
```

This supports intermediate steps without new building blocks. A deprotection, reduction, hydrolysis, protecting-group removal, linker truncation, or scaffold activation can be represented as a normal step whose reactants are prior products and/or static reagents.

## DELi product ambiguity handling

DELi has two distinct ambiguity situations.

### Multiple product fragments

This is when the reaction SMARTS has more than one product template, i.e. multiple fragments on the product side. DELi handles this at reaction initialization.

- If `pick_fragment` is not provided and RDKit reports more than one product template, DELi emits a `ReactionWarning`.
- It defaults to `pick_fragment = 0`, meaning the first product fragment is selected.
- If `pick_fragment` is provided, it is 1-indexed in configuration and converted to 0-index internally.
- If `pick_fragment` is outside the product-template range, DELi raises `ReactionParsingError`.

The useful concept here is that fragment selection is a declared reaction-step policy, not hidden in the enumeration loop.

### Multiple reaction outcomes

This is when RDKit returns multiple possible product tuples for the same reactant set, usually because the reaction SMARTS is broad or a reactant has multiple matching sites.

DELi handles this in `Reaction.react()`:

- If RDKit returns no products, it raises `ReactionError`.
- If RDKit returns multiple product outcomes and `fail_on_multiple_products` is true, it raises `ReactionError`.
- If RDKit returns multiple product outcomes and `ignore_multiple_products` is false, it emits `ReactionWarning`.
- In non-failing mode, it selects the first product outcome: `products[0][pick_fragment]`.
- The product property cache is updated before returning.

The reaction-step JSON can specify:

```json
{
  "rxn_smarts": "...",
  "reactants": ["product_1"],
  "pick_fragment": 1,
  "ignore_multiple_products": false,
  "fail_on_multiple_products": true
}
```

The policy is imperfect chemically because "first RDKit product" is not always a meaningful choice, but it is explicit and auditable. For production enumeration, `fail_on_multiple_products: true` is safer when reaction SMARTS are expected to be deterministic.

### Multiple alternative reactions

DELi also lets `rxn_smarts` be a list. In that case a `ReactionStep` tries each reaction in priority order and uses the first one that successfully produces a product. This is useful when one conceptual step may require multiple SMARTS patterns, such as primary amine first, secondary amine second.

This is different from multiple products. It is a priority list of alternative reaction definitions for one step.

### Enumeration failure handling

During library enumeration, DELi validates that each selected building block has a valid SMILES and can be converted to an RDKit Mol. If the reaction thread fails:

- `fail_on_error=True`: raise `EnumerationRunError`.
- `dropped_failed=True`: skip failed compounds.
- Default: yield a compound without an enumerated Mol/SMILES.

When exporting to file, failed compounds are written with `SMILES = NA`.

## Good design choices in DELT-Hit

DELT-Hit has several good ingredients worth preserving.

- A single Excel template is accessible to experimental users.
- The same source workbook drives demultiplexing and chemistry, reducing file scattering.
- `delt-hit init` creates an inspectable YAML config, which is better than using the workbook directly at runtime.
- Reaction graph PNGs are emitted, which helps users debug topology.
- Enumeration, property calculation, and representations are integrated into one CLI group.
- The `compounds`/`reactions` catalog is simple and spreadsheet-friendly.
- The current local move toward 1-based code indices is aligned with decoded count tables and removes an important source of confusion.

## Suboptimal design choices in DELT-Hit

The current enumeration design works for simple cases but becomes opaque as chemistry grows.

### Chemistry topology is split across row-level BB fields and `reaction_graph`

Each building-block row contains `educt`, `reaction`, and `product`, while additional reactions are in a separate `reaction_graph` sheet. This makes it hard to see the full reaction scheme in one place.

For example, the user must infer that a building block row defines three edges:

- building-block set to reaction
- educt to reaction
- reaction to product

Then the `reaction_graph` sheet defines other edges in the same graph language. That is powerful, but not self-documenting.

### Additional-reaction inclusion is implicit

The enumeration code adds connected components from the additional-reaction graph when a node in the selected BB subgraph is a sink in that additional subgraph. This is clever, but difficult to explain and debug. A chemist expects to say "after product_1, run deprotection", not "define a subgraph whose sink matches a selected node".

### Product ambiguity is not configurable

DELT-Hit `perform_reaction()` collects all RDKit products into a set of canonical SMILES and `complete_reaction_graph()` asserts that exactly one product was found. If there are zero products, it retries reactants in reversed order. If there are multiple products, enumeration fails unless `errors='ignore'` skips the compound.

There is no per-reaction policy equivalent to:

- choose product fragment N
- warn and choose first
- fail on multiple products
- ignore warnings
- try alternative reaction SMARTS in priority order

### Reactant order is handled by retrying reversed inputs

For two-reactant reactions, DELT-Hit tries one order and then the reversed order if no products are found. This can rescue some cases, but it does not scale cleanly to three or more reactants and can mask mismatches between graph ordering and reaction SMARTS.

### Building-block identity is numeric and order dependent

The enumeration output uses `code_*` based on row indices. That is convenient for matching decoded numeric barcodes, but less robust than stable building-block IDs. If rows are inserted or reordered, the chemistry identity can drift unless indices/codons are carefully controlled.

### Config validation is shallow

Current validation mainly asserts that all library reactions are present in the reaction catalog. It does not deeply validate:

- all graph nodes resolve to products, compounds, reactions, or building-block sets
- every reaction has the expected number of reactants
- every reaction graph is acyclic
- every combination has exactly one terminal product
- every building-block SMILES parses
- every reaction SMARTS parses before enumeration
- decoded count columns will join to enumerated code columns

### Enumeration concerns are concentrated in one large method

`Library.enumerate()` performs path resolution, config reading, graph construction, plotting, combinatorial expansion, graph selection, RDKit chemistry, debug image generation, and output writing. This makes it hard to test or reuse single-compound enumeration.

### Error handling can be surprising

Inside `complete_reaction_graph()`, some errors call `exit()` when `errors == 'raise'`, rather than raising a typed exception back to the CLI. This makes failures harder to catch in tests and harder to summarize for users.

### Some CLI/API naming and implementation details are opaque

Examples:

- `debug` is a string with values like `'False'`, `'all'`, `'valid'`, `'invalid'` rather than a typed option.
- `building_block_ids` filters building-block sheet names, not building-block row IDs.
- `represent(method='bert')` currently calls `run_morgan()` in the inspected code, so the implementation does not match the option name.

## Roadmap: cleaner DELT-Hit enumeration

The goal should not be to copy DELi. DELT-Hit can keep its Excel-first, end-to-end workflow while introducing a clearer internal chemistry model.

### Phase 1: Define explicit internal domain objects

Introduce typed objects or Pydantic/dataclass schemas for:

- `LibraryDefinition`
- `BuildingBlockSet`
- `BuildingBlock`
- `StaticCompound`
- `ReactionDefinition`
- `ReactionStep`
- `ReactionGraph`
- `EnumerationRecord`
- `EnumerationFailure`

Keep Excel and YAML as input/output formats, but parse them into these objects before enumeration. The enumerator should not work directly against raw dictionaries.

### Phase 2: Split barcode config from chemistry config

Keep a single Excel workbook if desired, but separate the generated YAML into clearer top-level sections:

```yaml
experiment: ...
barcodes:
  structure: ...
  constants: ...
  selections: ...
  building_block_tags: ...
chemistry:
  building_block_sets: ...
  static_compounds: ...
  reactions: ...
  reaction_steps: ...
analysis: ...
```

This makes it obvious which configuration is needed for demultiplexing and which is needed for enumeration.

### Phase 3: Replace implicit graph construction with explicit reaction steps

Add a new recommended chemistry-definition format:

```yaml
chemistry:
  building_block_sets:
    - id: B0
      id_column: id
      tag_column: codon
      smiles_column: smiles
    - id: B1
  static_compounds:
    scaffold_1:
      smiles: O=C(N)CCC(N=[N+]=[N-])C(O)=O
  reactions:
    ABF:
      smirks: "..."
      ambiguity:
        multiple_products: fail
        pick_fragment: 1
  steps:
    - id: bb1_coupling
      reaction: ABF
      reactants: [scaffold_1, B0]
      product: product_1
    - id: deprotect_1
      reaction: SR
      reactants: [product_1]
      product: product_1B
    - id: bb2_coupling
      reaction: ABF
      reactants: [product_1B, B1]
      product: product_2
```

The existing `reaction_graph`/`B*` sheet style can remain as a compatibility layer that compiles into this step model.

### Phase 4: Make ambiguity and product selection explicit

Add per-reaction or per-step policies:

```yaml
ambiguity:
  no_product: fail | skip
  multiple_products: fail | warn_first | first | all
  multiple_fragments: pick
  pick_fragment: 1
  canonicalize: true
```

Recommended defaults for production:

- fail on multiple products
- fail on no product
- require `pick_fragment` if multiple fragments are present

Recommended defaults for exploratory mode:

- warn and choose first
- record all warnings in an enumeration diagnostics table

### Phase 5: Add explicit reactant binding

Avoid relying on graph predecessor sort order or reversed-reactant retry. Each step should preserve reactant order as declared by the user because RDKit reaction SMARTS are order-sensitive.

For advanced cases, support optional named roles:

```yaml
reactants:
  - role: acid
    source: product_1
  - role: amine
    source: B1
```

Even if roles are not used by RDKit, they improve readability and validation.

### Phase 6: Support branch/subset rules clearly

For route-specific chemistry, introduce a subset column in building-block sheets, or a separate sheet that maps rows to subsets:

```yaml
building_block_sets:
  B1:
    subset_column: chemistry_class

steps:
  - id: amide_route
    applies_to:
      B1: amines
    reaction: amide_coupling
    reactants: [product_1, B1]
    product: product_2
  - id: ester_route
    applies_to:
      B1: alcohols
    reaction: esterification
    reactants: [product_1, B1]
    product: product_2
```

The core idea to borrow from DELi is not its exact `BB1:::subsetA` syntax, but the principle that branch applicability should be declarative, validated, and visible.

### Phase 7: Add single-compound enumeration

Expose an API and CLI like:

```bash
delt-hit library enumerate-one --config_path config.yaml --B0 12 --B1 48
```

or with stable IDs:

```bash
delt-hit library enumerate-one --config_path config.yaml --bb B0=A012 --bb B1=B048
```

This would make debugging much easier and directly supports decoded-hit structure lookup.

### Phase 8: Improve output schema

The output should include both numeric codes and stable building-block identifiers:

- `library_id`
- `compound_id`
- `smiles`
- `enumeration_status`
- `enumeration_error`
- `code_1`, `code_2`, etc.
- `B0_id`, `B0_smiles`
- `B1_id`, `B1_smiles`
- `final_product_node`
- optional `reaction_path_id`

This makes downstream analysis robust even if a building-block table is reordered later.

### Phase 9: Add validation and dry-run commands

Add:

```bash
delt-hit library validate --config_path config.yaml
delt-hit library graph --config_path config.yaml
delt-hit library enumerate-one --config_path config.yaml ...
```

Validation should report:

- missing columns
- invalid SMILES
- invalid SMARTS/SMIRKS
- unresolved reactants/products
- duplicate IDs/codons
- cycles in the reaction graph
- reaction arity mismatches
- ambiguous product policies
- count/enumeration join compatibility

The graph command should render the intended chemistry graph before enumeration, not only the assembled graph after the code has inferred extra edges.

## Roadmap: broader DELT-Hit framework improvements

### 1. Make configuration layers explicit

Current configuration mixes experiment metadata, barcode parsing, chemistry, selection metadata, and analysis. A clearer configuration hierarchy would reduce accidental coupling.

Recommended layers:

- `experiment`: paths, run name, compute settings
- `barcode_schema`: read regions, constants, selections, BB tag regions
- `chemistry`: library structure and enumeration
- `decoding`: cutadapt/error-correction settings
- `analysis`: groups and comparisons
- `outputs`: generated paths and artifact metadata

### 2. Treat Excel as an authoring format, not the canonical schema

Excel should remain supported, but the canonical schema should be typed YAML/JSON with a documented spec. The Excel parser should compile to the schema and produce detailed validation errors.

This preserves usability while making configs reviewable, testable, and versionable.

### 3. Add schema versioning and migrations

Add:

```yaml
schema_version: 1
```

Then provide migration functions as the config evolves. Without versioning, backward compatibility becomes guesswork.

### 4. Replace assertions with typed validation errors

Many parser checks use `assert`. Assertions can be disabled and often produce terse messages. Use explicit exceptions such as:

- `ConfigValidationError`
- `ChemistryValidationError`
- `DemultiplexValidationError`
- `EnumerationError`

Each error should include the sheet/section, row, field, observed value, and expected value.

### 5. Improve testing around realistic chemistry

Add unit tests for:

- one-step enumeration
- intermediate deprotection step
- static scaffold reaction
- two-reactant order validation
- multi-product ambiguity
- multi-fragment reaction with `pick_fragment`
- branch/subset route selection
- invalid graph with multiple sinks
- count table join to enumeration output

### 6. Normalize terminology

Current naming mixes `library`, `structure`, `catalog`, `whitelists`, `B*`, `educt`, `product`, `compound`, and `code_*`. Define a glossary and align file/schema fields:

- building block set
- building block
- tag/codon
- static compound
- reaction definition
- reaction step
- intermediate product
- terminal product
- decoded code
- compound ID

This will reduce user confusion more than almost any implementation change.

### 7. Make debugging artifacts first-class

Instead of writing only PNGs, write structured diagnostics:

- `reaction_graph.nodes.tsv`
- `reaction_graph.edges.tsv`
- `enumeration_failures.parquet`
- `reaction_warnings.parquet`
- `invalid_combinations.parquet`

PNG graphs are helpful for humans, but structured files are needed for large libraries.

### 8. Separate library enumeration from feature generation

Keep the CLI grouping, but internally separate:

- enumerate structures
- calculate properties
- calculate representations
- render QC plots

This will make each stage easier to cache, test, and rerun.

### 9. Make output paths predictable and configurable

Currently output location is tied to `experiment.save_dir / experiment.name`. That is reasonable, but a richer output manifest would help:

```yaml
outputs:
  library: ...
  properties: ...
  representations: ...
  demultiplex: ...
  analysis: ...
```

This helps dashboards and downstream tools avoid recomputing path conventions.

### 10. Improve CLI option types and names

Use typed booleans/enums instead of strings for options such as `debug`. Rename options where intent is ambiguous:

- `building_block_ids` -> `building_block_set_ids` if filtering B0/B1 sheets
- add `--bb-id` or `--code` for row-level filtering
- `errors` -> `--on-error fail|skip|record`

### 11. Add stable compound IDs

Generate a stable compound identifier such as:

```text
<library_name>-B0:<bb_id>-B1:<bb_id>
```

If only numeric codes exist, generate:

```text
<library_name>-code_1:<n>-code_2:<n>
```

This should be included in enumeration, counts joins, enrichment, and dashboard views.

## Recommended target design

The optimal design for DELT-Hit is a hybrid:

- Keep the Excel template for usability.
- Compile Excel into a versioned, typed chemistry/barcode schema.
- Define chemistry as explicit reaction steps with ordered reactants and named products.
- Keep graph visualization, but make graph construction transparent.
- Add stable building-block IDs and explicit compound IDs.
- Add per-reaction ambiguity policy.
- Provide single-compound enumeration for hit lookup.
- Record failures and warnings as data, not only logs.

This keeps what is good in DELT-Hit: integrated workflow, spreadsheet accessibility, full analysis pipeline, dashboard potential. It adds what is good in DELi at the conceptual level: reusable library definitions, explicit reaction steps, stable BB IDs, branch-aware enumeration, and clear `library_id + bb_ids -> structure` mapping.

## Immediate action items

1. Document the current Excel chemistry contract and generated YAML schema.
2. Add `delt-hit library validate`.
3. Extract RDKit reaction execution into a small tested module with ambiguity policies.
4. Add stable building-block IDs to `B*` sheets and enumeration outputs.
5. Add `enumerate-one` for debugging and decoded-hit structure lookup.
6. Add structured failure outputs for enumeration.
7. Create a new optional `reaction_steps` sheet/schema while keeping current `reaction_graph` support as a compatibility layer.
8. Write migration examples showing how current `reaction_graph` libraries map to explicit `reaction_steps`.
