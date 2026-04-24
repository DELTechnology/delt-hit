# DELT-Hit and DELi comparison for library enumeration

## Scope

This report compares DELT-Hit and DELi for the specific task of library enumeration: obtaining the chemical structure of each encoded compound from library definitions, building blocks, and reaction definitions. It also covers broader usability and framework-design observations because the enumeration configuration is tightly coupled to decoding and analysis in both tools.

The analysis is based on the local repositories:

- DELT-Hit: `/Users/adrianomartinelli/projects/delt-hit`
- DELi: `/Users/adrianomartinelli/projects/DELi`

## Executive summary

DELT-Hit is currently strongest as an integrated Excel-driven workflow. A single workbook can define the experiment, 
barcode structure, building blocks, reaction catalog, demultiplexing inputs, enumeration inputs, downstream properties,
and analysis. This is approachable for users who think in spreadsheets and want one artifact to configure a full run.

DELi is stronger as a library-definition and enumeration model. It separates library JSON, building-block CSV files,
reaction definitions, barcode schemas, and optional reusable data directories. 
Its enumeration model is explicit: a library owns ordered building-block sets and a reaction tree; a decoded compound
maps directly from `library_id + bb_ids` to a generated SMILES.

For robust enumeration, DELT-Hit would benefit from keeping its integrated workflow but extracting the chemistry model
into explicit, validated concepts: library definitions, building-block sets, reaction steps, named products, 
static compounds, branch/subset rules, ambiguity policies, and structured enumeration outputs.

## 1:1 comparison

| Area | DELT-Hit | DELi |
|---|---|---|
| Main scope | End-to-end DEL analysis workflow: Excel config, demultiplexing, enumeration, properties, representations, enrichment, dashboard | DEL definition, decoding, enumeration, UMI counting, analysis utilities |
| Primary input style | One Excel workbook converted to `config.yaml` | JSON library definitions plus CSV building-block files, usually under a DELi data directory |
| Library definition | Excel sheets: `reaction_graph`, `B*`, `compounds`, `reactions`, `structure`, `constant`, `selection` | Library JSON: `bb_sets`, `reactions`, optional `barcode_schema`, `scaffold`, `linker`, `doped` |
| Building blocks | Each `B*` Excel sheet has `codon`, `smiles`, `educt`, `reaction`, `product` | One CSV per BB set/cycle with `id`, optional `tag`, `smiles`, optional `subset_id` |
| Reaction definition | Flat catalog in the `reactions` sheet with `name`, `smirks`; graph edges decide how reactions connect | Reaction steps in JSON with `rxn_smarts`, `reactants`, optional `step_id`, `pick_fragment`, multiple-product flags |
| Reaction topology | NetworkX graph assembled from building-block rows and `reaction_graph` | Explicit reaction DAG/tree built from dependencies via `product_<step_id>` references |
| Intermediate non-BB steps | Supported via `reaction_graph` if connected correctly to products | Supported naturally as reaction steps whose reactants are prior products/static reagents |
| Branching chemistry | Possible but implicit; invalid multi-terminal combinations are skipped | More explicit support through reaction trees, BB subsets, pooled reactants, and branch validation |
| Static reagents/scaffolds | `compounds` sheet supplies named SMILES; referenced by graph edges | Literal SMILES in `reactants`, or reserved `scaffold`, `linker`, `truncated_linker` fields |
| Enumeration output | Parquet: `library.parquet` with `code_*` and final `smiles` | CSV/iterator: `DEL_ID`, `SMILES`, `LIB_ID`, `BB*_ID`, `BB*_SMILES` |
| Single-compound enumeration | Not exposed as a clean public primitive | First-class: `enumerate_by_bb_ids([...])` |
| Decode-to-structure mapping | Current local code uses 1-based building-block indices and `code_1`, `code_2`, etc.; mapping is numeric and spreadsheet-order dependent | Natural: decoded `bb_ids` are the same stable IDs used by building-block CSVs |
| Demultiplexing | Uses `cutadapt`; Excel `structure` and `whitelists` define regions | Built-in decode pipeline using barcode schemas, library tags, BB tags, UMI handling |
| Barcode schema | Spreadsheet `structure`, `selection`, `constant`, `B*` codons | JSON `barcode_schema` with reserved `library`, `bb1`, `bb2`, `umi`, etc. |
| Error correction | Region-level `max_error_rate`, `indels` passed to cutadapt | Barcode section error correction modes, e.g. hamming/levenshtein |
| UMI support | Not prominent in inspected enumeration/count path | Built into decoding/counting model |
| Tool/doped compounds | Not obvious in inspected path | Explicit support for tool compounds and doped compounds |
| Chemical validation | RDKit reaction application during enumeration | RDKit reaction parsing/application, product ambiguity handling, BB SMILES validation |
| Property calculation | Built-in RDKit descriptors via `library properties` | Enumeration gives structures; analysis utilities exist, but property pipeline is less central |
| Representations/fingerprints | Built-in Morgan fingerprints; current `bert` branch calls Morgan in inspected code | Not the main focus in enumeration API |
| Dashboard | Has dashboard CLI | No direct equivalent as a core feature |
| Configuration ergonomics | Easier for users who prefer Excel; one workbook contains everything | Better for version control, reuse, and validation; less friendly for spreadsheet-only users |
| Reusability | Config generated per experiment; chemistry and sequencing details are coupled | Strong separation: reusable library JSON, BB CSVs, reaction files, selection/decode configs |

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

