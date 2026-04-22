# Zero-Based vs One-Based Building-Block Codes

## Summary

There is an existing indexing mismatch between decoded count tables and enumerated
library tables.

Decoded demultiplex count codes are emitted as one-based values, while library
enumeration currently emits zero-based code values and zero-based code column
names. This can prevent direct joins or comparisons between selection counts and
the enumerated library.

## Observed Behavior

For a Cutadapt header fragment like:

```text
@read?1-B0.0?2-B1.3
```

the demultiplex postprocessor decodes the building-block suffixes as:

```python
barcodes == (1, 4)
```

The equivalent parsed building-block whitelist entries from Excel currently use:

```python
B0[0]["index"] == 0
B0[1]["index"] == 1
B0[2]["index"] == 2
```

Library enumeration then writes those values directly to the output library.

## Relevant Files

- `src/delt_hit/demultiplex/preprocess.py`
  - `write_fastq_files()` writes Cutadapt adapter records as
    `>{region.id}.{index}`.
  - The `index` here comes from Python `enumerate(region.codons)`, so adapter
    suffixes are zero-based.

- `src/delt_hit/demultiplex/postprocess.py`
  - `extract_ids()` parses selection IDs and building-block IDs from adapted read
    names.
  - Selection IDs are parsed directly.
  - Building-block IDs are parsed with `+ 1`:

    ```python
    barcodes = tuple(int(i.split('.')[-1]) + 1 for i in filter(lambda x: 'B' in x, adapters))
    ```

  - `save_counts()` writes decoded barcode tuples to `code_1`, `code_2`, ...
    count columns.

- `src/delt_hit/demultiplex/parser.py`
  - `whitelists_from_excel()` builds building-block whitelist records from Excel
    sheets.
  - It sets `df.index.name = 'index'` and then `reset_index()`, preserving
    pandas' default zero-based row index as each building block's `index`.

- `src/delt_hit/cli/library/api.py`
  - `Library.enumerate()` writes output records with:

    ```python
    record = {f'code_{i}': c['index'] for i, c in enumerate(comb)}
    ```

  - This creates zero-based column names (`code_0`, `code_1`, ...) and stores
    zero-based building-block indices.

- `documentation/cli.md`
  - The demultiplex documentation describes count output columns as `code_1`,
    `code_2`, ...

- `src/delt_hit/cli/analyse/api.py`
  - Analysis scripts assume `code_1` and `code_2` columns.

- `src/delt_hit/cli/dashboard/api.py`
  - Dashboard defaults and filtering examples use one-based-looking code values
    under `code_1`, `code_2`, ...

## Root Cause

There are two separate uses of Python `enumerate()` that encode different
semantics:

1. Demultiplex adapter names need zero-based suffixes because Cutadapt reports
   the exact adapter record names generated from Python lists.
2. User-facing count and library codes appear intended to be one-based
   identifiers (`code_1`, `code_2`, ...), matching documentation and analysis
   code.

The demultiplex path converts zero-based adapter suffixes to one-based barcode
codes during postprocessing. The library enumeration path does not perform the
same conversion and also names columns from a zero-based loop counter.

## Impact

The same physical building-block combination can be represented differently:

- Decoded counts: `code_1 = 1`, `code_2 = 4`
- Enumerated library: `code_0 = 0`, `code_1 = 3`

This is an off-by-one mismatch in values plus a column-name mismatch. Downstream
joins between counts and library records will either fail outright or silently
join against the wrong fields if columns are renamed ad hoc.

## Proposed Solution

Use one-based codes consistently for user-facing count and library outputs while
leaving Cutadapt adapter suffixes zero-based internally.

Recommended source changes:

1. In `src/delt_hit/demultiplex/parser.py`, convert building-block whitelist
   indices to one-based values when parsing Excel sheets:

   ```python
   df.index = df.index + 1
   df.index.name = 'index'
   df = df.reset_index()
   ```

2. In `src/delt_hit/cli/library/api.py`, make enumerated library code column
   names one-based:

   ```python
   record = {
       f'code_{code_position}': c['index']
       for code_position, c in enumerate(comb, start=1)
   }
   ```

3. Add regression tests covering:
   - building-block whitelist indices start at `1`;
   - `extract_ids()` converts `B*.0` to code value `1`;
   - `save_counts()` writes `code_1`, `code_2`, ... with one-based values;
   - library enumeration output uses `code_1`, `code_2`, ... and one-based
     building-block values.

## Migration Note

Existing generated config files may already contain zero-based building-block
`index` values. After applying the proposed source fix, configs should be
regenerated from Excel with `delt-hit init`, or existing configs should be
updated explicitly before re-running enumeration.

---

## Status (2026-04-22)

**Not yet started.** The issue is fully characterised; no source changes have
been made.

## Next Steps

1. **Fix `parser.py`** — shift whitelist indices to one-based on ingest
   (`df.index = df.index + 1`) so the `index` column stored in the config is
   already one-based.

2. **Fix `library/api.py`** — change `enumerate(comb)` → `enumerate(comb, start=1)`
   so output columns are `code_1`, `code_2`, ... and values match decoded counts.

3. **Add regression tests** — cover all four invariants listed in the Proposed
   Solution section above.

4. **Regenerate any existing configs** — run `delt-hit init` on affected projects
   after the fix lands; warn in the commit / PR that zero-based configs are
   incompatible with the new code.

5. **Open a PR** on a dedicated branch (e.g. `fix/one-based-codes`), separate
   from the bgzip/isal work.
