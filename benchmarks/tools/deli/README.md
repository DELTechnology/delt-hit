# DELi Benchmark Configuration

DELi decode settings used across benchmark datasets are defined in per-dataset
`decode_settings_v02.yaml` files. This document explains each setting and the
available options.

## Decode Settings

### `ignore_tool_compounds`
Skips "tool compounds" — single compounds added outside the main DEL library
structure (e.g., positive controls). Set to `true` for benchmarks to avoid
wrapping them as single-compound libraries.

### `demultiplexer_algorithm`
Controls how reads are matched to libraries.

| Value | Description |
|-------|-------------|
| `regex` | Fuzzy regex on static sections with wildcard gaps. Fast; handles ±5 bp INDELs. **Used in benchmarks.** |
| `cutadapt` | Cutadapt adapter logic (Front/Back/Anywhere). Good for terminal adapters; limited to 2 static sections. |
| `full` | Full semi-global alignment of the entire read. Most thorough but slowest. |

### `demultiplexer_mode`
Controls which barcode region is used for library assignment.

| Value | Description |
|-------|-------------|
| `library` | Uses only the library tag. One query per library. **Used in benchmarks.** |
| `single` | Picks one static section spanning the most libraries; resolves library via barcode calling. |
| `flanking` | Uses two static sections flanking the library tag. Most accurate. |

### `realign`
When `true`, performs a second semi-global pairwise alignment after initial
demultiplexing to recover reads with multiple INDELs. Increases runtime.
Set to `false` for benchmarks.

### `library_error_tolerance`
Allowed errors in static sections during demultiplexing. Accepts a **float** or
an **integer** with different meanings:

| Value type | Interpretation |
|------------|----------------|
| **float** `0 < v < 1` | Fraction of the static section length, e.g. `0.1` = 10%. Actual error budget = ⌊rate × length⌋, so the effective count varies with section length. |
| **int** `≥ 1` | Absolute number of allowed errors, unambiguous regardless of section length. |

> **Benchmarks use integers** — `0` for `err=0` datasets, `1` for `err=1`
> datasets. Avoid floats unless you specifically want length-proportional
> tolerance.

### `library_error_correction_mode_str`
Error correction strategy for library tag matching. Format: `<metric>:<threshold>`.
The threshold is always an **absolute integer** count of allowed errors.

| Value | Description |
|-------|-------------|
| `hamming_dist:1` | Substitutions only, max 1 error. Fast and strict. **Used in benchmarks.** |
| `levenshtein_dist:2` | Substitutions + INDELs, max 2 errors. |
| `levenshtein_dist:2,asymmetrical` | Levenshtein with asymmetric error weighting. |

### `default_error_correction_mode_str`
Same format as `library_error_correction_mode_str`; applied to building block
tag sections.

### `min_library_overlap`
Minimum base overlap between read and library sequence for a valid match.
Primarily relevant for the `cutadapt` algorithm. Higher = more specific.

### `revcomp`
When `true`, also searches the reverse complement of each read. Set to `false`
for synthetic benchmarks where strand orientation is fixed.

### `library_wiggle` / `wiggle`
Allow positional flexibility in library tag / building block tag matching.
Both disabled (`false`) in benchmarks for stricter, reproducible matching.

### `decode_matching_approach`
Strategy for resolving reads that match multiple libraries.

| Value | Description |
|-------|-------------|
| `first_perfect` | Assign to the first perfect match found. **Used in benchmarks.** |

### `min_read_length` / `max_read_length`
Reads outside `[50, 128]` bp are discarded before decoding.

---

## Dataset-Specific Notes

### `synthetic_4cycle_100m_err=1`
Uses the settings above with `hamming_dist:1` error correction. The `err=1`
suffix indicates the synthetic FASTQ was generated with a per-base error rate
of 1%, exercising error correction under realistic sequencing noise.

---

## Known Issues

### DELi fails to recover all reads when errors appear outside building block regions

**Affected datasets:** any `err=1` dataset where errors are injected into all
barcode regions (library tag, BB tags, closing tag).

**Symptom:** `counts match: False` — DELi decodes fewer compounds than expected
even with `library_error_tolerance: 1` and `default_error_correction_mode_str:
hamming_dist:1` correctly configured.

**Root cause (investigated 2026-04-26):** DELi uses a fuzzy regex
(`{e<=N}`) to locate the library tag and derive the positions of all downstream
BB sections by fixed offset arithmetic. When the library tag itself carries a
substitution error, the fuzzy match can anchor at a shifted position, causing
the inferred BB1 span to be off by one base. This leads to a failed barcode
lookup for BB1 even though the BB tag itself is within hamming distance 1 of a
known barcode. All 4 436 failures observed in the `err=1` run were at `bb1`
(the first section after the library tag); bb2–bb4 never fail because the read
is dropped as soon as bb1 fails.

**Verification:** Re-generating the dataset with errors injected only into BB
tags (`--errors-in-bb-only true`) leaves the library tag perfect, and DELi
still does not decode 100 % of reads (99 610 / 100 000 in the `_bbonly` run).
DELT-Hit returns correct counts for both datasets.

**Conclusion:** DELi has at least two sources of count loss:
1. Shifted section spans when the library tag is fuzzy-matched (errors in all
   regions).
2. A secondary loss of ~390 reads even when the library tag is error-free,
   indicating an additional bug in the BB span or barcode-calling logic.

DELi is not used for correctness-critical benchmarking until these issues are
resolved upstream.
