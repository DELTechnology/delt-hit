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
Fraction (or count) of allowed errors in static sections during demultiplexing.
`0.1` = 10% of section length. Higher values recover more reads but increase
false assignments.

### `library_error_correction_mode_str`
Error correction strategy for library tag matching. Format: `<metric>:<threshold>`.

| Value | Description |
|-------|-------------|
| `hamming_dist:1` | Substitutions only, max 1 error. Fast and strict. **Used in benchmarks.** |
| `levenshtein_dist:2` | Substitutions + INDELs, max 2 errors. |
| `levenshtein_dist:2,asymmetrical` | Levenshtein with asymmetric error weighting. |

### `default_error_correction_mode_str`
Same format as above; applied to non-library barcode sections (building block tags).

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
