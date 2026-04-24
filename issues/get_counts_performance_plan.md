# Plan: Cython C-level Chunk Parsing for get_counts

Context document: `issues/get_counts_performance.md`

## Goal

Replace `_count_bytes_chunk` (Python bytes ops, one allocation per token per line) with a Cython
C function `parse_chunk` that uses `libc.string.memchr` and manual digit parsing. Keep all other
infrastructure (isal decompression, `_iter_byte_chunks`, `multiprocessing.Pool`, merge loop)
unchanged. Also unify the serial path to use isal + the same C parser.

---

## Open Tasks

### 1. Add Cython to the build system

- [ ] Edit `pyproject.toml`: add `"Cython >= 3.0"` to `[build-system] requires`
- [ ] Create `setup.py` with the extension definition (see skeleton below)

```python
# setup.py
from setuptools import setup, Extension
from Cython.Build import cythonize

setup(
    ext_modules=cythonize(
        [Extension(
            "delt_hit.demultiplex._parse",
            ["src/delt_hit/demultiplex/_parse.pyx"],
            extra_compile_args=["-O3"],
        )],
        compiler_directives={"language_level": "3"},
    )
)
```

---

### 2. Write `src/delt_hit/demultiplex/_parse.pyx`

Drop-in replacement for `_count_bytes_chunk`. Same input (`bytes` chunk), same return type
(`dict[tuple, dict[tuple, int]]`).

```cython
# cython: language_level=3, boundscheck=False, wraparound=False
from libc.string cimport memchr

def parse_chunk(bytes chunk):
    cdef:
        const unsigned char* p = chunk
        const unsigned char* end = p + len(chunk)
        const unsigned char* nl
        const unsigned char* line_end
        const unsigned char* q
        const unsigned char* dot
        unsigned char first_char
        int val
        list sel_ids, bar_ids
        tuple sel_key, bar_key
        dict counts = {}
        dict bc_dict

    while p < end:
        nl = <const unsigned char*>memchr(p, b'\n'[0], end - p)
        line_end = nl if nl is not NULL else end

        if line_end > p:
            q = <const unsigned char*>memchr(p, b'?'[0], line_end - p)
            if q is not NULL:
                p = q + 1
                sel_ids = []
                bar_ids = []

                while p < line_end:
                    q = <const unsigned char*>memchr(p, b'?'[0], line_end - p)
                    if q is NULL:
                        q = line_end

                    if q > p:
                        first_char = p[0]
                        dot = q - 1
                        while dot > p and dot[0] != b'.'[0]:
                            dot -= 1
                        if dot[0] == b'.'[0]:
                            dot += 1
                            val = 0
                            while dot < q and b'0'[0] <= dot[0] <= b'9'[0]:
                                val = val * 10 + (dot[0] - b'0'[0])
                                dot += 1
                            if first_char == b'S'[0]:
                                sel_ids.append(val)
                            elif first_char == b'B'[0]:
                                bar_ids.append(val + 1)  # +1 matches Python impl

                    p = q + 1

                sel_key = tuple(sel_ids)
                bar_key = tuple(bar_ids)
                bc_dict = counts.get(sel_key)
                if bc_dict is None:
                    counts[sel_key] = {bar_key: 1}
                else:
                    bc_dict[bar_key] = bc_dict.get(bar_key, 0) + 1

        p = (nl + 1) if nl is not NULL else end

    return counts
```

---

### 3. Update `postprocess.py` — two minimal edits

**Edit A** — swap `_count_bytes_chunk` import at module level (try C, fall back to Python):

```python
try:
    from delt_hit.demultiplex._parse import parse_chunk as _count_bytes_chunk
except ImportError:
    pass  # keep existing Python _count_bytes_chunk as fallback
```

Place this after the existing `_count_bytes_chunk` definition so the Python version is the
fallback if the extension is not built.

**Edit B** — unify the serial path (`num_workers == 1`) to also use isal + `_count_bytes_chunk`
instead of stdlib gzip + `extract_ids`:

```python
if num_workers == 1:
    from isal import igzip
    counts: dict = {}
    with igzip.open(input_path, 'rb') as f:
        for chunk in tqdm(_iter_byte_chunks(input_path, chunk_size_bytes), ncols=100):
            for sel, bc_counts in _count_bytes_chunk(chunk).items():
                if sel not in counts:
                    counts[sel] = bc_counts
                else:
                    existing = counts[sel]
                    for bc, n in bc_counts.items():
                        existing[bc] = existing.get(bc, 0) + n
    return counts
```

This means the serial path gets both the isal decompression speedup and the C-level parsing.

---

### 4. Build and smoke-test

```bash
uv run python setup.py build_ext --inplace

uv run python -c "
from delt_hit.demultiplex._parse import parse_chunk
chunk = b'r1?S0.GAT.2?B1.GCC.5\nr2?S0.GAT.2?B1.GCC.5\n'
print(parse_chunk(chunk))
# Expected: {(2,): {(6,): 2}}
"
```

---

### 5. Correctness check against Python implementation

```bash
uv run python - <<'EOF'
from delt_hit.demultiplex._parse import parse_chunk
from delt_hit.demultiplex.postprocess import _count_bytes_chunk, _iter_byte_chunks
from pathlib import Path

path = Path("path/to/reads_with_adapters.gz")
for chunk in _iter_byte_chunks(path, 5_000_000):
    r_c  = parse_chunk(chunk)
    r_py = _count_bytes_chunk(chunk)  # Python fallback must still exist
    assert r_c == r_py, f"mismatch"
    break
print("First chunk matches")
EOF
```

---

### 6. Extend the benchmark to include the Cython path

Update `tests/benchmark_get_counts.py` to add a third timing row:

| Mode | Time | Speedup |
|---|---|---|
| Serial, stdlib gzip, Python extract_ids | baseline | 1× |
| Parallel N workers, isal, Python _count_bytes_chunk | measured | ~4.7× |
| Parallel N workers, isal, Cython parse_chunk | TBD | expected ~15–20× |

Run with:

```bash
uv run python tests/benchmark_get_counts.py
```

---

### 7. (Optional) Verify with a real production file

Run `get_counts` end-to-end on an actual `reads_with_adapters.gz` from a previous experiment and
compare counts to the saved `*_counts.txt` output. Confirms the +1 barcode offset, multi-selection
handling, and output format are correct at production scale.
