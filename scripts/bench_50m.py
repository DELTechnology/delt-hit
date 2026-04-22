"""Focused benchmark: batched bgz vs batched plain at 50M reads."""
import random, subprocess, tempfile, time, multiprocessing, os
from collections import defaultdict
from pathlib import Path
from multiprocessing import cpu_count
from isal import igzip

NUM_READS = 50_000_000
BATCH_SIZE = 150_000
WORKERS = [1, cpu_count()]


def _make_line(rng, num_selections=2, num_barcodes=3, sel_range=4, bc_range=100):
    sel_ids = [rng.randint(0, sel_range - 1) for _ in range(num_selections)]
    bc_ids  = [rng.randint(0, bc_range - 1)  for _ in range(num_barcodes)]
    sel_tokens = [f"X.S{i}.{s}" for i, s in enumerate(sel_ids, 1)]
    bc_tokens  = [f"X.B{i}.{b}" for i, b in enumerate(bc_ids, 1)]
    all_tokens = sel_tokens + bc_tokens
    rng.shuffle(all_tokens)
    return "read_" + str(rng.randint(0, 10**9)) + "?" + "?".join(all_tokens) + "\n"


def _count_batch_worker(lines):
    from delt_hit.demultiplex.postprocess import extract_ids
    counts = defaultdict(dict)
    for line in lines:
        sel, bc = extract_ids(line)
        inner = counts[sel]
        inner[bc] = inner.get(bc, 0) + 1
    return dict(counts)


def _batches(f, batch_size=BATCH_SIZE):
    batch = []
    for line in f:
        batch.append(line)
        if len(batch) == batch_size:
            yield batch
            batch = []
    if batch:
        yield batch


def bench_batched_bgz(path, nw):
    from delt_hit.demultiplex.postprocess import _merge
    counts = defaultdict(dict)
    t0 = time.perf_counter()
    with igzip.open(path, "rt") as f, multiprocessing.Pool(nw) as pool:
        for partial in pool.imap(_count_batch_worker, _batches(f)):
            _merge(counts, partial)
    return time.perf_counter() - t0


def bench_batched_plain(path, nw):
    from delt_hit.demultiplex.postprocess import _merge
    counts = defaultdict(dict)
    t0 = time.perf_counter()
    with open(path) as f, multiprocessing.Pool(nw) as pool:
        for partial in pool.imap(_count_batch_worker, _batches(f)):
            _merge(counts, partial)
    return time.perf_counter() - t0


if __name__ == "__main__":
    rng = random.Random(42)
    tmp_txt = Path(tempfile.mktemp(suffix=".tmp.txt"))
    bgz_path = Path(tempfile.mktemp(suffix=".bgz"))
    txt_path = Path(tempfile.mktemp(suffix=".txt"))

    print(f"Writing {NUM_READS:,} lines...", flush=True)
    with open(tmp_txt, "w") as fh:
        for _ in range(NUM_READS):
            fh.write(_make_line(rng))

    print("Compressing with bgzip...", flush=True)
    subprocess.run(["bgzip", "-f", "-@", str(cpu_count()), str(tmp_txt)], check=True)
    Path(str(tmp_txt) + ".gz").rename(bgz_path)
    tmp_txt.unlink(missing_ok=True)

    print("Decompressing to plain text...", flush=True)
    subprocess.run(["bgzip", "-d", "-k", "-@", str(cpu_count()), "-c", str(bgz_path)],
                   stdout=open(txt_path, "wb"), check=True)

    bgz_size = bgz_path.stat().st_size / 1e6
    txt_size = txt_path.stat().st_size / 1e6
    print(f"\nFile sizes: compressed={bgz_size:.1f} MB  uncompressed={txt_size:.1f} MB  ratio={txt_size/bgz_size:.2f}x\n")

    try:
        print("=" * 60)
        for nw in WORKERS:
            t = bench_batched_bgz(bgz_path, nw)
            print(f"batched bgz   w={nw:<3}: {t:.2f} s  ({NUM_READS/t/1e6:.2f} M reads/s)", flush=True)

        for nw in WORKERS:
            t = bench_batched_plain(txt_path, nw)
            print(f"batched plain w={nw:<3}: {t:.2f} s  ({NUM_READS/t/1e6:.2f} M reads/s)", flush=True)
        print("=" * 60)
    finally:
        bgz_path.unlink(missing_ok=True)
        txt_path.unlink(missing_ok=True)
