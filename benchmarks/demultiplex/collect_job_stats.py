#!/usr/bin/env python3
"""Collect Slurm job stats for benchmark runs and write job-stats.json next to timing.json.

For each tool/dataset directory that contains a timing.json, finds the most recent
COMPLETED Slurm job whose name matches '{tool}-{dataset}' and writes job-stats.json
with peak memory, CPU usage, wall time, and other useful metrics.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
TOOLS_ROOT = PROJECT_ROOT / "benchmarks" / "demultiplex" / "tools"

SACCT_FIELDS = [
    "JobID",
    "JobName",
    "State",
    "Elapsed",
    "ElapsedRaw",
    "CPUTimeRAW",
    "NCPUS",
    "NNodes",
    "MaxRSS",
    "MaxVMSize",
    "AveRSS",
    "Start",
    "End",
    "ExitCode",
    "NodeList",
]


def run_sacct(user: str, start_date: str = "2020-01-01") -> list[dict]:
    cmd = [
        "sacct",
        "-u", user,
        f"--format={','.join(SACCT_FIELDS)}",
        "-S", start_date,
        "-P",
        "--noheader",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    rows = []
    for line in result.stdout.splitlines():
        parts = line.split("|")
        if len(parts) != len(SACCT_FIELDS):
            continue
        rows.append(dict(zip(SACCT_FIELDS, parts)))
    return rows


def parse_max_rss_kb(value: str) -> int | None:
    """Convert sacct memory string (e.g. '14284692K') to bytes."""
    if not value:
        return None
    value = value.strip()
    multipliers = {"K": 1024, "M": 1024**2, "G": 1024**3, "T": 1024**4}
    if value and value[-1] in multipliers:
        return int(float(value[:-1]) * multipliers[value[-1]])
    if value.isdigit():
        return int(value)
    return None


def elapsed_to_seconds(elapsed: str) -> int | None:
    """Convert HH:MM:SS or D-HH:MM:SS to seconds."""
    if not elapsed:
        return None
    try:
        days = 0
        if "-" in elapsed:
            day_part, elapsed = elapsed.split("-", 1)
            days = int(day_part)
        parts = elapsed.split(":")
        h, m, s = int(parts[0]), int(parts[1]), int(parts[2])
        return days * 86400 + h * 3600 + m * 60 + s
    except Exception:
        return None


def build_job_index(rows: list[dict]) -> dict[str, list[dict]]:
    """Group batch step rows by JobName, keeping only COMPLETED parent jobs."""
    # Collect parent job IDs that completed successfully
    completed_ids: set[str] = set()
    parent_rows: dict[str, dict] = {}
    for row in rows:
        job_id = row["JobID"]
        if "." in job_id:
            continue  # skip steps here
        if row["State"] == "COMPLETED":
            completed_ids.add(job_id)
            parent_rows[job_id] = row

    # Collect batch step rows for completed jobs (they hold MaxRSS)
    batch_rows: dict[str, dict] = {}
    for row in rows:
        job_id = row["JobID"]
        if not job_id.endswith(".batch"):
            continue
        parent_id = job_id.split(".")[0]
        if parent_id in completed_ids:
            batch_rows[parent_id] = row

    # Index by job name -> list of (job_id, merged record) sorted by job_id (ascending = chronological)
    index: dict[str, list[dict]] = {}
    for job_id, parent in parent_rows.items():
        name = parent["JobName"]
        batch = batch_rows.get(job_id, {})
        record = {
            "job_id": job_id,
            "state": parent["State"],
            "start": parent["Start"],
            "end": parent["End"],
            "elapsed": parent["Elapsed"],
            "elapsed_s": elapsed_to_seconds(parent["Elapsed"]),
            "elapsed_raw_s": int(parent["ElapsedRaw"]) if parent["ElapsedRaw"].isdigit() else None,
            "ncpus": int(parent["NCPUS"]) if parent["NCPUS"].isdigit() else None,
            "nnodes": int(parent["NNodes"]) if parent["NNodes"].isdigit() else None,
            "cpu_time_raw_s": int(parent["CPUTimeRAW"]) if parent["CPUTimeRAW"].isdigit() else None,
            "node": parent["NodeList"],
            "exit_code": parent["ExitCode"],
            # Memory from batch step
            "peak_rss_bytes": parse_max_rss_kb(batch.get("MaxRSS", "")),
            "peak_vmsize_bytes": parse_max_rss_kb(batch.get("MaxVMSize", "")),
            "avg_rss_bytes": parse_max_rss_kb(batch.get("AveRSS", "")),
        }
        # Derived: CPU efficiency (wall * ncpus vs cpu_time)
        if record["elapsed_s"] and record["ncpus"] and record["cpu_time_raw_s"]:
            total_possible = record["elapsed_s"] * record["ncpus"]
            record["cpu_efficiency"] = round(record["cpu_time_raw_s"] / total_possible, 4) if total_possible else None
        else:
            record["cpu_efficiency"] = None

        index.setdefault(name, []).append(record)

    # Sort each list by job_id numerically so the last entry is the most recent
    for name in index:
        index[name].sort(key=lambda r: int(r["job_id"]))

    return index


def get_current_user() -> str:
    result = subprocess.run(["whoami"], capture_output=True, text=True, check=True)
    return result.stdout.strip()


def main() -> None:
    user = get_current_user()
    print(f"Querying sacct for user {user} ...")
    rows = run_sacct(user)
    job_index = build_job_index(rows)
    print(f"Found {len(job_index)} distinct completed job names.")

    written = 0
    missing = 0
    for timing_path in sorted(TOOLS_ROOT.glob("*/**/timing.json")):
        tool_dir = timing_path.parent
        # Infer tool name from grandparent dir name and dataset from parent dir name
        parts = tool_dir.relative_to(TOOLS_ROOT).parts
        tool_name = parts[0]
        dataset_name = parts[-1]
        job_name = f"{tool_name}-{dataset_name}"

        records = job_index.get(job_name)
        if not records:
            print(f"  [MISS] No completed job found for {job_name}")
            missing += 1
            continue

        # Take the most recent completed run
        stats = records[-1]
        out_path = tool_dir / "job-stats.json"
        out_path.write_text(json.dumps(stats, indent=2) + "\n")
        print(f"  [OK]   {job_name}  ->  {out_path.relative_to(TOOLS_ROOT)}")
        written += 1

    print(f"\nDone: {written} written, {missing} missing.")


if __name__ == "__main__":
    main()
