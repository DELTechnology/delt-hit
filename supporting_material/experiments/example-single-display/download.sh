#!/usr/bin/env bash

set -euo pipefail

# Run after activating the delt-hit Conda environment.
python - <<'PYTHON'
import json
import shutil
from pathlib import Path
from urllib.request import urlopen
from tqdm import tqdm

article_id = "31198468"
file_id = "61487743"
out = Path("campaign.fastq.gz")

if out.is_file() and out.stat().st_size:
    print(f"{out} already exists, skipping")
else:
    with urlopen(f"https://api.figshare.com/v2/articles/{article_id}", timeout=60) as response:
        meta = json.load(response)
    target = next(f for f in meta["files"] if str(f["id"]) == file_id)
    print("File:", target["name"])
    print("Size:", target["size"])
    partial = out.with_suffix(out.suffix + ".part")
    with urlopen(target["download_url"], timeout=120) as response, partial.open("wb") as dest:
        with tqdm.wrapattr(response, "read", total=target["size"], desc=out.name) as source:
            shutil.copyfileobj(source, dest, length=1024 * 1024)
    if partial.stat().st_size != target["size"]:
        raise RuntimeError("Downloaded file size differs from the archive metadata")
    partial.replace(out)
    print(f"Saved to {out}")
PYTHON
