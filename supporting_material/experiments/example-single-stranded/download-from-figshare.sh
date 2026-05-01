ARTICLE_ID=31198468
FILE_ID=61487743
OUT=campaign.fastq.gz

uv run python - <<'PY'
import requests

article_id = "31198468"
file_id = "61487743"
out = "campaign.fastq.gz"

api_url = f"https://api.figshare.com/v2/articles/{article_id}"
meta = requests.get(api_url, timeout=60).json()

files = meta["files"]
target = next(f for f in files if str(f["id"]) == file_id)

print("File:", target["name"])
print("Size:", target["size"])
print("Download URL:", target["download_url"])

with requests.get(target["download_url"], stream=True, timeout=120) as r:
    r.raise_for_status()
    with open(out, "wb") as f:
        for chunk in r.iter_content(chunk_size=1024 * 1024):
            if chunk:
                f.write(chunk)

print(f"Saved to {out}")
PY