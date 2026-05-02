# Downloading archives from Polybox / public shares

## Direct-download rule

Always use the direct download endpoint:

- Share page: `https://polybox.ethz.ch/index.php/s/<TOKEN>`
- Download URL: `https://polybox.ethz.ch/index.php/s/<TOKEN>/download`

Example:

```bash
curl -L --fail --continue-at - --output out.fastq.gz \
  "https://polybox.ethz.ch/index.php/s/<TOKEN>/download"
```

## Validation rule

Do not assume a public-share link returns the real payload. Validate that the download starts correctly:

```bash
curl -I -L "https://polybox.ethz.ch/index.php/s/<TOKEN>/download"
file out
xxd -l 16 out
```

Expected magic bytes:

- `50 4b 03 04` -> ZIP
- `1f 8b` -> GZIP / BGZF
- `fd 37 7a 58 5a 00` -> XZ
- `42 5a 68` -> BZIP2

## Validated links for this repository

### Favalli

- Lane 1 FASTQ: `https://polybox.ethz.ch/index.php/s/HbtDGTznwWGxwKr/download`
  Saved as `supporting_material/experiments/favalli/20190812.A-1907_NF2GB2_s1_R1.fastq.gz`
- Lane 2 FASTQ: `https://polybox.ethz.ch/index.php/s/9bX3YnzyDy8ez8g/download`
  Saved as `supporting_material/experiments/favalli/20190812.A-1907_NF2GB2_s2_R1.fastq.gz`
- Lane 1 FASTA: `https://polybox.ethz.ch/index.php/s/9xZkcW5jWtLHiFo/download`
  Saved as `supporting_material/experiments/favalli/20190812.A-1907_NF2GB2_s1_R1.fasta.gz`
- Lane 2 FASTA: `https://polybox.ethz.ch/index.php/s/XY2ieRTwdM337Jc/download`
  Saved as `supporting_material/experiments/favalli/20190812.A-1907_NF2GB2_s2_R1.fasta.gz`
- Published evaluation table: `https://polybox.ethz.ch/index.php/s/cgxBWkHnpkT6rBi/download`
  Saved as `supporting_material/experiments/favalli/published/1907_NF2GB2_s1_R1_260424JS_2026_4_24_16_20_51_eval.txt`

Run:

```bash
cd supporting_material/experiments/favalli
bash download.sh
```

### Johanna

- FASTQ: `https://polybox.ethz.ch/index.php/s/bXTtkZ2666eZoDT/download`
  Saved as `supporting_material/experiments/johanna/407991_1-2601_JPAG_L1c_NF_501701_S3_R1_001.fastq.gz`

Run:

```bash
cd supporting_material/experiments/johanna
bash download.sh
```

### Pure-DEL

- Lane 1 FASTQ: `https://polybox.ethz.ch/index.php/s/H8YB6wdxpwWgzLB/download`
  Saved as `supporting_material/experiments/pure-del/328321_1-230906_MK_TSD501701_S1_R1_001_SP_fl1.fastq.gz`
- Lane 2 FASTQ: `https://polybox.ethz.ch/index.php/s/wK6t8qjDTay8rAg/download`
  Saved as `supporting_material/experiments/pure-del/328321_2-230906_MK_TSD504704_S2_R1_001_UP_fl2.fastq.gz`
- Published counts ZIP: `https://zenodo.org/records/11122901/files/Keller_Petrov_etal_Supplementary_Material2_Sequencing_Data.zip?download=1`
  Extracted into `supporting_material/experiments/pure-del/published/`

Run:

```bash
cd supporting_material/experiments/pure-del
bash download.sh
```

## Notes on extraction

- The Polybox FASTQ links validated above return real `application/gzip` payloads with BGZF-compatible gzip magic bytes.
- The Favalli published-data link returns a real text payload and is saved into `supporting_material/experiments/favalli/published/`.
- The Zenodo published-data link returns a real ZIP archive containing the selection tables. The `pure-del/download.sh` helper unzips it and copies the `selection_*_.txt` files into `published/`.
