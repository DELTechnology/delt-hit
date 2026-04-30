# Downloading archives from Polybox / public shares

## 1. Always use the direct download endpoint
Wrong (HTML page):
https://polybox.ethz.ch/index.php/s/<TOKEN>

Right (actual archive):
https://polybox.ethz.ch/index.php/s/<TOKEN>/download

Download:
wget -c --content-disposition "https://polybox.ethz.ch/index.php/s/<TOKEN>/download"

## 2. Check what the file actually is
file download

If unclear:
xxd -l 16 download

Magic bytes:
50 4b 03 04  -> zip
1f 8b        -> gzip
fd 37 7a...  -> xz
42 5a 68     -> bzip2

## 3. Extract based on type
ZIP:
unzip download

TAR:
tar -xf download

TAR.GZ:
tar -xzf download

TAR.XZ:
tar -xJf download

## Rule of thumb
Never trust file extension.
Trust magic bytes (file or xxd).