#!/usr/bin/env bash
set -euo pipefail
rel="${1:-}"; outdir="data/gencode"; mkdir -p "$outdir"
base="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human"
try_release(){ local r="$1"; local dir="$base/release_${r}"; curl -fsI "$dir/" >/dev/null || return 1
  echo "$r" > "$outdir/RELEASE.txt"
  curl -fsSL "$dir/gencode.v${r}.primary_assembly.annotation.gff3.gz" -o "$outdir/gencode.v${r}.primary_assembly.annotation.gff3.gz"
  curl -fsSL "$dir/GRCh38.primary_assembly.genome.fa.gz" -o "$outdir/GRCh38.primary_assembly.genome.fa.gz"
}
if [[ -n "$rel" ]]; then try_release "$rel"; else
  html=$(curl -fsSL "$base/latest_release/")
  ver=$(echo "$html" | sed -n 's/.*gencode\.v\([0-9][0-9]*\)\.primary_assembly\.annotation\.gff3\.gz.*/\1/p' | head -n1)
  [[ -n "$ver" ]] && try_release "$ver" || for v in {50..38}; do try_release "$v" && break; done
fi
ls -lh "$outdir" | cat
