#!/usr/bin/env bash
set -euo pipefail
IN=${1:-results/nuc/transcripts.fa}
OUTDIR=${2:-results/nuc_clusters}
MIN_ID=${3:-0.98}
MIN_COV=${4:-0.9}
THREADS=${THREADS:-16}
mkdir -p "$OUTDIR"
command -v mmseqs >/dev/null 2>&1 || { echo "mmseqs not found" >&2; exit 2; }
mmseqs createdb "$IN" "$OUTDIR/db" >/dev/null
mmseqs cluster "$OUTDIR/db" "$OUTDIR/clu" "$OUTDIR/tmp" --min-seq-id "$MIN_ID" -c "$MIN_COV" --cov-mode 1 --threads "$THREADS" >/dev/null
mmseqs createtsv "$OUTDIR/db" "$OUTDIR/db" "$OUTDIR/copylike_cluster.tsv" "$OUTDIR/clu" >/dev/null || true
mmseqs createseqfiledb "$OUTDIR/db" "$OUTDIR/clu" "$OUTDIR/seqdb" >/dev/null || true
mmseqs result2flat "$OUTDIR/db" "$OUTDIR/db" "$OUTDIR/seqdb" "$OUTDIR/copylike_all_seqs.fasta" >/dev/null || true
rm -rf "$OUTDIR/db"* "$OUTDIR/tmp" "$OUTDIR/seqdb"* || true
echo "Clusters written to $OUTDIR"
