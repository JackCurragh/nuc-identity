# nuc-identity

Standalone pipeline to find exact and near-identical nucleotide sequences across loci from a genome + GFF3/GTF or from GENCODE transcripts, then report cross-gene duplicates.

Features
- Fetch GENCODE data (GTF, transcripts.fa.gz; optional genome+GFF3)
- Build spliced cDNA (or CDS) sequences with metadata
- Cluster nucleotide sequences (MMseqs2): exact and near-identical
- Exact-duplicate scan (hash-based) as a fast approximation
- Reports: cross-gene duplicate groups and pairs + HTML summary

Quickstart (local)
```bash
make env
make fetch-gencode            # downloads v49 by default
make all-cdna                 # transcripts -> clusters -> report
# OR from genome+GFF3 splicing
make all-genome NUC_MODE=cdna # or NUC_MODE=cds
```

HPC (SLURM) example
```bash
# Use your mmseqs2 module and set threads
module load mmseqs2
export THREADS=32
sbatch -J nucid -c 32 --mem=120G --time=24:00:00 \
  --wrap "cd $PWD && make -j 1 all-cdna"
```

Outputs
- `data/gencode/` – downloads/cache
- `results/nuc/transcripts.fa`, `transcript_meta.tsv` – spliced sequences + metadata
- `results/nuc_clusters/` – `exact_cluster.tsv`, `copylike_cluster.tsv`, cross-gene TSVs, HTML summary

Requirements
- Python 3.10+, Biopython, pandas
- MMseqs2 (for near-identical clustering). Exact duplicates work without it.

Notes
- BLASTn is not splice-aware for genome↔cDNA; this pipeline reconstructs spliced sequences first.
- Tune identity/coverage in `scripts/cluster_nucleotides.sh` or pass env vars.

