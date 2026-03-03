SHELL := /bin/bash
PY := python3
VENV := .venv
GENCODE_RELEASE ?=

.PHONY: env fetch-gencode build-cdna build-genome cluster report all-cdna all-genome exact-only

help:
	@echo "Targets: env, fetch-gencode, build-cdna, build-genome, cluster, report, all-cdna, all-genome, exact-only"

env:
	python3 -m venv $(VENV)
	. $(VENV)/bin/activate && pip install --upgrade pip && pip install -r requirements.txt

fetch-gencode:
	./scripts/download_gencode.sh $(GENCODE_RELEASE)
	./scripts/download_gencode_genome.sh $(GENCODE_RELEASE) || true

build-cdna: data/gencode/gencode.v$(shell cat data/gencode/RELEASE.txt 2>/dev/null).annotation.gtf.gz data/gencode/gencode.v$(shell cat data/gencode/RELEASE.txt 2>/dev/null).transcripts.fa.gz
	$(PY) scripts/prepare_transcripts_from_gencode.py --gtf $< --fa $(word 2,$^) --out-fasta results/nuc/transcripts.fa --out-meta results/nuc/transcript_meta.tsv

build-genome:
	$(PY) scripts/extract_spliced_sequences.py --gff data/gencode/gencode.v$(shell cat data/gencode/RELEASE.txt 2>/dev/null).primary_assembly.annotation.gff3.gz --genome data/gencode/GRCh38.primary_assembly.genome.fa.gz --mode $${NUC_MODE:-cdna} --out-fasta results/nuc/transcripts.fa --out-meta results/nuc/transcript_meta.tsv

cluster: results/nuc/transcripts.fa
	./scripts/cluster_nucleotides.sh results/nuc/transcripts.fa results/nuc_clusters $${MIN_ID:-0.98} $${MIN_COV:-0.9}

report:
	$(PY) scripts/report_nuc_clusters.py --meta results/nuc/transcript_meta.tsv --copylike results/nuc_clusters/copylike_cluster.tsv --exact results/nuc_clusters/exact_cluster.tsv --outdir results/nuc_clusters

exact-only:
	$(PY) scripts/find_exact_nuc_duplicates.py --fasta results/nuc/transcripts.fa --out-pairs results/nuc_clusters/exact_hash_cross_gene_pairs.tsv --out-summary results/nuc_clusters/exact_hash_groups.tsv

all-cdna: fetch-gencode build-cdna cluster report
all-genome: fetch-gencode build-genome cluster report

