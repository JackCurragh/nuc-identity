SHELL := /bin/bash
PY := python3
VENV := .venv
GENCODE_RELEASE ?=

.PHONY: env fetch-gencode build-cdna build-genome cluster cluster-cdhit report all-cdna all-genome exact-only

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
	# Prefer MMseqs2; fall back to CD-HIT-EST if absent
	if command -v mmseqs >/dev/null 2>&1; then \
	  ./scripts/cluster_nucleotides.sh results/nuc/transcripts.fa results/nuc_clusters $${MIN_ID:-0.98} $${MIN_COV:-0.9}; \
	else \
	  $(MAKE) cluster-cdhit; \
	fi

cluster-cdhit: results/nuc/transcripts.fa
	@which cd-hit-est >/dev/null 2>&1 || { echo "cd-hit-est not found. Install via mamba/conda: mamba install -c bioconda cd-hit" >&2; exit 2; }
	mkdir -p results/nuc_clusters
	cd-hit-est -i results/nuc/transcripts.fa -o results/nuc_clusters/copylike -c $${MIN_ID:-0.98} -aS $${MIN_COV:-0.9} -aL $${MIN_COV:-0.9} -T $${THREADS:-16} -M $${MEM_MB:-0}
	$(PY) scripts/cdhit_clstr_to_pairs.py results/nuc_clusters/copylike.clstr results/nuc_clusters/copylike_cluster.tsv

report:
	$(PY) scripts/report_nuc_clusters.py --meta results/nuc/transcript_meta.tsv --copylike results/nuc_clusters/copylike_cluster.tsv --exact results/nuc_clusters/exact_cluster.tsv --outdir results/nuc_clusters

exact-only:
	$(PY) scripts/find_exact_nuc_duplicates.py --fasta results/nuc/transcripts.fa --out-pairs results/nuc_clusters/exact_hash_cross_gene_pairs.tsv --out-summary results/nuc_clusters/exact_hash_groups.tsv

all-cdna: fetch-gencode build-cdna cluster report
all-genome: fetch-gencode build-genome cluster report
