#!/usr/bin/env python3
import argparse, gzip
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Tuple
from Bio import SeqIO

def open_text(p: Path):
    return gzip.open(p, 'rt') if str(p).endswith('.gz') else open(p, 'rt')

def parse_attrs(s: str) -> Dict[str,str]:
    out={}
    for kv in s.strip().strip(';').split(';'):
        kv=kv.strip();
        if not kv: continue
        if '=' in kv:
            k,v=kv.split('=',1); out[k.strip()]=v.strip()
        elif ' ' in kv:
            k,v=kv.split(' ',1); out[k.strip()]=v.strip().strip('"')
    return out

def load_genome(fa: Path) -> Dict[str,str]:
    seqs={}
    opener = gzip.open if str(fa).endswith('.gz') else open
    with opener(fa,'rt') as h:
        for rec in SeqIO.parse(h,'fasta'):
            seqs[rec.id]=str(rec.seq)
    return seqs

def revcomp(s:str)->str:
    return s.translate(str.maketrans('ACGTacgtNn','TGCAtgcaNn'))[::-1]

def extract_segments(genome:Dict[str,str], chrom:str, segs:List[Tuple[int,int]], strand:str)->str:
    g=genome.get(chrom); 
    if g is None: return ''
    parts=[]
    for s,e in segs:
        s0=max(1,s)-1; e0=min(len(g),e)
        if e0<=s0: continue
        parts.append(g[s0:e0])
    seq=''.join(parts)
    return revcomp(seq) if strand=='-' else seq

def extract_from_gff(gff: Path, genome: Dict[str,str], mode:str):
    tx_chrom={}; tx_strand={}; tx_gene={}; tx_segs=defaultdict(list)
    with open_text(gff) as f:
        for line in f:
            if not line or line.startswith('#'): continue
            p=line.rstrip('\n').split('\t');
            if len(p)<9: continue
            chrom,source,feat,start,end,score,strand,frame,attrs=p
            s=int(start); e=int(end); ad=parse_attrs(attrs)
            if feat.lower() in ('mrna','transcript'):
                tid=ad.get('ID') or ad.get('transcript_id') or ''
                if tid:
                    tx_chrom[tid]=chrom; tx_strand[tid]=strand
                    tx_gene[tid]=(ad.get('gene_id',''), ad.get('gene_name', ad.get('Name','')))
                continue
            if (mode=='cdna' and feat.lower()=='exon') or (mode=='cds' and feat=='CDS'):
                parent=ad.get('Parent') or ad.get('transcript_id') or ''
                if parent:
                    tx_chrom.setdefault(parent, chrom); tx_strand.setdefault(parent, strand)
                    if ad.get('gene_id') or ad.get('gene_name'):
                        tx_gene.setdefault(parent,(ad.get('gene_id',''), ad.get('gene_name','')))
                    tx_segs[parent].append((s,e))
    fasta=[]; meta=[]
    for tid,segs in tx_segs.items():
        segs.sort(key=lambda x:x[0]); chrom=tx_chrom.get(tid,''); strand=tx_strand.get(tid,'+')
        gid,gname=tx_gene.get(tid,('',''))
        seq=extract_segments(genome, chrom, segs, strand)
        if not seq: continue
        span_start=min(s for s,_ in segs); span_end=max(e for _,e in segs)
        header=f"{tid}|{gid}|{gname}|{chrom}:{span_start}-{span_end}:{strand}|len={len(seq)}"
        fasta.append(f">{header}\n{seq}\n")
        meta.append('\t'.join([tid,gid,gname,chrom,str(span_start),str(span_end),strand,mode,str(len(seq))]))
    return fasta, meta

def main():
    ap=argparse.ArgumentParser(description='Extract spliced cDNA/CDS from GFF/GTF + genome')
    ap.add_argument('--gff', required=True); ap.add_argument('--genome', required=True)
    ap.add_argument('--mode', choices=['cdna','cds'], default='cdna')
    ap.add_argument('--out-fasta', default='results/nuc/transcripts.fa')
    ap.add_argument('--out-meta', default='results/nuc/transcript_meta.tsv')
    a=ap.parse_args()
    Path(a.out_fasta).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out_meta).parent.mkdir(parents=True, exist_ok=True)
    genome=load_genome(Path(a.genome))
    fa, meta=extract_from_gff(Path(a.gff), genome, a.mode)
    with open(a.out_fasta,'w') as fo: fo.write(''.join(fa))
    with open(a.out_meta,'w') as mo:
        mo.write('transcript_id\tgene_id\tgene_name\tchrom\tstart\tend\tstrand\tmode\tlength\n')
        mo.write('\n'.join(meta)+('\n' if meta else ''))
    print(f"Wrote {len(fa)} sequences -> {a.out_fasta}")

if __name__=='__main__':
    main()

