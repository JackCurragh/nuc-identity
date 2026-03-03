#!/usr/bin/env python3
import argparse, gzip
from pathlib import Path
from Bio import SeqIO

def parse_gtf_attrs(s:str):
    out={}
    for kv in s.strip().strip(';').split(';'):
        kv=kv.strip();
        if not kv: continue
        if ' ' in kv:
            k,v=kv.split(' ',1); out[k.strip()]=v.strip().strip('"')
    return out

def build_tx_meta(gtf:Path):
    tx2gene={}; tx2loc={}
    op=gzip.open if str(gtf).endswith('.gz') else open
    with op(gtf,'rt') as f:
        for line in f:
            if not line or line.startswith('#'): continue
            p=line.rstrip('\n').split('\t')
            if len(p)<9: continue
            chrom,source,feat,start,end,score,strand,frame,attrs=p
            if feat!='transcript': continue
            ad=parse_gtf_attrs(attrs)
            tid=ad.get('transcript_id'); gid=ad.get('gene_id'); gname=ad.get('gene_name')
            if not tid: continue
            tx2gene[tid]=(gid or '', gname or '')
            tx2loc[tid]=(chrom,strand,int(start),int(end))
    return tx2gene, tx2loc

def main():
    ap=argparse.ArgumentParser(description='Build spliced FASTA from GENCODE transcripts + GTF metadata')
    ap.add_argument('--gtf', required=True); ap.add_argument('--fa', required=True)
    ap.add_argument('--out-fasta', default='results/nuc/transcripts.fa')
    ap.add_argument('--out-meta', default='results/nuc/transcript_meta.tsv')
    a=ap.parse_args()
    Path(a.out_fasta).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out_meta).parent.mkdir(parents=True, exist_ok=True)
    tx2gene, tx2loc=build_tx_meta(Path(a.gtf))
    n=0
    op=gzip.open if str(a.fa).endswith('.gz') else open
    with op(a.fa,'rt') as fh, open(a.out_fasta,'w') as fo, open(a.out_meta,'w') as mo:
        mo.write('transcript_id\tgene_id\tgene_name\tchrom\tstart\tend\tstrand\tmode\tlength\n')
        for rec in SeqIO.parse(fh,'fasta'):
            tid=rec.id.split('|')[0]
            gid,gname=tx2gene.get(tid,('',''))
            chrom,strand,start,end=tx2loc.get(tid,('', '+', 0, 0))
            seq=str(rec.seq)
            header=f"{tid}|{gid}|{gname}|{chrom}:{start}-{end}:{strand}|len={len(seq)}"
            fo.write(f">{header}\n{seq}\n")
            mo.write('\t'.join([tid,gid,gname,chrom,str(start),str(end),strand,'cdna',str(len(seq))])+'\n')
            n+=1
    print(f"Wrote {n} sequences -> {a.out_fasta}")

if __name__=='__main__':
    main()

