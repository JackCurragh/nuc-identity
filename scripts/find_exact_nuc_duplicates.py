#!/usr/bin/env python3
import argparse, gzip, hashlib
from pathlib import Path

def iter_fa(p: Path):
    op=gzip.open if str(p).endswith('.gz') else open
    with op(p,'rt') as fh:
        h=None; parts=[]
        for line in fh:
            if line.startswith('>'):
                if h is not None:
                    yield h, ''.join(parts)
                h=line[1:].strip().split()[0]; parts=[]
            else:
                parts.append(line.strip())
        if h is not None:
            yield h, ''.join(parts)

def parse_header(h:str):
    pr=h.split('|');
    return (pr[0] if pr else h, pr[1] if len(pr)>1 else '', pr[2] if len(pr)>2 else '', pr[3] if len(pr)>3 else '')

def main():
    ap=argparse.ArgumentParser(description='Exact duplicate nucleotide sequences across loci (hash-based)')
    ap.add_argument('--fasta', default='results/nuc/transcripts.fa')
    ap.add_argument('--out-pairs', default='results/nuc_clusters/exact_hash_cross_gene_pairs.tsv')
    ap.add_argument('--out-summary', default='results/nuc_clusters/exact_hash_groups.tsv')
    a=ap.parse_args()
    Path(a.out_pairs).parent.mkdir(parents=True, exist_ok=True)

    groups={}
    n=0
    for hdr, seq in iter_fa(Path(a.fasta)):
        n+=1
        h=hashlib.md5(seq.encode('utf-8')).hexdigest()
        tid,gid,gname,loc=parse_header(hdr)
        groups.setdefault(h, []).append((tid,gid,gname,loc))
    out_pairs=[]; out_groups=[]
    for h, items in groups.items():
        genes={it[1] for it in items if it[1]}
        if len(items)>1 and len(genes)>1:
            out_groups.append((h, len(items), len(genes)))
            for i in range(len(items)):
                for j in range(i+1, len(items)):
                    t1,g1,n1,l1=items[i]; t2,g2,n2,l2=items[j]
                    if g1!=g2:
                        out_pairs.append((h,t1,g1,n1,l1,t2,g2,n2,l2))
    with open(a.out_pairs,'w') as fo:
        fo.write('hash\tt1\tg1\tn1\tloc1\tt2\tg2\tn2\tloc2\n')
        for r in out_pairs: fo.write('\t'.join(map(str,r))+'\n')
    with open(a.out_summary,'w') as fo:
        fo.write('hash\tn_transcripts\tn_genes\n')
        for h,nt,ng in sorted(out_groups, key=lambda x:(-x[1],-x[2])):
            fo.write(f"{h}\t{nt}\t{ng}\n")
    print(f"Exact cross-gene groups: {len(out_groups)}; pairs: {len(out_pairs)}")

if __name__=='__main__':
    main()

