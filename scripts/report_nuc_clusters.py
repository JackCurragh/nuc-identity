#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd

def load_pairs(tsv: Path) -> pd.DataFrame:
    if not tsv.exists():
        return pd.DataFrame(columns=['rep','member'])
    df = pd.read_csv(tsv, sep='\t', header=None, names=['rep','member'], dtype=str)
    return df

def load_meta(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=['transcript_id','gene_id','gene_name','chrom','start','end','strand','mode','length'])
    df = pd.read_csv(path, sep='\t', dtype=str)
    for c in ['transcript_id','gene_id','gene_name','chrom','start','end','strand','mode','length']:
        if c not in df.columns:
            df[c] = ''
    return df[['transcript_id','gene_id','gene_name','chrom','start','end','strand','mode','length']]

def build_members(pairs: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    m = pairs.merge(meta, how='left', left_on='member', right_on='transcript_id')
    reps = pairs[['rep']].drop_duplicates().rename(columns={'rep':'transcript_id'})
    reps = reps.merge(meta, on='transcript_id', how='left')
    reps['member'] = reps['transcript_id']
    reps['rep'] = reps['transcript_id']
    cols = ['rep','member','transcript_id','gene_id','gene_name','chrom','start','end','strand','mode','length']
    return pd.concat([m[cols], reps[cols]], ignore_index=True).drop_duplicates()

def summarize(members: pd.DataFrame) -> pd.DataFrame:
    agg = members.groupby('rep').agg(
        n_transcripts=('member','size'),
        n_genes=('gene_id', lambda s: s.dropna().astype(str).nunique()),
        gene_names=('gene_name', lambda s: '|'.join(sorted(set([x for x in s.dropna().astype(str) if x]))[:60])),
        chroms=('chrom', lambda s: '|'.join(sorted(set([x for x in s.dropna().astype(str) if x]))[:20]))
    ).reset_index()
    return agg

def cross_gene_views(members: pd.DataFrame, outdir: Path, prefix: str):
    if members.empty:
        return
    multi = members.groupby('rep').filter(lambda df: df['gene_id'].dropna().astype(str).nunique() > 1)
    if multi.empty:
        return
    multi.to_csv(outdir / f'{prefix}_cross_gene_members.tsv', sep='\t', index=False)
    pairs = []
    for rep, grp in multi.groupby('rep'):
        rows = grp[['transcript_id','gene_id','gene_name','chrom','start','end','strand']].drop_duplicates()
        rows = rows.sort_values('transcript_id')
        arr = rows.to_dict('records')
        for i in range(len(arr)):
            for j in range(i+1, len(arr)):
                if arr[i]['gene_id'] != arr[j]['gene_id']:
                    pairs.append({
                        'rep': rep,
                        't1': arr[i]['transcript_id'], 'g1': arr[i]['gene_id'], 'n1': arr[i]['gene_name'], 'loc1': f"{arr[i]['chrom']}:{arr[i]['start']}-{arr[i]['end']}:{arr[i]['strand']}",
                        't2': arr[j]['transcript_id'], 'g2': arr[j]['gene_id'], 'n2': arr[j]['gene_name'], 'loc2': f"{arr[j]['chrom']}:{arr[j]['start']}-{arr[j]['end']}:{arr[j]['strand']}",
                    })
    if pairs:
        pd.DataFrame(pairs).to_csv(outdir / f'{prefix}_cross_gene_pairs.tsv', sep='\t', index=False)

def write_html(copy_sum: pd.DataFrame, out_html: Path):
    html = []
    html.append('<html><head><meta charset="utf-8"><style>body{font:14px system-ui} table{border-collapse:collapse} td,th{border:1px solid #ddd;padding:4px} h2{margin-top:1.5em}</style></head><body>')
    html.append('<h1>Nucleotide Similarity Clusters</h1>')
    if not copy_sum.empty:
        html.append(f'<h2>Near-identical (copy-like)</h2>')
        html.append(f'<p>Total clusters: <b>{copy_sum.shape[0]}</b></p>')
        html.append(copy_sum.sort_values(['n_transcripts','n_genes'], ascending=False).head(100).to_html(index=False, escape=False))
    html.append('</body></html>')
    out_html.write_text('\n'.join(html))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--meta', default='results/nuc/transcript_meta.tsv')
    ap.add_argument('--copylike', default='results/nuc_clusters/copylike_cluster.tsv')
    ap.add_argument('--outdir', default='results/nuc_clusters')
    args = ap.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    pairs_copy = load_pairs(Path(args.copylike))
    meta = load_meta(Path(args.meta))
    if not pairs_copy.empty:
        members_copy = build_members(pairs_copy, meta)
        members_copy.to_csv(outdir / 'copylike_members.tsv', sep='\t', index=False)
        copy_sum = summarize(members_copy)
        copy_sum.to_csv(outdir / 'copylike_summary.tsv', sep='\t', index=False)
        cross_gene_views(members_copy, outdir, prefix='copylike')
    else:
        copy_sum = pd.DataFrame()
    write_html(copy_sum, out_html=outdir / 'nucleotide_clusters_report.html')
    print('Wrote', outdir)

if __name__ == '__main__':
    main()

