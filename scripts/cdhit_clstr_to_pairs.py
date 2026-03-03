#!/usr/bin/env python3
import sys
import re

# Convert CD-HIT .clstr to rep\tmember TSV. Assumes headers have no spaces.

def main():
    if len(sys.argv) < 3:
        print("Usage: cdhit_clstr_to_pairs.py in.clstr out.tsv", file=sys.stderr)
        sys.exit(2)
    inp, outp = sys.argv[1], sys.argv[2]
    rep = None
    pairs = []
    with open(inp, 'r') as f:
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            if line.startswith('>Cluster'):
                rep = None
                continue
            # Example entry: "0\t123nt, >seq_id... *" or "1\t120nt, >seq2... at 99.7%"
            m = re.search(r'>\s*([^\.\s]+)\.\.', line)
            if not m:
                continue
            sid = m.group(1)
            is_rep = line.strip().endswith('*')
            if is_rep:
                rep = sid
            if rep is None:
                # First member seen before a representative; treat as rep
                rep = sid
            else:
                # Emit rep, member (skip self if the starred line is encountered)
                if not is_rep:
                    pairs.append((rep, sid))
    with open(outp, 'w') as fo:
        for r, m in pairs:
            fo.write(f"{r}\t{m}\n")

if __name__ == '__main__':
    main()

