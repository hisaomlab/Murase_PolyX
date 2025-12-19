#!/usr/bin/env python3

import argparse
import gzip
from collections import defaultdict

AA_LIST = "ACDEFGHIKLMNPQRSTVWY"


def open_fasta(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def parse_fasta_longest(path):
    seqs = {}
    gid = None
    buf = []

    with open_fasta(path) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if gid:
                    s = "".join(buf)
                    if gid not in seqs or len(s) > len(seqs[gid]):
                        seqs[gid] = s
                gid = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)

        if gid:
            s = "".join(buf)
            if gid not in seqs or len(s) > len(seqs[gid]):
                seqs[gid] = s

    return seqs


def analyze_polyx(seqs, threshold):
    max_len = {aa: 0 for aa in AA_LIST}
    hits = defaultdict(list)

    for gid, seq in seqs.items():
        prev, run = None, 0
        for aa in seq:
            if aa == prev:
                run += 1
            else:
                if prev in max_len:
                    if run > max_len[prev]:
                        max_len[prev] = run
                    if run >= threshold:
                        hits[prev].append((gid, run))
                prev, run = aa, 1
        if prev in max_len:
            if run > max_len[prev]:
                max_len[prev] = run
            if run >= threshold:
                hits[prev].append((gid, run))

    return max_len, hits


def write_summary(max_len, hits, n_genes, out):
    out.write("AA\tmax_polyX_length\t#poly10_genes\t#poly10_per_gene\n")
    for aa in AA_LIST:
        genes = {g for g, _ in hits.get(aa, [])}
        out.write(f"{aa}\t{max_len[aa]}\t{len(genes)}\t{len(genes)/n_genes:.5f}\n")


def write_fasta(seqs, hits, out):
    seen = set()
    for aa in AA_LIST:
        for gid, run in hits.get(aa, []):
            if (gid, aa) in seen:
                continue
            seen.add((gid, aa))
            out.write(f">{gid}\t{aa * run}\n{seqs[gid]}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fasta")
    ap.add_argument("--threshold", type=int, default=10)
    ap.add_argument("--label", default=None, help="label to avoid overwrite")
    args = ap.parse_args()

    label = f"_{args.label}" if args.label else ""
    summary_out = f"polyx_summary{label}.tsv"
    fasta_out = f"poly10{label}.fa"

    seqs = parse_fasta_longest(args.fasta)
    max_len, hits = analyze_polyx(seqs, args.threshold)

    with open(summary_out, "w") as o:
        write_summary(max_len, hits, len(seqs), o)
    with open(fasta_out, "w") as o:
        write_fasta(seqs, hits, o)


if __name__ == "__main__":
    main()

