#!/usr/bin/env python3

import argparse
import base64
import gzip
import io
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def infer_sample_id(path):
    name = path.name
    for suffix in [".meth.pairs.gz", ".meth.pairs", ".pairs.gz", ".pairs"]:
        if name.endswith(suffix):
            return name[:-len(suffix)]
    return path.stem


def fig_to_base64():
    buf = io.BytesIO()
    plt.tight_layout()
    plt.savefig(buf, format="png", dpi=150)
    plt.close()
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pairs", required=True, type=Path)
    parser.add_argument("--out", default="pair_level_qc.html", type=Path)
    parser.add_argument("--cis-distance", default=10000, type=int)
    args = parser.parse_args()

    sample_id = infer_sample_id(args.pairs)

    columns = None
    idx = None

    total_pairs = 0
    cis_pairs = 0
    trans_pairs = 0
    cis_long_pairs = 0

    cis_distances = []
    distance_bins = Counter()
    fragment_lengths = Counter()
    read_lengths = Counter()

    with open_text(args.pairs) as f:
        for line in f:
            if line.startswith("#columns:"):
                columns = line.strip().split(":", 1)[1].split()
                idx = {c: i for i, c in enumerate(columns)}
                required = ["chrom1", "pos1", "chrom2", "pos2", "pos51", "pos52", "pos31", "pos32", "seq1", "seq2"]
                missing = [x for x in required if x not in idx]
                if missing:
                    raise ValueError(f"Missing required columns: {missing}")
                continue

            if line.startswith("#"):
                continue

            if idx is None:
                raise ValueError("Missing #columns header in pairs file")

            fields = line.rstrip("\n").split("\t")
            total_pairs += 1

            chrom1 = fields[idx["chrom1"]]
            chrom2 = fields[idx["chrom2"]]
            pos1 = int(fields[idx["pos1"]])
            pos2 = int(fields[idx["pos2"]])

            if chrom1 == chrom2:
                cis_pairs += 1
                dist = abs(pos2 - pos1)
                cis_distances.append(dist)

                if dist >= args.cis_distance:
                    cis_long_pairs += 1

                if dist < 1_000:
                    distance_bins["<1kb"] += 1
                elif dist < 10_000:
                    distance_bins["1-10kb"] += 1
                elif dist < 100_000:
                    distance_bins["10-100kb"] += 1
                elif dist < 1_000_000:
                    distance_bins["100kb-1Mb"] += 1
                else:
                    distance_bins[">=1Mb"] += 1
            else:
                trans_pairs += 1

            pos51 = int(fields[idx["pos51"]])
            pos52 = int(fields[idx["pos52"]])
            pos31 = int(fields[idx["pos31"]])
            pos32 = int(fields[idx["pos32"]])

            fragment_lengths[abs(pos31 - pos51) + 1] += 1
            fragment_lengths[abs(pos32 - pos52) + 1] += 1

            seq1 = fields[idx["seq1"]]
            seq2 = fields[idx["seq2"]]
            if seq1 and seq1 != ".":
                read_lengths[len(seq1)] += 1
            if seq2 and seq2 != ".":
                read_lengths[len(seq2)] += 1

    cis_trans_ratio = cis_pairs / trans_pairs if trans_pairs > 0 else "NA"
    pct_cis_long_among_cis = 100 * cis_long_pairs / cis_pairs if cis_pairs > 0 else 0
    pct_cis_long_among_total = 100 * cis_long_pairs / total_pairs if total_pairs > 0 else 0
    pct_cis = 100 * cis_pairs / total_pairs if total_pairs > 0 else 0
    pct_trans = 100 * trans_pairs / total_pairs if total_pairs > 0 else 0

    modal_read_length = read_lengths.most_common(1)[0][0] if read_lengths else "NA"

    # Distance distribution plot
    ordered_bins = ["<1kb", "1-10kb", "10-100kb", "100kb-1Mb", ">=1Mb"]
    plt.figure(figsize=(7, 4))
    plt.bar(ordered_bins, [distance_bins[b] for b in ordered_bins])
    plt.ylabel("Number of cis pairs")
    plt.xlabel("Cis genomic distance")
    plt.title("Cis distance distribution")
    distance_plot = fig_to_base64()

    # Fragment length distribution plot
    frag_items = sorted(fragment_lengths.items())
    frag_x = [x for x, y in frag_items if x <= 1000]
    frag_y = [y for x, y in frag_items if x <= 1000]
    plt.figure(figsize=(7, 4))
    plt.plot(frag_x, frag_y)
    plt.ylabel("Count")
    plt.xlabel("Fragment length")
    plt.title("Fragment length distribution (<=1000 bp)")
    frag_plot = fig_to_base64()

    def fmt(x):
        if isinstance(x, float):
            return f"{x:.2f}"
        return str(x)

    html = f"""<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>Pair-level QC</title>
<style>
body {{ font-family: Arial, sans-serif; margin: 40px; }}
table {{ border-collapse: collapse; margin-bottom: 30px; }}
th, td {{ border: 1px solid #ddd; padding: 8px 12px; }}
th {{ background: #f2f2f2; }}
img {{ max-width: 800px; display: block; margin-bottom: 30px; }}
</style>
</head>
<body>
<h1>Pair-level QC Report</h1>
<h2>Sample: {sample_id}</h2>

<h2>Summary metrics</h2>
<table>
<tr><th>Metric</th><th>Value</th></tr>
<tr><td>Total pairs / non-duplicate pairs</td><td>{total_pairs:,}</td></tr>
<tr><td>Cis pairs</td><td>{cis_pairs:,} ({fmt(pct_cis)}%)</td></tr>
<tr><td>Trans pairs</td><td>{trans_pairs:,} ({fmt(pct_trans)}%)</td></tr>
<tr><td>Cis/trans ratio</td><td>{fmt(cis_trans_ratio)}</td></tr>
<tr><td>Cis pairs &gt;= {args.cis_distance:,} bp</td><td>{cis_long_pairs:,}</td></tr>
<tr><td>% cis pairs &gt;= {args.cis_distance:,} bp among cis pairs</td><td>{fmt(pct_cis_long_among_cis)}%</td></tr>
<tr><td>% cis pairs &gt;= {args.cis_distance:,} bp among total pairs</td><td>{fmt(pct_cis_long_among_total)}%</td></tr>
<tr><td>Modal read length</td><td>{modal_read_length}</td></tr>
</table>

<h2>Cis distance distribution</h2>
<img src="data:image/png;base64,{distance_plot}">

<h2>Fragment length distribution</h2>
<img src="data:image/png;base64,{frag_plot}">

<h2>Distance bin counts</h2>
<table>
<tr><th>Distance bin</th><th>Cis pairs</th></tr>
{''.join(f"<tr><td>{b}</td><td>{distance_bins[b]:,}</td></tr>" for b in ordered_bins)}
</table>

</body>
</html>
"""

    args.out.write_text(html)


if __name__ == "__main__":
    main()
