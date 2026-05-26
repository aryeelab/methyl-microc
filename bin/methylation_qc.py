#!/usr/bin/env python3

import argparse
import base64
import gzip
import io
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def fig_to_base64(fig):
    buf = io.BytesIO()
    fig.tight_layout()
    fig.savefig(buf, format="png", dpi=150)
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("utf-8")


def accumulate_meth_string(meth, frag_len, sums, counts):
    if len(meth) != frag_len:
        return False

    for i, ch in enumerate(meth):
        if ch == "1":
            sums[i] += 1
            counts[i] += 1
        elif ch == "0":
            counts[i] += 1

    return True


def infer_sample_id(path):
    name = path.name
    for suffix in [".meth.pairs.gz", ".meth.pairs"]:
        if name.endswith(suffix):
            return name[:-len(suffix)]
    return path.stem


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pairs", required=True, type=Path)
    parser.add_argument("--out", default="methylation_level_qc.html", type=Path)
    parser.add_argument("--frag-len", default=150, type=int)
    parser.add_argument("--max-data-lines", default=10_000_000, type=int)
    args = parser.parse_args()

    sample_id = infer_sample_id(args.pairs)
    frag_len = args.frag_len

    sums1 = np.zeros(frag_len, dtype=np.int64)
    cnts1 = np.zeros(frag_len, dtype=np.int64)
    sums2 = np.zeros(frag_len, dtype=np.int64)
    cnts2 = np.zeros(frag_len, dtype=np.int64)

    n_frag1 = 0
    n_frag2 = 0
    n_data_lines = 0

    cols = None
    idx = None

    with open_text(args.pairs) as fh:
        for line in fh:
            if line.startswith("#columns:"):
                cols = line.strip().split(":", 1)[1].split()
                idx = {c: i for i, c in enumerate(cols)}
                if "meth1" not in idx or "meth2" not in idx:
                    raise ValueError("Expected meth1/meth2 columns in #columns header")
                continue

            if line.startswith("#") or not line.strip():
                continue

            if idx is None:
                raise ValueError("Missing #columns header in methylated pairs file")

            n_data_lines += 1
            if n_data_lines > args.max_data_lines:
                break

            parts = line.rstrip("\n").split("\t")

            if accumulate_meth_string(parts[idx["meth1"]], frag_len, sums1, cnts1):
                n_frag1 += 1
            if accumulate_meth_string(parts[idx["meth2"]], frag_len, sums2, cnts2):
                n_frag2 += 1

    avg1 = sums1 / np.where(cnts1 == 0, np.nan, cnts1)
    avg2 = sums2 / np.where(cnts2 == 0, np.nan, cnts2)
    avg_all = (sums1 + sums2) / np.where((cnts1 + cnts2) == 0, np.nan, (cnts1 + cnts2))

    mean_meth1 = np.nanmean(avg1)
    mean_meth2 = np.nanmean(avg2)
    mean_meth_all = np.nanmean(avg_all)

    pos = np.arange(1, frag_len + 1)

    fig, ax = plt.subplots(figsize=(10, 3))
    ax.plot(pos, avg1, alpha=0.7, label=f"Fragment 1 (n={n_frag1})")
    ax.plot(pos, avg2, alpha=0.7, label=f"Fragment 2 (n={n_frag2})")
    ax.set_xlabel("Position along fragment (bp)")
    ax.set_xticks(np.arange(10, frag_len + 1, 10))
    ax.set_ylabel("Average methylation")
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, alpha=0.2)
    ax.legend(loc="lower right", fontsize=9)

    ax2 = ax.twinx()
    ax2.plot(pos, cnts1 + cnts2, color="gray", alpha=0.25, lw=1)
    ax2.set_ylabel("CpG calls contributing")

    ax.set_title(f"Positional methylation bias (fragment length = {frag_len} bp)")
    methylation_plot = fig_to_base64(fig)

    html = f"""<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>Methylation-level QC</title>
<style>
body {{ font-family: Arial, sans-serif; margin: 40px; }}
table {{ border-collapse: collapse; margin-bottom: 30px; }}
th, td {{ border: 1px solid #ddd; padding: 8px 12px; }}
th {{ background: #f2f2f2; }}
img {{ max-width: 1000px; display: block; margin-bottom: 30px; }}
</style>
</head>
<body>
<h1>Methylation-level QC Report</h1>
<h2>Sample: {sample_id}</h2>

<h2>Summary metrics</h2>
<table>
<tr><th>Metric</th><th>Value</th></tr>
<tr><td>Fragment length used</td><td>{frag_len} bp</td></tr>
<tr><td>Maximum data lines scanned</td><td>{args.max_data_lines:,}</td></tr>
<tr><td>Data lines actually scanned</td><td>{min(n_data_lines, args.max_data_lines):,}</td></tr>
<tr><td>Fragment 1 used</td><td>{n_frag1:,}</td></tr>
<tr><td>Fragment 2 used</td><td>{n_frag2:,}</td></tr>
<tr><td>Mean methylation, Fragment 1</td><td>{mean_meth1:.4f}</td></tr>
<tr><td>Mean methylation, Fragment 2</td><td>{mean_meth2:.4f}</td></tr>
<tr><td>Mean methylation, combined</td><td>{mean_meth_all:.4f}</td></tr>
</table>

<h2>Positional methylation bias</h2>
<img src="data:image/png;base64,{methylation_plot}">

</body>
</html>
"""

    args.out.write_text(html)


if __name__ == "__main__":
    main()
