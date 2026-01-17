#!/usr/bin/env python3
"""Sanity-check meth1/meth2 strings against reference CpG sites for one .pairsam record.

This script is intended for quick validation during tutorial/README usage.
It extracts either the Nth non-header record (default: first) or the first
record matching a given readID, fetches the reference
sequence spanning pos5..pos3 (inclusive) for each side, orients it 5'->3'
strand-aware, and checks:
  1) len(meth) == abs(pos3-pos5)+1
  2) CpG cytosines in the oriented reference align to meth calls

Requires:
  - pairsam columns: chrom{1,2},pos5{1,2},pos3{1,2},strand{1,2},cigar{1,2},seq{1,2},meth{1,2}
  - samtools in PATH (used for faidx access; will create FASTA .fai if missing)
"""

from __future__ import annotations

import argparse
import gzip
import os
import subprocess
import sys
from typing import Iterable


def _open_maybe_gzip(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def revcomp(seq: str) -> str:
    tab = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(tab)[::-1]


def ensure_fai(fasta: str, samtools: str) -> None:
    fai = fasta + ".fai"
    if os.path.exists(fai):
        return
    subprocess.run([samtools, "faidx", fasta], check=True)


def fetch_fasta(samtools: str, fasta: str, chrom: str, start_1based: int, end_1based: int) -> str:
    region = f"{chrom}:{start_1based}-{end_1based}"
    out = subprocess.check_output([samtools, "faidx", fasta, region], text=True)
    lines = out.splitlines()
    return "".join(lines[1:]).upper()


def _select_data_record(
    pairs_path: str,
    *,
    record_index_1based: int = 1,
    read_id: str | None = None,
) -> tuple[list[str], list[str]]:
    if record_index_1based < 1:
        raise SystemExit("ERROR: --record must be >= 1")

    cols: list[str] | None = None
    readid_idx: int | None = None
    k = 0

    with _open_maybe_gzip(pairs_path) as fh:
        for line in fh:
            if line.startswith("#columns:"):
                cols = line.rstrip("\n").split("\t")[1:]
                if read_id is not None:
                    if "readID" not in cols:
                        raise SystemExit("ERROR: --readID requires a 'readID' column")
                    readid_idx = cols.index("readID")
                continue

            if not line.strip() or line.startswith("#"):
                continue
            if cols is None:
                raise SystemExit("ERROR: no #columns header found before first data record")

            parts = line.rstrip("\n").split("\t")
            k += 1
            if read_id is not None:
                if readid_idx is not None and parts[readid_idx] == read_id:
                    return cols, parts
            else:
                if k == record_index_1based:
                    return cols, parts

    if cols is None:
        raise SystemExit("ERROR: no #columns header found")
    if read_id is not None:
        raise SystemExit(f"ERROR: no data record found with readID={read_id!r}")
    raise SystemExit(f"ERROR: no data record at index {record_index_1based}")


def main(argv: Iterable[str]) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", default="results/pairs/test_sample.meth.pairsam.gz", help="Input .pairsam(.gz) with meth1/meth2")
    ap.add_argument("--fasta", default="references/chr22/chr22.fa", help="Reference FASTA (faidx-indexed)")
    sel = ap.add_mutually_exclusive_group()
    sel.add_argument("--record", type=int, default=1, help="1-based index of the non-header record to validate")
    sel.add_argument("--readID", dest="read_id", default=None, help="Validate the first record whose readID matches exactly")
    ap.add_argument("--samtools", default="samtools", help="samtools executable (default: samtools)")
    args = ap.parse_args(list(argv))

    cols, first = _select_data_record(args.pairs, record_index_1based=args.record, read_id=args.read_id)

    need = [
        "readID",
        "chrom1",
        "pos51",
        "pos31",
        "strand1",
        "cigar1",
        "seq1",
        "meth1",
        "chrom2",
        "pos52",
        "pos32",
        "strand2",
        "cigar2",
        "seq2",
        "meth2",
    ]
    missing = [c for c in need if c not in cols]
    if missing:
        raise SystemExit(f"ERROR: missing required columns: {missing}")
    idx = {c: cols.index(c) for c in need}

    ensure_fai(args.fasta, args.samtools)

    read_id = first[idx["readID"]]
    print(f"readID: {read_id}\n")

    def report_end(end: int) -> bool:
        chrom = first[idx[f"chrom{end}"]]
        pos5 = int(first[idx[f"pos5{end}"]])
        pos3 = int(first[idx[f"pos3{end}"]])
        strand = first[idx[f"strand{end}"]]
        cigar = first[idx[f"cigar{end}"]]
        seq = first[idx[f"seq{end}"]]
        meth = first[idx[f"meth{end}"]]

        start = min(pos5, pos3)
        endp = max(pos5, pos3)
        fraglen = abs(pos3 - pos5) + 1

        ref_genomic = fetch_fasta(args.samtools, args.fasta, chrom, start, endp)
        ref = revcomp(ref_genomic) if strand == "-" else ref_genomic

        ok = True
        if len(ref) != fraglen:
            print(f"ERROR end{end}: len(ref)={len(ref)} != fraglen={fraglen}")
            ok = False
        if len(meth) != fraglen:
            print(f"ERROR end{end}: len(meth)={len(meth)} != fraglen={fraglen}")
            ok = False

        cpg_idx = [i for i in range(len(ref) - 1) if ref[i : i + 2] == "CG"]
        cpg_set = set(cpg_idx)
        bad_cpg = [i for i in cpg_idx if meth[i] not in "01."]
        bad_non = [i for i in range(len(ref)) if (i not in cpg_set) and meth[i] in "01"]

        print(f"=== FRAGMENT {end} ===")
        print(f"chrom:   {chrom}")
        print(f"strand:  {strand}")
        print(f"pos5:    {pos5}")
        print(f"pos3:    {pos3}")
        print(f"span:    {chrom}:{start}-{endp} (len={fraglen})")
        print(f"cigar:   {cigar}")
        print(f"seq:     {seq}")
        print(f"ref:     {ref}")
        print(f"meth:    {meth}")

        marks = [" "] * len(ref)
        for i in cpg_idx:
            marks[i] = "^"  # points at the 'C' in each 'CG'
        print("cpg-C:   " + "".join(marks))

        print(f"CpG sites in span: {len(cpg_idx)}")
        if bad_cpg:
            print(f"ERROR end{end}: CpG positions with unexpected meth char (expected 0/1/.): {bad_cpg[:10]}")
            ok = False
        if bad_non:
            print(f"ERROR end{end}: non-CpG positions with 0/1 in meth string: {bad_non[:10]}")
            ok = False
        print("OK" if ok else "FAIL")
        print()
        return ok

    ok1 = report_end(1)
    ok2 = report_end(2)
    if not (ok1 and ok2):
        return 1

    print("All checks passed: meth strings match abs(pos3-pos5)+1 and align to CpG sites in the reference span.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))

