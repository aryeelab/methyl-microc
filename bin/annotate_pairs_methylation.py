#!/usr/bin/env python3

from __future__ import annotations

import argparse
import gzip
import os
import re
import sys
from collections import OrderedDict

_CIGAR_RE = re.compile(r"([0-9]+)([MIDNSHP=X])")
CIGAR_CACHE = {}


def parse_cigar(cigar: str):
    if cigar in CIGAR_CACHE:
        return CIGAR_CACHE[cigar]

    if cigar == "*" or not cigar:
        ops = []
    else:
        ops = [(int(n), op) for n, op in _CIGAR_RE.findall(cigar)]
        if "".join(f"{n}{op}" for n, op in ops) != cigar:
            raise ValueError(f"Invalid CIGAR: {cigar}")

    CIGAR_CACHE[cigar] = ops
    return ops


class FastaFai:
    def __init__(self, fasta_path: str, fai_path: str, cache_window: int = 2_000_000, max_windows: int = 32):
        self.fasta_path = fasta_path
        self.fai_path = fai_path
        self.cache_window = cache_window
        self.max_windows = max_windows
        self._fp = open(fasta_path, "rb")
        self._idx = self._load_fai(fai_path)
        self._cache = OrderedDict()

    def close(self):
        if self._fp and not self._fp.closed:
            self._fp.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()
        return False

    @staticmethod
    def _load_fai(fai_path: str):
        idx = {}
        with open(fai_path, "rt") as f:
            for line in f:
                if not line.strip():
                    continue
                name, ln, offset, line_bases, line_width, *_ = line.rstrip("\n").split("\t")
                idx[name] = {
                    "len": int(ln),
                    "offset": int(offset),
                    "line_bases": int(line_bases),
                    "line_width": int(line_width),
                }
        return idx

    def fetch_raw(self, name: str, start_1based: int, end_1based: int) -> str:
        if start_1based > end_1based:
            return ""
        if name not in self._idx:
            raise KeyError(f"Contig {name!r} not found in FASTA index")

        info = self._idx[name]
        clen = info["len"]
        start = max(1, start_1based)
        end = min(clen, end_1based)
        if start > end:
            return ""

        lb = info["line_bases"]
        lw = info["line_width"]
        offset = info["offset"]

        i0 = start - 1
        start_byte = offset + (i0 // lb) * lw + (i0 % lb)
        n_bases = end - start + 1
        lines_crossed = ((i0 % lb) + n_bases - 1) // lb
        n_bytes = n_bases + lines_crossed * (lw - lb)

        self._fp.seek(start_byte)
        raw = self._fp.read(n_bytes)
        seq = raw.replace(b"\n", b"").replace(b"\r", b"")[:n_bases]
        return seq.decode("ascii").upper()

    def fetch(self, name: str, start_1based: int, end_1based: int) -> str:
        if start_1based > end_1based:
            return ""
        if name not in self._idx:
            return "N" * max(0, end_1based - start_1based + 1)

        clen = self._idx[name]["len"]
        start = max(1, start_1based)
        end = min(clen, end_1based)
        if start > end:
            return ""

        block = (start - 1) // self.cache_window
        block_start = block * self.cache_window + 1
        block_end = min(clen, block_start + self.cache_window - 1)

        key = (name, block_start, block_end)
        if key not in self._cache:
            self._cache[key] = self.fetch_raw(name, block_start, block_end)
            self._cache.move_to_end(key)
            while len(self._cache) > self.max_windows:
                self._cache.popitem(last=False)
        else:
            self._cache.move_to_end(key)

        seq = self._cache[key]
        s = start - block_start
        e = end - block_start + 1
        return seq[s:e]


def build_ref_base_map(ref_start_1based: int, cigar: str, seq_aln: str) -> dict[int, str]:
    seq_aln = (seq_aln or "").upper()
    ref_pos = ref_start_1based
    q_pos = 0
    out = {}

    for length, op in parse_cigar(cigar):
        if op in ("M", "=", "X"):
            sub = seq_aln[q_pos:q_pos + length]
            for i, base in enumerate(sub):
                out[ref_pos + i] = base
            ref_pos += length
            q_pos += length
        elif op in ("I", "S"):
            q_pos += length
        elif op in ("D", "N"):
            ref_pos += length
        elif op in ("H", "P"):
            continue
        else:
            raise ValueError(f"Unsupported CIGAR op: {op!r} in {cigar!r}")

    return out


def methyl_string_for_side(fasta, chrom, strand, pos5, pos3, cigar, seq):
    if chrom == "!" or pos5 <= 0 or pos3 <= 0 or cigar in ("", "*") or not seq:
        return "."

    step = 1 if pos3 >= pos5 else -1
    frag_len = abs(pos3 - pos5) + 1
    left = min(pos5, pos3)
    right = max(pos5, pos3)

    ref_start = left - 1
    ref_end = right + 1
    ref_seq = fasta.fetch(chrom, ref_start, ref_end)

    seq = seq.upper()

    simple_match = cigar in (f"{frag_len}M", f"{frag_len}=", f"{frag_len}X") and len(seq) >= frag_len

    if not simple_match:
        ref_base_map = build_ref_base_map(left, cigar, seq)
    else:
        ref_base_map = None

    chars = []
    p = pos5

    for offset in range(frag_len):
        i = p - ref_start

        if strand == "+":
            b0 = ref_seq[i] if 0 <= i < len(ref_seq) else "N"
            b1 = ref_seq[i + 1] if 0 <= i + 1 < len(ref_seq) else "N"

            if not (b0 == "C" and b1 == "G"):
                chars.append("-")
            else:
                if simple_match:
                    rb = seq[offset]
                else:
                    rb = ref_base_map.get(p, "")
                if rb == "C":
                    chars.append("1")
                elif rb == "T":
                    chars.append("0")
                else:
                    chars.append(".")

        else:
            b_1 = ref_seq[i - 1] if 0 <= i - 1 < len(ref_seq) else "N"
            b0 = ref_seq[i] if 0 <= i < len(ref_seq) else "N"

            if not (b_1 == "C" and b0 == "G"):
                chars.append("-")
            else:
                if simple_match:
                    rb = seq[offset]
                else:
                    rb = ref_base_map.get(p, "")
                if rb == "G":
                    chars.append("1")
                elif rb == "A":
                    chars.append("0")
                else:
                    chars.append(".")

        p += step

    return "".join(chars)


def open_text_maybe_gz(path, mode):
    if path is None or path == "-":
        return sys.stdin if "r" in mode else sys.stdout
    if path.endswith(".gz"):
        return gzip.open(path, mode + "t")
    return open(path, mode + "t")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--fai", default=None)
    ap.add_argument("--input", default="-")
    ap.add_argument("--output", default="-")
    ap.add_argument("--progress-every", type=int, default=1_000_000)
    args = ap.parse_args(argv)

    fai = args.fai or (args.fasta + ".fai")
    if not os.path.exists(args.fasta):
        ap.error(f"FASTA not found: {args.fasta}")
    if not os.path.exists(fai):
        ap.error(f"FASTA index not found: {fai}")

    with FastaFai(args.fasta, fai) as fasta, open_text_maybe_gz(args.input, "r") as inp, open_text_maybe_gz(args.output, "w") as out:
        columns = None
        header_lines = []
        first_data = None

        for line in inp:
            if line.startswith("#"):
                if line.startswith("#columns:"):
                    cols = line.rstrip("\n").split()[1:]
                    if "meth1" not in cols:
                        cols = cols + ["meth1", "meth2"]
                    columns = cols
                    header_lines.append("#columns:\t" + "\t".join(cols) + "\n")
                else:
                    header_lines.append(line)
            else:
                first_data = line
                break

        if columns is None:
            raise RuntimeError("Missing #columns header")

        out.writelines(header_lines)

        idx = {c: i for i, c in enumerate(columns)}
        required = [
            "chrom1", "chrom2", "strand1", "strand2",
            "pos51", "pos31", "pos52", "pos32",
            "cigar1", "cigar2", "seq1", "seq2",
        ]
        missing = [c for c in required if c not in idx]
        if missing:
            raise RuntimeError("Input pairs is missing required columns: " + ",".join(missing))

        n = 0

        def handle_data_line(data_line):
            fields = data_line.rstrip("\n").split("\t")

            chrom1 = fields[idx["chrom1"]]
            strand1 = fields[idx["strand1"]]
            pos51 = int(fields[idx["pos51"]])
            pos31 = int(fields[idx["pos31"]])
            cigar1 = fields[idx["cigar1"]]
            seq1 = fields[idx["seq1"]]

            chrom2 = fields[idx["chrom2"]]
            strand2 = fields[idx["strand2"]]
            pos52 = int(fields[idx["pos52"]])
            pos32 = int(fields[idx["pos32"]])
            cigar2 = fields[idx["cigar2"]]
            seq2 = fields[idx["seq2"]]

            m1 = methyl_string_for_side(fasta, chrom1, strand1, pos51, pos31, cigar1, seq1)
            m2 = methyl_string_for_side(fasta, chrom2, strand2, pos52, pos32, cigar2, seq2)

            if "meth1" in idx and idx["meth1"] < len(fields):
                fields[idx["meth1"]] = m1
                fields[idx["meth2"]] = m2
                return "\t".join(fields) + "\n"

            return "\t".join(fields + [m1, m2]) + "\n"

        if first_data is not None and first_data.strip():
            out.write(handle_data_line(first_data))
            n += 1

        for line in inp:
            if not line.strip() or line.startswith("#"):
                continue

            out.write(handle_data_line(line))
            n += 1

            if args.progress_every > 0 and n % args.progress_every == 0:
                print(f"[annotate_pairs_methylation] processed {n:,} pairs", file=sys.stderr, flush=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
