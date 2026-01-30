#!/usr/bin/env python3
"""Annotate a .pairs file with per-fragment CpG methylation strings.

The methylation string is defined over the *fragment* for each side, i.e. the
reference interval from pos5->pos3 (5' to 3', strand-aware), NOT over the full
read sequence.

Output columns appended:
  - meth1 : methylation string for side 1
  - meth2 : methylation string for side 2

Alphabet:
  - '-' : position is not a CpG cytosine on the fragment strand
  - '1' : CpG, methylated (C on + strand or G on - strand)
  - '0' : CpG, unmethylated (T on + strand or A on - strand)
  - '.' : CpG, but no clear call (deletion, N, other base, missing data)

Requires pairs to contain: pos51,pos31,pos52,pos32,strand1,strand2,cigar1,
cigar2,seq1,seq2.
"""

from __future__ import annotations

import argparse
import gzip
import os
import re
import sys
from collections import OrderedDict


_CIGAR_RE = re.compile(r"([0-9]+)([MIDNSHP=X])")


def revcomp(seq: str) -> str:
    tab = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(tab)[::-1]


def parse_cigar(cigar: str):
    if cigar == "*" or not cigar:
        return []
    parts = _CIGAR_RE.findall(cigar)
    if not parts or "".join(n + op for n, op in parts) != cigar:
        raise ValueError(f"Invalid CIGAR: {cigar!r}")
    return [(int(n), op) for n, op in parts]


def reverse_cigar(cigar: str) -> str:
    """Reverse CIGAR operations for use with a reverse-complemented query.

    Many aligners keep SEQ in the original read orientation even for reverse-
    strand alignments, while CIGAR is still reported in reference left->right
    order. If we reverse-complement the sequence to orient it left->right, we
    must also reverse the order of CIGAR operations so that leading/trailing
    clips/indels apply to the correct end of the query.
    """

    ops = parse_cigar(cigar)
    if not ops:
        return cigar
    return "".join(f"{n}{op}" for n, op in reversed(ops))


def build_ref_base_map(ref_start_1based: int, cigar: str, seq_aln: str) -> dict[int, str]:
    """Return map: reference_pos(1-based) -> aligned base (uppercase).

    Only positions that consume both reference+query (M/=X) are included.
    Deletions/skips will be absent from the map.
    """
    seq_aln = (seq_aln or "").upper()
    ref_pos = ref_start_1based
    q_pos = 0
    out: dict[int, str] = {}
    for length, op in parse_cigar(cigar):
        if op in ("M", "=", "X"):
            for i in range(length):
                if q_pos + i < len(seq_aln):
                    out[ref_pos + i] = seq_aln[q_pos + i]
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


class FastaFai:
    """Minimal FASTA random access via .fai (no external deps)."""

    def __init__(self, fasta_path: str, fai_path: str, cache_contigs: int = 8):
        self.fasta_path = fasta_path
        self.fai_path = fai_path
        self._fp = open(fasta_path, "rb")
        self._idx = self._load_fai(fai_path)
        self._cache = OrderedDict()  # contig -> (start,end,seq)
        self._cache_contigs = cache_contigs

    def close(self):
        fp = getattr(self, "_fp", None)
        if fp is not None and not fp.closed:
            fp.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()
        return False

    def __del__(self):
        # Best-effort cleanup; avoids ResourceWarnings in some contexts.
        try:
            self.close()
        except Exception:
            pass

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

    def contig_len(self, name: str) -> int:
        return self._idx[name]["len"]

    def fetch(self, name: str, start_1based: int, end_1based: int) -> str:
        """Fetch [start,end] inclusive, 1-based. Returns uppercase sequence."""
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

        # byte offset of start base (0-based within contig)
        i0 = start - 1
        start_byte = offset + (i0 // lb) * lw + (i0 % lb)
        n_bases = end - start + 1
        # number of newline bytes inside the span
        lines_crossed = ((i0 % lb) + n_bases - 1) // lb
        n_bytes = n_bases + lines_crossed * (lw - lb)
        self._fp.seek(start_byte)
        raw = self._fp.read(n_bytes)
        seq = raw.replace(b"\n", b"").replace(b"\r", b"")[:n_bases]
        return seq.decode("ascii").upper()

    def _fetch_cached(self, name: str, start_1based: int, end_1based: int, pad: int = 2000) -> tuple[int, str]:
        """Return (window_start, window_seq) for a cached window covering [start,end]."""
        info = self._idx.get(name)
        if not info:
            raise KeyError(f"Contig {name!r} not found in FASTA index")
        clen = info["len"]
        start = max(1, start_1based)
        end = min(clen, end_1based)
        if start > end:
            return 1, ""

        rec = self._cache.get(name)
        if rec is not None:
            w0, w1, wseq = rec
            if start >= w0 and end <= w1:
                self._cache.move_to_end(name)
                return w0, wseq

        w0 = max(1, start - pad)
        w1 = min(clen, end + pad)
        wseq = self.fetch(name, w0, w1)
        self._cache[name] = (w0, w1, wseq)
        self._cache.move_to_end(name)
        while len(self._cache) > self._cache_contigs:
            self._cache.popitem(last=False)
        return w0, wseq

    def base(self, name: str, pos_1based: int) -> str:
        if name not in self._idx:
            return "N"
        clen = self._idx[name]["len"]
        if pos_1based < 1 or pos_1based > clen:
            return "N"
        w0, wseq = self._fetch_cached(name, pos_1based, pos_1based)
        return wseq[pos_1based - w0] if wseq else "N"


def methyl_string_for_side(
    fasta: FastaFai,
    chrom: str,
    strand: str,
    pos5: int,
    pos3: int,
    cigar: str,
    seq: str,
) -> str:
    if chrom == "!" or pos5 <= 0 or pos3 <= 0 or not cigar or cigar == "*" or not seq or seq == "*":
        return "."

    step = 1 if pos3 >= pos5 else -1
    frag_len = abs(pos3 - pos5) + 1

    left = min(pos5, pos3)
    # pairtools `parse --add-columns ... cigar,seq` reports `seq{1,2}` already
    # oriented in reference left->right order (for '-' alignments this is the
    # reverse-complement of the original read sequence). The corresponding
    # `cigar{1,2}` is reported in the same orientation.
    #
    # Therefore we can project `seq` onto the reference span using `cigar`
    # directly, independent of strand.
    ref_base_map = build_ref_base_map(left, cigar, seq)

    out = []
    for i in range(frag_len):
        p = pos5 + i * step
        if strand == "+":
            if fasta.base(chrom, p) != "C" or fasta.base(chrom, p + 1) != "G":
                out.append("-")
                continue
            b = ref_base_map.get(p)
            if b == "C":
                out.append("1")
            elif b == "T":
                out.append("0")
            else:
                out.append(".")
        else:
            if fasta.base(chrom, p) != "G" or fasta.base(chrom, p - 1) != "C":
                out.append("-")
                continue
            b = ref_base_map.get(p)
            if b == "G":
                out.append("1")
            elif b == "A":
                out.append("0")
            else:
                out.append(".")
    return "".join(out)


def _open_text_maybe_gz(path: str | None, mode: str):
    if path is None or path == "-":
        return sys.stdin if "r" in mode else sys.stdout
    if path.endswith(".gz"):
        return gzip.open(path, mode + "t")
    return open(path, mode + "t")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="Append meth1/meth2 methylation strings to a pairs file (streaming)."
    )
    ap.add_argument("--fasta", required=True, help="Reference FASTA (unconverted), must be indexed")
    ap.add_argument(
        "--fai",
        default=None,
        help="FASTA index (.fai). Default: --fasta + .fai",
    )
    ap.add_argument("--input", default="-", help="Input .pairs(.gz); default stdin")
    ap.add_argument("--output", default="-", help="Output .pairs(.gz); default stdout")
    args = ap.parse_args(argv)

    fai = args.fai or (args.fasta + ".fai")
    if not os.path.exists(args.fasta):
        ap.error(f"FASTA not found: {args.fasta}")
    if not os.path.exists(fai):
        ap.error(f"FASTA index not found: {fai} (run: samtools faidx {args.fasta})")

    with FastaFai(args.fasta, fai) as fasta, _open_text_maybe_gz(args.input, "r") as inp, _open_text_maybe_gz(
        args.output, "w"
    ) as out:
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
            raise RuntimeError("Missing #columns: header; not a valid .pairs file")

        out.writelines(header_lines)

        idx = {c: i for i, c in enumerate(columns)}
        required = [
            "chrom1",
            "chrom2",
            "strand1",
            "strand2",
            "pos51",
            "pos31",
            "pos52",
            "pos32",
            "cigar1",
            "cigar2",
            "seq1",
            "seq2",
        ]
        missing = [c for c in required if c not in idx]
        if missing:
            raise RuntimeError(
                "Input pairs is missing required columns: "
                + ",".join(missing)
                + "\nHint: generate with pairtools parse --add-columns pos5,pos3,cigar,seq (optionally add --drop-sam to remove sam1/sam2)"
            )

        def handle_data_line(data_line: str):
            fields = data_line.rstrip("\n").split("\t")
            # If the header already had meth1/meth2, preserve existing columns by not duplicating.
            if len(fields) != len(columns) and len(fields) + 2 == len(columns):
                # Common case: header updated by us but body hasn't been yet.
                pass
            elif len(fields) != len(columns) - 2 and len(fields) != len(columns):
                raise RuntimeError(
                    f"Column count mismatch: header has {len(columns)} columns, row has {len(fields)} fields"
                )

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

            # If input already had meth columns, overwrite them; otherwise append.
            if "meth1" in idx and idx["meth1"] < len(fields):
                fields[idx["meth1"]] = m1
                fields[idx["meth2"]] = m2
                return "\t".join(fields) + "\n"
            return "\t".join(fields + [m1, m2]) + "\n"

        if first_data is not None and first_data.strip():
            out.write(handle_data_line(first_data))
        for line in inp:
            if not line.strip() or line.startswith("#"):
                # pairs bodies should not contain '#', but be tolerant.
                continue
            out.write(handle_data_line(line))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
