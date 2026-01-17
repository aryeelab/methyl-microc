import importlib.util
import os
import pathlib
import tempfile
import unittest


def _load_module():
    repo_root = pathlib.Path(__file__).resolve().parents[1]
    mod_path = repo_root / "bin" / "annotate_pairsam_methylation.py"
    spec = importlib.util.spec_from_file_location("annotate_pairsam_methylation", mod_path)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


M = _load_module()


def _write_unwrapped_fasta_with_fai(tmpdir: str, name: str, seq: str):
    fasta_path = os.path.join(tmpdir, "ref.fa")
    fai_path = fasta_path + ".fai"
    seq = seq.upper()
    with open(fasta_path, "wb") as f:
        header = f">{name}\n".encode("ascii")
        f.write(header)
        offset = f.tell()
        f.write(seq.encode("ascii") + b"\n")

    # Unwrapped FASTA: all bases on one line.
    line_bases = len(seq)
    line_width = line_bases + 1
    with open(fai_path, "wt") as f:
        f.write(f"{name}\t{len(seq)}\t{offset}\t{line_bases}\t{line_width}\n")
    return fasta_path, fai_path


class TestAnnotatePairsamMethylation(unittest.TestCase):
    def test_fasta_fai_fetch_and_base(self):
        with tempfile.TemporaryDirectory() as d:
            fasta, fai = _write_unwrapped_fasta_with_fai(d, "chr1", "ACGCGT")
            with M.FastaFai(fasta, fai) as ff:
                self.assertEqual(ff.fetch("chr1", 1, 6), "ACGCGT")
                self.assertEqual(ff.base("chr1", 2), "C")
                self.assertEqual(ff.base("chr1", 999), "N")

    def test_methyl_string_plus_fragment(self):
        # chr1: A C G C G T
        # fragment is pos5=2 -> pos3=5 on '+' strand (positions 2..5)
        # CpGs at pos2 (C+next G) and pos4 (C+next G)
        # aligned bases at pos2=C (meth), pos4=T (unmeth)
        with tempfile.TemporaryDirectory() as d:
            fasta, fai = _write_unwrapped_fasta_with_fai(d, "chr1", "ACGCGT")
            with M.FastaFai(fasta, fai) as ff:
                s = M.methyl_string_for_side(
                    ff, "chr1", "+", 2, 5, "4M", "CGTG"  # ref-pos 2..5
                )
                self.assertEqual(s, "1-0-")

    def test_methyl_string_minus_fragment(self):
        # Same reference, but fragment is pos5=5 -> pos3=2 on '-' strand.
        # CpG cytosines on minus correspond to reference 'G' preceded by 'C': positions 5 and 3.
        # Along 5'->3' (5,4,3,2): pos5 base A -> unmeth (0); pos3 base G -> meth (1)
        with tempfile.TemporaryDirectory() as d:
            fasta, fai = _write_unwrapped_fasta_with_fai(d, "chr1", "ACGCGT")
            with M.FastaFai(fasta, fai) as ff:
                s = M.methyl_string_for_side(
                    ff,
                    "chr1",
                    "-",
                    5,
                    2,
                    "4M",
                    "CGCA",  # pairtools-style: already oriented left->right on reference
                )
                self.assertEqual(s, "0-1-")

    def test_methyl_string_minus_fragment_softclip_requires_cigar_reverse(self):
        # Same as test_methyl_string_minus_fragment, but the read has a leading
        # soft-clip. In pairtools-style pairsam, the seq/cigar are already
        # oriented left->right on the reference, so the soft-clip for a '-'
        # alignment appears on the opposite end compared to the original read.
        with tempfile.TemporaryDirectory() as d:
            fasta, fai = _write_unwrapped_fasta_with_fai(d, "chr1", "ACGCGT")
            with M.FastaFai(fasta, fai) as ff:
                s = M.methyl_string_for_side(
                    ff,
                    "chr1",
                    "-",
                    5,
                    2,
                    "4M1S",
                    "CGCAT",  # pairtools-style: seq and cigar already reversed for '-'
                )
                self.assertEqual(s, "0-1-")

    def test_unclear_call_deletion(self):
        # chr1: A C G
        # CpG at pos2, but read has a deletion at ref pos2 => '.'
        with tempfile.TemporaryDirectory() as d:
            fasta, fai = _write_unwrapped_fasta_with_fai(d, "chr1", "ACG")
            with M.FastaFai(fasta, fai) as ff:
                s = M.methyl_string_for_side(ff, "chr1", "+", 2, 3, "1D1M", "G")
                self.assertEqual(s, ".-")


if __name__ == "__main__":
    unittest.main()
