import importlib.util
import pathlib
import unittest


def _load_module():
    repo_root = pathlib.Path(__file__).resolve().parents[1]
    mod_path = repo_root / "bin" / "validate_pairsam_methylation.py"
    spec = importlib.util.spec_from_file_location("validate_pairsam_methylation", mod_path)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


V = _load_module()


class TestValidatePairsamMethylation(unittest.TestCase):
    def test_orient_for_display_minus_matches_revcomp_reference(self):
        # If seq/cigar are already reference-oriented (pairtools-style), then for
        # '-' strand display we must flip at the *end* so seq(aln) matches the
        # reverse-complemented reference.
        ref_lr = "ACGTA"
        seq = "ACGTA"
        cigar = "5M"
        ref_disp, seq_disp = V.orient_for_display(ref_lr, seq, cigar, "-")
        self.assertEqual(ref_disp, "TACGT")
        self.assertEqual(seq_disp, "TACGT")

    def test_orient_for_display_plus_is_identity(self):
        ref_lr = "ACGTA"
        seq = "ACGTA"
        cigar = "5M"
        ref_disp, seq_disp = V.orient_for_display(ref_lr, seq, cigar, "+")
        self.assertEqual(ref_disp, "ACGTA")
        self.assertEqual(seq_disp, "ACGTA")

    def test_count_matches_to_ref_or_bs(self):
        # ref_disp:    A C G C
        # ref(bs):     A T G T
        # seq_disp:    A T G C  => matches all positions to either ref or ref(bs)
        ref_disp = "ACGC"
        seq_disp = "ATGC"
        n_match, n_bases = V.count_matches_to_ref_or_bs(seq_disp, ref_disp)
        self.assertEqual((n_match, n_bases), (4, 4))

        # gaps and ambiguous bases are ignored in n_bases
        ref_disp = "ACGC"
        seq_disp = "A-NC"
        n_match, n_bases = V.count_matches_to_ref_or_bs(seq_disp, ref_disp)
        # positions: A (match), '-' (skip), N (skip), C (match)
        self.assertEqual((n_match, n_bases), (2, 2))


if __name__ == "__main__":
    unittest.main()

