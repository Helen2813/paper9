import csv
import tempfile
import unittest
from pathlib import Path
from p9darc.core import (read_pseudobulk, read_gene_z, donor_paired_profiles,
                         select_markers, diagnostic_ranks, make_report, leave_one_donor_out, donor_balance)


class DirectionAuditTests(unittest.TestCase):
    def make_rows(self, n=8):
        rows = []
        for d in range(n):
            for i in range(30):
                g = "G{}".format(i)
                rows.append(dict(gene=g, donor=str(d), cell_type="stroma", cpm=100 if i < 15 else 1))
                rows.append(dict(gene=g, donor=str(d), cell_type="vascular", cpm=1 if i < 15 else 100))
        return rows

    def test_donor_paired_directions(self):
        prof = donor_paired_profiles(self.make_rows(), 5)
        self.assertGreater(prof["stroma"]["G1"]["delta"], 0)
        self.assertLess(prof["stroma"]["G20"]["delta"], 0)
        self.assertEqual(prof["stroma"]["G1"]["n_donors"], 8)
        markers = select_markers(prof, 10)
        self.assertIn("G1", markers["stroma"]["UP"])
        self.assertIn("G20", markers["stroma"]["DOWN"])

    def test_insufficient_donors_withholds_profiles(self):
        self.assertEqual(donor_paired_profiles(self.make_rows(3), 5), {})

    def test_no_calibrated_claim_from_toy(self):
        rows = self.make_rows()
        z = {"G{}".format(i): 3.0 if i < 15 else 0.3 for i in range(30)}
        report = make_report(rows, z, {"stroma": "stromal", "vascular": "vascular"}, 10, 5, True)
        self.assertEqual(report["inferential_status"], "NOT_CALIBRATED")
        self.assertFalse(report["has_formal_genetic_inference"])
        self.assertEqual(report["rankings"]["UP"][0]["cell_type"], "stroma")
        self.assertEqual(sum(report["lodo"]["top_cell_counts"]["UP"].values()), 8)
        self.assertNotIn("pvalue", report)

    def test_donor_balance_dominance_detected(self):
        rows = [dict(gene="G1", donor="D1", cell_type="stroma", cpm=10, n_cells=391),
                dict(gene="G1", donor="D2", cell_type="stroma", cpm=10, n_cells=7),
                dict(gene="G1", donor="D3", cell_type="stroma", cpm=10, n_cells=2)]
        b = donor_balance(rows)["stroma"]
        self.assertGreater(b["largest_donor_share"], 0.97)
        self.assertLess(b["effective_donor_count_from_cell_weights"], 1.1)
        self.assertTrue(b["donor_imbalance_warning"])

    def test_duplicate_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/"pb.tsv"
            with p.open("w", newline="") as f:
                w = csv.writer(f, delimiter="\t")
                w.writerow(["gene", "donor", "cell_type", "cpm"])
                w.writerow(["G1", "D1", "A", "3"])
                w.writerow(["G1", "D1", "A", "3"])
            with self.assertRaises(ValueError):
                read_pseudobulk(p)

    def test_negative_cpm_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/"pb.tsv"
            p.write_text("gene\tdonor\tcell_type\tcpm\nG1\tD1\tA\t-3\n")
            with self.assertRaises(ValueError):
                read_pseudobulk(p)

    def test_nonoverlap_no_genetic_ranks(self):
        rows = self.make_rows()
        prof = donor_paired_profiles(rows, 5)
        ranks = diagnostic_ranks(select_markers(prof, 10), {"NONMATCH": 3})
        self.assertTrue(all(not val for val in ranks.values()))


if __name__ == "__main__":
    unittest.main()
