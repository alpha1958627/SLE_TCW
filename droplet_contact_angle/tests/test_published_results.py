"""Verify the released data against raw profiles, metadata, and solver identity."""
from pathlib import Path
import math
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'analysis'))
from published_results import load_measurements, verify_assets


class PublishedResults(unittest.TestCase):
    def test_raw_profiles_reproduce_six_points_and_metadata(self):
        rows = load_measurements()
        self.assertEqual(len(rows), 6)
        self.assertEqual({(r['nominal_angle_deg'], r['level']) for r in rows},
                         {(angle, level) for angle in (15, 25, 35) for level in (13, 14)})
        for row in rows:
            k = row['k']
            # Integral of Pi from hstar to infinity gives C=1/4 for n=5,m=3.
            self.assertAlmostEqual(math.degrees(math.acos(1-k*.01/4)),
                                   row['nominal_angle_deg'], places=6)
        passed = [r['case'] for r in rows if r['radius_shape_gate_passed']]
        self.assertEqual(passed, ['k35_L13'])

    def test_assets_and_exact_tested_solver_hash(self):
        self.assertGreater(verify_assets(), 30)


if __name__ == '__main__':
    unittest.main()
