import pathlib
import unittest


ROOT = pathlib.Path(__file__).resolve().parents[1]


class YNeighborIndexTest(unittest.TestCase):
    def test_both_solvers_use_explicit_y_offsets(self):
        for name in ("sle_L31_C035.c", "sle_L31_C070.c"):
            source = (ROOT / "solvers" / name).read_text(encoding="utf-8")
            for expression in (
                "alpha.y[0,1]",
                "gamma.y[0,1]",
                "pres[0,-1]",
                "f[0,-1]",
                "f_old[0,-1]",
            ):
                with self.subTest(solver=name, expression=expression):
                    self.assertIn(expression, source)
            self.assertNotIn("alpha.y[1]", source)
            self.assertNotIn("gamma.y[1]", source)


if __name__ == "__main__":
    unittest.main()
