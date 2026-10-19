import pathlib
import unittest


ROOT = pathlib.Path(__file__).resolve().parents[1]
V3 = ROOT.parent / "film_rupture" / "v3_tiered"


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

    def test_v3_tiered_sources_keep_fixes(self):
        sources = sorted(V3.glob("source_v3_*.c"))
        self.assertGreaterEqual(len(sources), 8)
        for path in sources:
            source = path.read_text(encoding="utf-8")
            for expression in (
                "alpha.y[0,1]",
                "gamma.y[0,1]",
                "pres[0,-1]",
                "f[0,-1]",
                "f_old[0,-1]",
                "(ix%nidx+nidx)%nidx",          # MPI periodic noise wrap
                "hstar = 0.0;",                 # no-repulsion physics (77ce5cb)
                "hmin < 0.3 ? 1e-2 : 1e-1",     # tiered dt cap
                "sect_hmin_%06.2f_r%d.dat",     # per-0.01 section output
                "sect_t_%09.0f_r%d.dat",        # periodic section output
            ):
                with self.subTest(solver=path.name, expression=expression):
                    self.assertIn(expression, source)
            self.assertNotIn("alpha.y[1]", source)
            self.assertNotIn("gamma.y[1]", source)


if __name__ == "__main__":
    unittest.main()
