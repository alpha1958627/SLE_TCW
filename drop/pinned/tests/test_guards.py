"""Regression checks for the actual stencil and headerless-restart failure."""
from pathlib import Path
import tempfile
import unittest
from audit import rows
from check_source import check
from native_dump import self_check

WORK = Path(__file__).parents[1]/'runs'
WORK.mkdir(exist_ok=True)

class Guards(unittest.TestCase):
    def test_real_source_and_negative_control(self):
        source = Path(__file__).parents[1]/'contact_angle.c'
        check(source)
        with tempfile.TemporaryDirectory(dir=str(WORK)) as directory:
            old = Path(directory)/'control.c'
            old.write_bytes(source.read_bytes().replace(b'M.y[0,1]', b'M.y[1]'))
            with self.assertRaises(AssertionError):
                check(old)

    def test_restore_csv_header_is_optional_but_shape_is_checked(self):
        record = '40010,400.1,23.77,13.53,0,0,1,0,1,0,0\n'
        header = ('step,time,V_total,V_excess,rel_drift,max_rel_drift,picard_iter,'
                  'picard_res,mg_i,mg_resa_after,mg_resb_before\n')
        with tempfile.TemporaryDirectory(dir=str(WORK)) as directory:
            path = Path(directory)/'volume.csv'
            path.write_text(record)
            raw = rows(path)
            path.write_text(header+record)
            self.assertEqual(rows(path), raw)
            path.write_text(record.rstrip()+',unexpected\n')
            with self.assertRaises(AssertionError):
                rows(path)

    def test_native_tree_parser(self):
        self_check()


if __name__ == '__main__':
    unittest.main()
