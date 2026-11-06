"""The aggregation functions, and the format of the output CSV files.

Run with: python -m unittest discover -s tests -p "test_*.py"
"""
import os
import tempfile
import unittest

import pyarrow.parquet as pq

from pygridmap import gridtiler

# a 2x2 block of cells, aggregated to a single cell with a=2
ROWS = "x,y,pop,age,type\n0,0,1,10,A\n1000,0,2,20,B\n0,1000,4,30,B\n1000,1000,9,,A\n"


class TestFunctions(unittest.TestCase):

    def test_median(self):
        self.assertEqual(gridtiler.aggregation_median(["1", "2", "4", "9"]), 3)
        self.assertEqual(gridtiler.aggregation_median(["9", "1", "4"]), 4)
        self.assertEqual(gridtiler.aggregation_median(["", "2"]), 2)

    def test_min_max(self):
        # numbers, not texts: "10" < "9"
        self.assertEqual(gridtiler.aggregation_min(["10", "9"]), 9)
        self.assertEqual(gridtiler.aggregation_max(["10", "9"]), 10)
        self.assertEqual(gridtiler.aggregation_min(["", "0", "5"]), 0)

    def test_mode(self):
        self.assertEqual(gridtiler.aggregation_mode(["A", "B", "B"]), "B")
        self.assertEqual(gridtiler.aggregation_mode(["A", "B", "B", "A"]), "A")
        self.assertEqual(gridtiler.aggregation_mode(["", "", "A"]), "A")

    def test_all_empty(self):
        for f in [gridtiler.aggregation_median, gridtiler.aggregation_min, gridtiler.aggregation_max, gridtiler.aggregation_mode]:
            self.assertIsNone(f(["", None]), f.__name__)

    def test_callable_as_aggregation_fun(self):
        # grid_aggregation calls the functions with two arguments
        for f in [gridtiler.aggregation_median, gridtiler.aggregation_min, gridtiler.aggregation_max, gridtiler.aggregation_mode]:
            self.assertEqual(f(["1", "1"], 4), f(["1", "1"]), f.__name__)


class TestFiles(unittest.TestCase):

    def setUp(self):
        self.dir = tempfile.TemporaryDirectory()
        self.input = os.path.join(self.dir.name, "in.csv")
        with open(self.input, 'w', newline='') as f:
            f.write(ROWS)
        self.output = os.path.join(self.dir.name, "out.csv")

    def tearDown(self):
        self.dir.cleanup()

    def read(self):
        with open(self.output, 'rb') as f:
            return f.read()

    def test_aggregation(self):
        gridtiler.grid_aggregation(self.input, 1000, self.output, 2, aggregation_fun={
            "pop": gridtiler.aggregation_median,
            "age": gridtiler.aggregation_max,
            "type": gridtiler.aggregation_mode,
        })
        self.assertEqual(self.read(), b"x,y,pop,age,type\r\n0,0,3,30,A\r\n")

    def test_aggregation_all_empty(self):
        with open(self.input, 'w', newline='') as f:
            f.write("x,y,a\n0,0,\n1000,0,\n")
        gridtiler.grid_aggregation(self.input, 1000, self.output, 2, aggregation_fun={"a": gridtiler.aggregation_mode})
        self.assertEqual(self.read(), b"x,y,a\r\n0,0,\r\n")

    def test_aggregation_rounding(self):
        gridtiler.grid_aggregation(self.input, 1000, self.output, 2, aggregation_rounding=0, aggregation_fun={
            "pop": gridtiler.aggregation_average_3,
            "age": gridtiler.aggregation_min,
            "type": gridtiler.aggregation_mode,
        })
        self.assertEqual(self.read(), b"x,y,pop,age,type\r\n0,0,4,10,A\r\n")

    def test_transformation_line_endings(self):
        gridtiler.grid_transformation(self.input, lambda c: None, self.output)
        self.assertEqual(self.read().splitlines(keepends=True)[0], b"x,y,pop,age,type\r\n")
        self.assertNotIn(b"\r\r", self.read())

    def test_tiling_after_aggregation(self):
        gridtiler.grid_aggregation(self.input, 1000, self.output, 2, aggregation_fun={"type": gridtiler.aggregation_mode})
        tiles = os.path.join(self.dir.name, "tiles")
        gridtiler.grid_tiling(self.output, tiles, 2000)
        with open(os.path.join(tiles, "0", "0.csv"), 'rb') as f:
            self.assertEqual(f.read(), b"x,y,pop,age,type\r\n0,0,16,60,A\r\n")

    def test_parquet_default_compression_is_readable_by_gridviz_parquet(self):
        tiles = os.path.join(self.dir.name, "tiles")
        gridtiler.grid_tiling(self.input, tiles, 1000, format="parquet")
        meta = pq.read_metadata(os.path.join(tiles, "0", "0.parquet"))
        self.assertEqual(meta.num_rows, 4)
        self.assertEqual(meta.row_group(0).column(0).compression, "SNAPPY")


if __name__ == "__main__":
    unittest.main()
