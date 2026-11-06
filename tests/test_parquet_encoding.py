"""The parquet output, and the encoding of the output files.

Run with: python -m unittest discover -s tests -p "test_*.py"
"""
import os
import tempfile
import unittest

import pyarrow.parquet as pq

from pygridmap import gridtiler

# with non-ASCII values, empty values, and two tiles
ROWS = "x,y,pop,share,type\n0,0,1,0.5,A\n1000,0,2,1.5,Ä\n0,1000,4,,東京\n130000,0,7,3,\n"


class TestFiles(unittest.TestCase):

    def setUp(self):
        self.dir = tempfile.TemporaryDirectory()
        self.input = os.path.join(self.dir.name, "in.csv")
        with open(self.input, 'w', encoding='utf-8', newline='') as f:
            f.write(ROWS)
        self.tiles = os.path.join(self.dir.name, "tiles")
        self.output = os.path.join(self.dir.name, "out.csv")

    def tearDown(self):
        self.dir.cleanup()

    def tile(self, name="0.parquet", x="0"):
        return os.path.join(self.tiles, x, name)

    def files(self):
        return sorted(os.path.relpath(os.path.join(r, f), self.tiles).replace("\\", "/") for r, _, fs in os.walk(self.tiles) for f in fs)

    def read_bytes(self, path):
        with open(path, 'rb') as f:
            return f.read()

    def test_csv_tile_is_utf8(self):
        gridtiler.grid_tiling(self.input, self.tiles, 1000)
        text = self.read_bytes(self.tile("0.csv")).decode("utf-8")
        self.assertIn("Ä", text)
        self.assertIn("東京", text)

    def test_aggregation_is_utf8(self):
        gridtiler.grid_aggregation(self.input, 1000, self.output, 2, aggregation_fun={"type": lambda vs, _: "".join(vs)})
        self.assertIn("AÄ東京", self.read_bytes(self.output).decode("utf-8"))

    def test_transformation_is_utf8(self):
        gridtiler.grid_transformation(self.input, lambda c: None, self.output)
        text = self.read_bytes(self.output).decode("utf-8")
        self.assertIn("東京", text)

    def test_parquet_replaces_the_csv_tiles(self):
        gridtiler.grid_tiling(self.input, self.tiles, 1000, format="parquet")
        self.assertEqual(self.files(), ["0/0.parquet", "1/0.parquet", "info.json"])

    def test_parquet_values(self):
        gridtiler.grid_tiling(self.input, self.tiles, 1000, format="parquet")
        rows = pq.read_table(self.tile()).to_pylist()
        self.assertEqual(rows, [
            {"x": 0, "y": 0, "pop": 1, "share": 0.5, "type": "A"},
            {"x": 1, "y": 0, "pop": 2, "share": 1.5, "type": "Ä"},
            {"x": 0, "y": 1, "pop": 4, "share": None, "type": "東京"},
        ])

    def test_parquet_compressions(self):
        expected = {"snappy": "SNAPPY", "gzip": "GZIP", "brotli": "BROTLI", "lz4": "LZ4", "zstd": "ZSTD", None: "UNCOMPRESSED"}
        for codec, name in expected.items():
            with self.subTest(codec=codec):
                gridtiler.grid_tiling(self.input, self.tiles, 1000, format="parquet", parquet_compression=codec, clean_output_folder=True)
                self.assertEqual(pq.read_metadata(self.tile()).row_group(0).column(0).compression, name)
                self.assertEqual(pq.read_table(self.tile()).to_pylist()[1]["type"], "Ä")

    def test_parquet_unknown_compression(self):
        with self.assertRaises(Exception):
            gridtiler.grid_tiling(self.input, self.tiles, 1000, format="parquet", parquet_compression="lzo")

    def test_parquet_file_extension(self):
        gridtiler.grid_tiling(self.input, self.tiles, 1000, format="parquet", file_extension="pq")
        self.assertEqual(self.files(), ["0/0.pq", "1/0.pq", "info.json"])

    def test_csv_to_parquet_can_keep_the_csv(self):
        gridtiler.grid_tiling(self.input, self.tiles, 1000)
        gridtiler.csv_to_parquet(self.tiles, clean=False)
        self.assertEqual(self.files(), ["0/0.csv", "0/0.parquet", "1/0.csv", "1/0.parquet", "info.json"])


if __name__ == "__main__":
    unittest.main()
