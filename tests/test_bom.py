"""A BOM at the start of the input CSV file must not end up in a column name of the output.

Run with: python -m unittest discover -s tests -p "test_*.py"
"""
import csv
import glob
import os
import tempfile
import unittest

from pygridmap import gridtiler

# a CSV exported as "UTF-8 with BOM" (as Excel and many Windows tools do), with the value column first
ROWS = "number_of_addresses,x,y\r\n5,1000,2000\r\n7,1250,2000\r\n"


def header_of(path):
    with open(path, 'rb') as f:
        return f.readline().rstrip(b"\r\n")


class TestInputWithBOM(unittest.TestCase):

    def setUp(self):
        self.dir = tempfile.TemporaryDirectory()
        self.input = os.path.join(self.dir.name, "in.csv")
        with open(self.input, 'w', encoding='utf-8-sig', newline='') as f:
            f.write(ROWS)

    def tearDown(self):
        self.dir.cleanup()

    def test_tiling(self):
        out = os.path.join(self.dir.name, "tiles")
        gridtiler.grid_tiling(self.input, out, 250, crs="EPSG:3035")
        tiles = glob.glob(os.path.join(out, "**", "*.csv"), recursive=True)
        self.assertTrue(tiles)
        for tile in tiles:
            self.assertEqual(header_of(tile), b"x,y,number_of_addresses")

    def test_transformation(self):
        out = os.path.join(self.dir.name, "out.csv")
        gridtiler.grid_transformation(self.input, lambda c: None, out)
        self.assertEqual(header_of(out), b"x,y,number_of_addresses")

    def test_aggregation(self):
        out = os.path.join(self.dir.name, "agg.csv")
        gridtiler.grid_aggregation(self.input, 250, out, 2)
        with open(out, newline='', encoding='utf-8') as f:
            columns = next(csv.reader(f))
        self.assertIn("number_of_addresses", columns)
        self.assertFalse([c for c in columns if "﻿" in c or "ï»¿" in c])


class TestInputWithoutBOM(unittest.TestCase):

    def test_tiling_is_unchanged(self):
        with tempfile.TemporaryDirectory() as d:
            src = os.path.join(d, "in.csv")
            with open(src, 'w', encoding='utf-8', newline='') as f:
                f.write(ROWS)
            gridtiler.grid_tiling(src, os.path.join(d, "tiles"), 250, crs="EPSG:3035")
            tile = glob.glob(os.path.join(d, "tiles", "**", "*.csv"), recursive=True)[0]
            self.assertEqual(header_of(tile), b"x,y,number_of_addresses")


if __name__ == "__main__":
    unittest.main()
