"""base, gridding and overlay, against known answers.

Run with: python -m unittest discover -s tests -p "test_*.py"
"""
import os
import unittest

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import Polygon, box

from pygridmap.base import FrameProcessor, GridProcessor
from pygridmap.gridding import GridMaker, grid_maker
from pygridmap.overlay import GridOverlay, area_interpolate

CRS = "EPSG:3035"
BBOX = [0, 0, 100000, 80000]
CELL = 10000
FULL = [box(x, y, x + CELL, y + CELL) for x in range(0, 100000, CELL) for y in range(0, 80000, CELL)]


# at the top level, as the workers of the multiprocessing pool have to import them
def frame_sum(df):
    return df.assign(c = df["a"] + df["b"])


def row_sum(row):
    return row["a"] + row["b"]


def process_id(_):
    import time
    time.sleep(0.5) # so that a worker can not process all tasks alone
    return os.getpid()


def cells_of(geometries):
    return {tuple(float(round(v, 3)) for v in g.bounds) for g in geometries}


class TestBase(unittest.TestCase):

    def setUp(self):
        self.df = pd.DataFrame({"a": range(10), "b": range(10)})
        self.a = gpd.GeoDataFrame(geometry = [box(0, 0, 2, 2)], crs = CRS)
        self.b = gpd.GeoDataFrame(geometry = [box(1, 1, 3, 3)], crs = CRS)

    def test_frame_processor(self):
        r = FrameProcessor(cores = 2)(self.df, frame_sum)
        self.assertEqual(r["c"].tolist(), (self.df["a"] + self.df["b"]).tolist())

    def test_frame_processor_on_row(self):
        r = FrameProcessor(cores = 2).on_row(self.df, row_sum)
        self.assertEqual(r.tolist(), (self.df["a"] + self.df["b"]).tolist())

    def test_frame_processor_more_cores_than_rows(self):
        r = FrameProcessor(cores = 3)(self.df.head(2), frame_sum)
        self.assertEqual(len(r), 2)

    def test_get_bbox(self):
        self.assertEqual(GridProcessor.get_bbox(self.a, self.b), [0, 0, 2, 2])
        self.assertEqual(tuple(GridProcessor.get_bbox(self.a, self.b, bbox = False)), (0, 0, 3, 3))
        self.assertEqual(len(GridProcessor.get_bbox(self.a, self.b.to_crs(4326), bbox = False)), 4)

    def test_intersection(self):
        self.assertEqual(GridProcessor.intersection([0, 0, 2, 2], [1, 1, 3, 3]), [1, 1, 2, 2])
        self.assertEqual(GridProcessor.intersection(box(0, 0, 2, 2), [1, 1, 3, 3]).area, 1)
        self.assertEqual(GridProcessor.intersection(self.a, self.b).area.sum(), 1)
        self.assertEqual(GridProcessor.intersection(self.a, self.b, rtree = True).area.sum(), 1)
        self.assertIsNone(GridProcessor.intersection(self.a, [5, 5, 6, 6], rtree = True))

    def test_unary_union(self):
        self.assertEqual(GridProcessor.unary_union(self.a).area, 4)
        self.assertEqual(GridProcessor.unary_union([0, 0, 2, 2], [1, 1, 3, 3]), [0, 0, 3, 3])
        self.assertEqual(GridProcessor.unary_union(box(0, 0, 2, 2), self.b).area, 7)
        self.assertEqual(GridProcessor.unary_union(self.a, self.b).area, 7)


class TestGridMaker(unittest.TestCase):

    def setUp(self):
        # two features: the second one is in another corner
        self.mask = gpd.GeoDataFrame({"id": [1, 2]}, crs = CRS, geometry = [
            Polygon([(0, 0), (100000, 0), (0, 80000)]),
            Polygon([(100000, 80000), (100000, 60000), (70000, 80000)]),
        ])
        self.intersects = cells_of(c for c in FULL if any(c.intersects(p) for p in self.mask.geometry))
        self.within = cells_of(c for c in FULL if any(c.within(p) for p in self.mask.geometry))

    def assert_cells(self, grid, expected):
        self.assertEqual(len(grid), len(expected), "duplicates or missing cells")
        self.assertEqual(cells_of(grid.geometry), expected)

    def test_cells_without_trimming(self):
        for mode, tile in [("seq", 1), ("prll", 4), ("prll", [3, 4]), ("qtree", [4, 4]), ("qtree", [8, 8])]:
            for cores in (1, 2):
                with self.subTest(mode = mode, tile = tile, cores = cores):
                    grid = GridMaker(mode = mode, cell = CELL, tile = tile, cores = cores)(BBOX, mask = True, crs = CRS, trim = False, crop = True)
                    self.assert_cells(grid, cells_of(FULL))
                    self.assertEqual(str(grid.crs), CRS)

    def test_mask_all_features(self):
        for mode, tile in [("seq", 1), ("prll", 4), ("prll", [2, 2]), ("qtree", [8, 8]), ("qtree", [4, 4]), ("qtree", [2, 2])]:
            for cores in (1, 2):
                with self.subTest(mode = mode, tile = tile, cores = cores):
                    make = GridMaker(mode = mode, cell = CELL, tile = tile, cores = cores)
                    self.assert_cells(make(BBOX, mask = self.mask, crs = CRS, trim = True, crop = True), self.intersects)
                    self.assert_cells(make(BBOX, mask = self.mask, crs = CRS, interior = True, crop = True), self.within)

    def test_flags_and_positions(self):
        for mode, tile in [("prll", 4), ("qtree", [8, 8])]:
            with self.subTest(mode = mode):
                grid = GridMaker(mode = mode, cell = CELL, tile = tile, cores = 1)(BBOX, mask = self.mask, crs = CRS, crop = True)
                bounds = grid.geometry.bounds
                self.assertTrue((grid[GridMaker.COL_X] == bounds.minx).all())
                self.assertTrue((grid[GridMaker.COL_Y] == bounds.miny).all())
                self.assertEqual(int(grid[GridMaker.COL_INTERSECTS].sum()), len(self.intersects))
                self.assertEqual(int(grid[GridMaker.COL_WITHIN].sum()), len(self.within))

    def test_not_square_cells(self):
        expected = cells_of(box(x, y, x + 10000, y + 20000) for x in range(0, 100000, 10000) for y in range(0, 80000, 20000))
        for mode, tile in [("prll", 2), ("qtree", [4, 8])]:
            with self.subTest(mode = mode):
                grid = GridMaker(mode = mode, cell = [20000, 10000], tile = tile, cores = 1)(BBOX, mask = True, crs = CRS, trim = False, crop = True)
                self.assert_cells(grid, expected)
                mask = gpd.GeoDataFrame(geometry = [box(0, 0, 100000, 80000)], crs = CRS)
                self.assert_cells(GridMaker(mode = mode, cell = [20000, 10000], tile = tile, cores = 1)(BBOX, mask = mask, crs = CRS, crop = True), expected)

    def test_no_mask_gives_the_tiles(self):
        grid = GridMaker(mode = "prll", cell = CELL, tile = 4, cores = 1)(BBOX, crs = CRS, trim = False, crop = True)
        self.assertEqual(len(grid), 4)
        self.assertAlmostEqual(grid.area.sum(), 100000 * 80000)

    def test_crs_formats(self):
        for crs in (3035, "3035", CRS, self.mask.crs, None):
            with self.subTest(crs = crs):
                grid = GridMaker(mode = "seq", cell = CELL, tile = 1, cores = 1)(BBOX, mask = self.mask, crs = crs)
                self.assertEqual(str(grid.crs), CRS)

    def test_sorted(self):
        make = lambda sort: GridMaker(mode = "prll", cell = CELL, tile = 4, cores = 1, sorted = sort)(BBOX, mask = True, crs = CRS, trim = False, crop = True)
        X, Y, T = GridMaker.COL_X, GridMaker.COL_Y, GridMaker.COL_TILE
        for sort, cols in [("tile", [T]), ("tilerc", [T, Y, X]), ("tilecr", [T, X, Y]), ("rc", [Y, X]), ("cr", [X, Y]), (True, [T])]:
            with self.subTest(sort = sort):
                grid = make(sort)
                self.assert_cells(grid, cells_of(FULL))
                self.assertEqual(grid[cols].values.tolist(), grid[cols].sort_values(cols).values.tolist())

    def test_xypos(self):
        for xypos, dx, dy in [("LLc", 0, 0), ("LRc", CELL, 0), ("URc", CELL, CELL), ("ULc", 0, CELL), ("CC", CELL / 2, CELL / 2), ("centre", CELL / 2, CELL / 2)]:
            with self.subTest(xypos = xypos):
                make = GridMaker(mode = "seq", cell = CELL, tile = 1, cores = 1, xypos = xypos)
                self.assertEqual(make.xypos, xypos)
                grid = make(BBOX, mask = True, crs = CRS, trim = False, crop = True)
                self.assert_cells(grid, cells_of(FULL))
                bounds = grid.geometry.bounds
                self.assertTrue(np.allclose(grid[GridMaker.COL_X], bounds.minx + dx))
                self.assertTrue(np.allclose(grid[GridMaker.COL_Y], bounds.miny + dy))

    def test_drop(self):
        make = GridMaker(mode = "seq", cell = CELL, tile = 1, cores = 1)
        self.assertEqual(list(make(BBOX, mask = True, crs = CRS, trim = False, drop = True).columns), ["geometry"])
        self.assertNotIn(GridMaker.COL_TILE, make(BBOX, mask = True, crs = CRS, trim = False, drop = GridMaker.COL_TILE).columns)

    def test_mask_outside(self):
        far = gpd.GeoDataFrame(geometry = [box(1e7, 1e7, 1e7 + 1, 1e7 + 1)], crs = CRS)
        for mode, tile in [("prll", 4), ("qtree", [4, 4])]:
            with self.subTest(mode = mode):
                grid = GridMaker(mode = mode, cell = CELL, tile = tile, cores = 1)(BBOX, mask = far, crs = CRS, trim = True)
                self.assertEqual(len(grid), 0)
                self.assertEqual(str(grid.crs), CRS)

    def test_qtree_requires_a_power_of_2(self):
        with self.assertRaises(IOError):
            GridMaker(mode = "qtree", cell = CELL, tile = [3, 3], cores = 1)(BBOX, mask = True, crs = CRS)

    def test_grid_maker(self):
        self.assert_cells(grid_maker(self.mask, CELL, cores = 1, tile = 4), self.intersects)
        self.assert_cells(grid_maker(self.mask, CELL, cores = 2, tile = 4), self.intersects)
        self.assert_cells(grid_maker(self.mask, CELL, cores = 1, tile = 4, interior = True), self.within)
        self.assert_cells(grid_maker(self.mask, CELL, bbox = BBOX, cores = 1, tile = 4), self.intersects)
        self.assertEqual(str(grid_maker(self.mask, CELL, crs = CRS, cores = 1, tile = 1).crs), CRS)

    def test_tasks_run_in_parallel(self):
        ids = GridProcessor.run_tasks(process_id, [(i,) for i in range(4)], 2)
        self.assertEqual(len(ids), 4)
        self.assertEqual(len(set(ids)), 2)
        self.assertNotIn(os.getpid(), ids)
        # no pool is created when it is not necessary
        self.assertEqual(GridProcessor.run_tasks(process_id, [(0,)], 2), [os.getpid()])
        self.assertEqual(set(GridProcessor.run_tasks(process_id, [(0,), (1,)], 1)), {os.getpid()})

    def test_tasks_errors_are_raised(self):
        for cores in (1, 2):
            with self.subTest(cores = cores):
                with self.assertRaises(TypeError):
                    GridProcessor.run_tasks(frame_sum, [(1, 2), (3, 4)], cores)


class TestOverlay(unittest.TestCase):

    def setUp(self):
        # two sources, each one split equally between two cells
        self.src = gpd.GeoDataFrame({"pop": [100.0, 40.0], "other": [1.0, 2.0]}, crs = CRS, geometry = [box(0, 0, 200, 100), box(0, 100, 200, 200)])
        self.tgt = gpd.GeoDataFrame({"name": list("abcd")}, crs = CRS, geometry = [box(x, y, x + 100, y + 100) for x in (0, 100) for y in (0, 100)])
        self.expected = {"a": 50.0, "b": 20.0, "c": 50.0, "d": 20.0}

    def values(self, result, col = "pop"):
        self.assertEqual(result[GridOverlay.COL_GRIDX].is_unique if GridOverlay.COL_GRIDX in result else True, True)
        return dict(zip(result["name"], result[col]))

    def test_area_interpolate(self):
        for kw in [dict(tile = 1, cores = 1), dict(tile = 1, cores = 2), dict(tile = 4, cores = 1), dict(tile = 4, cores = 2),
                   dict(tile = 9, cores = 2), dict(tile = 4, cores = 2, memory_split = True), dict(tile = 1, cores = 2, memory_split = True), dict(cores = 2)]:
            with self.subTest(**kw):
                result = area_interpolate(self.src, self.tgt, ["pop"], **kw)
                self.assertEqual(len(result), 4)
                self.assertEqual(self.values(result), self.expected)

    def test_area_interpolate_column_as_text(self):
        self.assertEqual(self.values(area_interpolate(self.src, self.tgt, "pop", tile = 1, cores = 1)), self.expected)

    def test_input_data_unchanged(self):
        columns, target_columns = list(self.src.columns), list(self.tgt.columns)
        area_interpolate(self.src, self.tgt, ["pop"], tile = 4, cores = 1)
        self.assertEqual(list(self.src.columns), columns)
        self.assertEqual(list(self.tgt.columns), target_columns)
        self.assertEqual(self.src["pop"].tolist(), [100.0, 40.0])

    def test_no_cell_counted_twice(self):
        # regular grid, the cells of which are on the borders of the tiles
        region = gpd.GeoDataFrame({"pop": [1000.0, 600.0]}, crs = CRS, geometry = [box(0, 0, 50000, 80000), box(50000, 0, 100000, 80000)])
        grid = GridMaker(mode = "prll", cell = CELL, tile = 4, cores = 1)(BBOX, mask = True, crs = CRS, trim = False, crop = True, drop = True)
        for tile, cores in [(1, 1), (1, 2), (4, 1), (4, 2), (9, 2), (16, 2), ([30000, 30000], 2)]:
            with self.subTest(tile = tile, cores = cores):
                result = area_interpolate(region, grid, ["pop"], tile = tile, cores = cores)
                self.assertEqual(len(result), 80)
                self.assertAlmostEqual(result["pop"].sum(), 1600)
                self.assertTrue(result[GridOverlay.COL_GRIDX].is_unique)

    def test_grid_with_columns_and_filtered_grid(self):
        region = gpd.GeoDataFrame({"pop": [1600.0]}, crs = CRS, geometry = [box(0, 0, 100000, 80000)])
        grid = GridMaker(mode = "prll", cell = CELL, tile = 4, cores = 1)(BBOX, mask = True, crs = CRS, trim = False, crop = True)
        kept = grid[grid[GridMaker.COL_X] >= 50000]
        self.assertNotEqual(kept.index.tolist(), list(range(len(kept))))
        for tile in (1, 4):
            with self.subTest(tile = tile):
                result = area_interpolate(region, kept, ["pop"], tile = tile, cores = 1)
                self.assertEqual(len(result), len(kept))
                self.assertEqual(sorted(result[GridOverlay.COL_GRIDX]), sorted(kept.index))
                self.assertEqual(cells_of(result.geometry), cells_of(kept.geometry))
                self.assertAlmostEqual(result["pop"].sum(), 1600 * len(kept) / 80)

    def test_tile_column(self):
        region = gpd.GeoDataFrame({"pop": [1600.0]}, crs = CRS, geometry = [box(0, 0, 100000, 80000)])
        grid = GridMaker(mode = "prll", cell = CELL, tile = 4, cores = 1)(BBOX, mask = True, crs = CRS, trim = False, crop = True)
        for how in ("column", True):
            with self.subTest(how = how):
                proc = GridOverlay(mode = "prll", how = "intersection", cores = 1, tile = GridOverlay.COL_TILE)
                result = proc(region, grid, rule = "sum", columns = ["pop"])
                self.assertEqual(len(result), 80)
                self.assertAlmostEqual(result["pop"].sum(), 1600)

    def overlay(self, rule = None, how = "intersection", **kw):
        return GridOverlay(mode = "prll", how = how, cores = 1, tile = 1)(self.src, self.tgt, rule = rule, columns = ["pop"] if rule else None, **kw)

    def test_overlay_pieces(self):
        result = self.overlay()
        self.assertEqual(len(result), 4)
        self.assertAlmostEqual(result.area.sum(), 40000)
        self.assertTrue(np.allclose(result[GridOverlay.COL_PCT], 0.5))
        self.assertAlmostEqual(GridOverlay(how = "union", cores = 1, tile = 1)(self.src, self.tgt).area.sum(), 40000)
        self.assertEqual(len(GridOverlay(mode = "seq", cores = 1, tile = 1)(self.src, self.tgt)), 4)

    def test_rules(self):
        self.assertEqual(self.values(self.overlay("sum")), self.expected)
        straddling = gpd.GeoDataFrame({"name": ["x"]}, crs = CRS, geometry = [box(0, 50, 100, 200)])
        def apply(rule):
            return GridOverlay(cores = 1, tile = 1)(self.src, straddling, rule = rule, columns = ["pop"])["pop"].iloc[0]
        self.assertAlmostEqual(apply("sum"), 100 * 0.25 + 40 * 0.5)
        self.assertEqual(apply("max"), 40) # the source overlapped the most, in share of its area
        self.assertEqual(apply("min"), 100)
        self.assertEqual(sorted(apply("list")), [40, 100])

    def test_cover_and_drop(self):
        result = self.overlay(cover = True)
        self.assertIn(GridOverlay.COL_POLIDX, result)
        self.assertNotIn(GridOverlay.COL_POLIDX, self.src)
        result = self.overlay("sum", drop = True)
        self.assertEqual(self.values(result), self.expected)
        self.assertNotIn("other", result)
        self.assertNotIn(GridOverlay.COL_TILE, result)
        result = self.overlay(drop = True)
        self.assertIn(GridOverlay.COL_AREA, result)
        self.assertIn(GridOverlay.COL_PCT, result)
        self.assertNotIn("other", result)

    def test_geoseries(self):
        result = GridOverlay(cores = 1, tile = 1)(self.src.geometry, self.tgt.geometry)
        self.assertEqual(len(result), 4)
        self.assertAlmostEqual(result.area.sum(), 40000)

    def test_nothing_overlaps(self):
        far = gpd.GeoDataFrame({"pop": [1.0]}, crs = CRS, geometry = [box(1e6, 1e6, 1e6 + 1, 1e6 + 1)])
        result = area_interpolate(far, self.tgt, ["pop"], tile = 4, cores = 1)
        self.assertEqual(len(result), 0)
        self.assertEqual(str(result.crs), CRS)

    def test_errors_are_not_hidden(self):
        with self.assertRaises(IOError):
            area_interpolate(self.src, self.tgt, ["unknown"], tile = 1, cores = 1)


if __name__ == "__main__":
    unittest.main()
