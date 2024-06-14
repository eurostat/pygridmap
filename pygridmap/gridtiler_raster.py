#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
.. _gridtiler_raster

.. Links
.. _Eurostat: http://ec.europa.eu/eurostat/web/main
.. |Eurostat| replace:: `Eurostat <Eurostat_>`_

Tile gridded data from raster files for visualisation with GridViz javascript library.

**Dependencies**

*require*:      :mod:`os`, :mod:`csv`, :mod:`math`, :mod:`json`, :mod:`rasterio`, :mod:`rasterio.transform `, :mod:`panda`

**Contents**
"""

# *credits*:      `jgaffuri <julien.gaffuri@ec.europa.eu>`_ 
# *since*:        May 2024

#%% Settings     

import rasterio
from math import ceil,floor
import os
import csv
import json
import pandas as pd
import concurrent.futures
from datetime import datetime


def tiling_raster(rasters, output_folder, crs="", tile_size_cell=128, format="csv", parquet_compression="snappy", num_processors_to_use=1):
    """Tile gridded statistics from raster files. Note: all raster files should be based on the same gridded system: same resolution, same size, same origin point.

    Args:
        rasters (dict): A dictionnary with all data on the attributes and the raster file they are retrieved from.
        output_folder (str): The path to the output folder where to store the tiles.
        crs (str, optional): A text describing the grid CRS. Defaults to "".
        tile_size_cell (int, optional): The size of a tile, in number of cells. Defaults to 128.
        format (str, optional): The output file encodings format, either "csv" of "parquet". Defaults to "csv".
        parquet_compression (str, optional): The parquet compression. Be aware gridviz-parquet supports only snappy encodings, currently. Defaults to "snappy".

    Returns:
        _type_: _description_
    """

    #prepare variable
    resolution = None
    bounds = None
    width = None
    height = None

    #prepare and load raster file data
    for label in rasters:
        raster = rasters[label]
        #open file
        src = rasterio.open(raster["file"])

        #get base information on the rasters
        if resolution==None: resolution = src.res[0]
        if bounds==None: bounds = src.bounds
        if width==None: width = src.width
        if height==None: height = src.height

        raster["src"] = src
        raster["nodata"] = src.meta["nodata"]
        raster["data"] = src.read(raster["band"])
        if not "no_data_values" in raster: raster["no_data_values"] = []

    #bounds
    x_min, y_min, x_max, y_max = bounds.left, bounds.bottom, bounds.right, bounds.top

    #get keys
    keys = rasters.keys()

    #function to make cell template
    def build_cell(x, y):
        c = { "x":x, "y":y }
        for k in keys: c[k] = None
        return c

    #function to make a tile
    def make_tile(xyt):
        [xt, yt] = xyt
        print(datetime.now(), "tile", xt, yt)

        #prepare tile cells
        cells_index = {}

        #prepare raster data query window
        min_col = xt * tile_size_cell
        #min_row = yt * tile_size_cell
        min_row = height - (yt+1) * tile_size_cell
        #col first, then row: window(col, row, w, h)
        window = rasterio.windows.Window(min_col, min_row, tile_size_cell, tile_size_cell)
        #print(min_col, min_row)

        #handle every raster
        for key in keys:
            #print(key)

            #get raster
            raster = rasters[key]
            src = raster["src"]

            #get tile data for key
            data = src.read(raster["band"], window=window)

            #get dimensions of the data matrix
            data_height, data_width = data.shape
            dh = tile_size_cell - data_height

            #make cells
            for col in range(0, data_width):
                for row in range(0, data_height):

                    #get value
                    value = data[row, col]

                    #if no value, skip
                    if value == raster["nodata"] or value in raster["no_data_values"]: continue

                    #get cell from index. if it does not exists, create it
                    if col in cells_index: col_ = cells_index[col]
                    else: col_ = {}; cells_index[col] = col_
                    if row in col_: cell = col_[row]
                    else: cell = build_cell(col, tile_size_cell - row - 1 - dh); col_[row] = cell

                    #set cell value
                    cell[key] = value

        #get cells as a list
        cells = [cell for col in cells_index.values() for cell in col.values()]
        del cells_index

        print(datetime.now(), "tile", xt, yt, "-", len(cells), "cells")

        #if no cell within tile, skip
        if len(cells) == 0: return

        #remove column with all values null
        #check columns
        for key in keys:
            #check if cells all have key as column
            toRemove = True
            for c in cells:
                if c[key]==None: continue
                toRemove = False
                break
            #remove column
            if toRemove:
                for c in cells: del c[key]

        #sort cells
        cells = sorted(cells, key=lambda d: (d['x'], d['y']))

        #make csv header, ensuring x and y are first columns
        headers = list(cells[0].keys())
        headers.remove("x")
        headers.remove("y")
        headers.insert(0, "x")
        headers.insert(1, "y")

        #create output folder, if it does not already exists
        fo = output_folder + "/" + str(xt) + "/"
        if not os.path.exists(fo): os.makedirs(fo)

        #save as CSV file
        cfp = fo + str(yt) + ".csv"
        with open(cfp, 'w', newline='') as csv_file:
            #get writer
            writer = csv.DictWriter(csv_file, fieldnames=headers)
            #write the header
            writer.writeheader()

            #write the cell rows
            for c in cells:
                writer.writerow(c)

        if format == "csv": return

        #csv to parquet

        #load csv file            
        df = pd.read_csv(cfp)
        #save as parquet            
        df.to_parquet(fo + str(yt) + ".parquet", engine='pyarrow', compression=parquet_compression, index=False)
        #delete csv file
        os.remove(cfp)






    #tile frame caracteristics
    tile_size_geo = resolution * tile_size_cell
    tile_min_x = 0 #floor( (x_min - x_origin) / tile_size_geo )
    tile_min_y = 0 #floor( (y_min - y_origin) / tile_size_geo )
    tile_max_x = floor( (x_max - x_min) / tile_size_geo )
    tile_max_y = floor( (y_max - y_min) / tile_size_geo )


    #write info.json file
    data = {
        "dims": [],
        "crs": crs,
        "tileSizeCell": tile_size_cell,
        "originPoint": {
            "x": x_min,
            "y": y_min
        },
        "resolutionGeo": resolution,
        "tilingBounds": {
            "xMin": 0,
            "yMin": 0,
            "xMax": tile_max_x,
            "yMax": tile_max_y
        }
    }

    if not os.path.exists(output_folder): os.makedirs(output_folder)

    with open(output_folder + '/info.json', 'w') as json_file:
        json.dump(data, json_file, indent=3)



    #make list of tiles x,y
    pairs = []
    for xt in range(tile_min_x, tile_max_x+1):
        for yt in range(tile_min_y, tile_max_y+1):
            pairs.append([xt, yt])

    #make tiles, in parallel
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_processors_to_use) as executor:
        { executor.submit(make_tile, tile): tile for tile in pairs }
















def tiling_raster_generic(rasters, output_folder, resolution_out, x_min, y_min, x_max, y_max, x_origin=None, y_origin=None, crs="", tile_size_cell=128, format="csv", parquet_compression="snappy", num_processors_to_use=1):
    """Tile gridded statistics from raster files. Here, the input file do not have to be necessary based on the same grid frame. But they need to be in the same CRS.

    Args:
        rasters (dict): A dictionnary with all data on the attributes and the raster file they are retrieved from.
        output_folder (str): The path to the output folder where to store the tiles.
        resolution_out (float): The resolution of the output grid in the CRS UoM (usually meters).
        x_min (float): The extent to be tiled
        y_min (float): The extent to be tiled
        x_max (float): The extent to be tiled
        y_max (float): The extent to be tiled
        x_origin (float, optional): The origin position - if not specified, x_min is used. Defaults to None.
        y_origin (float, optional): The origin position - if not specified, y_min is used. Defaults to None.
        crs (str, optional): A text describing the grid CRS. Defaults to "".
        tile_size_cell (int, optional): The size of a tile, in number of cells. Defaults to 128.
        format (str, optional): The output file encodings format, either "csv" of "parquet". Defaults to "csv".
        parquet_compression (str, optional): The parquet compression. Be aware gridviz-parquet supports only snappy encodings, currently. Defaults to "snappy".

    Returns:
        _type_: _description_
    """

    #set origin, if not specified
    if x_origin==None: x_origin=x_min
    if y_origin==None: y_origin=y_min

    #prepare and load raster file data
    for label in rasters:
        raster = rasters[label]
        #open file
        src = rasterio.open(raster["file"])
        raster["src"] = src
        #raster["mmap"] = src.read(1, out_shape=(1, src.height, src.width), masked=True)
        raster["nodata"] = src.meta["nodata"]
        raster["data"] = src.read(raster["band"])
        if not "no_data_values" in raster: raster["no_data_values"] = []

    #tile frame caracteristics
    tile_size_geo = resolution_out * tile_size_cell
    tile_min_x = floor( (x_min - x_origin) / tile_size_geo )
    tile_min_y = floor( (y_min - y_origin) / tile_size_geo )
    tile_max_x = ceil( (x_max - x_origin) / tile_size_geo )
    tile_max_y = ceil( (y_max - y_origin) / tile_size_geo )

    r2 = resolution_out/2

    #get keys
    keys = rasters.keys()

    #function to make cell template
    def build_cell():
        c = {}
        for k in keys: c[k] = None
        return c

    #function to make a tile
    def make_tile(xyt):
        [xt, yt] = xyt
        print(datetime.now(), "tile", xt, yt)

        #prepare tile cells
        cells = []

        for xtc in range(0, tile_size_cell):
            for ytc in range(0, tile_size_cell):
                #print("cell", xtc, ytc)

                #new cell
                cell = None

                #compute geo coordinate
                xc = x_origin + xt * tile_size_geo + xtc*resolution_out
                yc = y_origin + yt * tile_size_geo + ytc*resolution_out

                #check limits
                if xc<x_min: continue
                if xc>x_max: continue
                if yc<y_min: continue
                if yc>y_max: continue

                #get values
                for key in keys:

                    #get value
                    raster = rasters[key]
                    src = raster["src"]
                    #v = raster["mmap"][src.index(xc+r2, yc+r2)]
                    row, col = src.index(xc+r2, yc+r2)
                    if col>=src.width or col<0: continue
                    if row>=src.height or row <0: continue
                    v = raster["data"][row,col]

                    if v == raster["nodata"] or v in raster["no_data_values"]: continue

                    #if not built, build cell
                    if cell == None: cell = build_cell()
                    cell[key] = v

                #no value found: skip
                if cell == None: continue

                #set cell x,y within its tile
                cell["x"] = xtc
                cell["y"] = ytc

                #store cell
                cells.append(cell)

        print(datetime.now(), "tile", xt, yt, "-", len(cells), "cells")

        #if no cell within tile, skip
        if len(cells) == 0: return

        #remove column with all values null
        #check columns
        for key in keys:
            #check if cells all have key as column
            toRemove = True
            for c in cells:
                if c[key]==None: continue
                toRemove = False
                break
            #remove column
            if toRemove:
                for c in cells: del c[key]

        #make csv header, ensuring x and y are first columns
        headers = list(cells[0].keys())
        headers.remove("x")
        headers.remove("y")
        headers.insert(0, "x")
        headers.insert(1, "y")

        #create output folder, if it does not already exists
        fo = output_folder + "/" + str(xt) + "/"
        if not os.path.exists(fo): os.makedirs(fo)

        #save as CSV file
        cfp = fo + str(yt) + ".csv"
        with open(cfp, 'w', newline='') as csv_file:
            #get writer
            writer = csv.DictWriter(csv_file, fieldnames=headers)
            #write the header
            writer.writeheader()

            #write the cell rows
            for c in cells:
                writer.writerow(c)

        if format == "csv": return xyt

        #csv to parquet

        #load csv file            
        df = pd.read_csv(cfp)
        #save as parquet            
        df.to_parquet(fo + str(yt) + ".parquet", engine='pyarrow', compression=parquet_compression, index=False)
        #delete csv file
        os.remove(cfp)

        return xyt





    #make list of tiles x,y
    pairs = []
    for xt in range(tile_min_x, tile_max_x):
        for yt in range(tile_min_y, tile_max_y):
            pairs.append([xt, yt])

    #minimum and maximum tile x,y, for info.json file
    min_tx=None
    min_ty=None
    max_tx=None
    max_ty=None

    #make tiles, in parallel
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_processors_to_use) as executor:
        tasks_to_do = {executor.submit(make_tile, tile): tile for tile in pairs}

        # merge task outputs
        for task_output in concurrent.futures.as_completed(tasks_to_do):
            xyt = task_output.result()
            if xyt!=None:
                [xt,yt]=xyt
                #store extreme positions, for info.json file
                if min_tx == None or xt<min_tx: min_tx = xt
                if min_ty == None or yt<min_ty: min_ty = yt
                if max_tx == None or xt>max_tx: max_tx = xt
                if max_ty == None or yt>max_ty: max_ty = yt


    #write info.json
    data = {
        "dims": [],
        "crs": crs,
        "tileSizeCell": tile_size_cell,
        "originPoint": {
            "x": x_origin,
            "y": y_origin
        },
        "resolutionGeo": resolution_out,
        "tilingBounds": {
            "yMin": min_ty,
            "yMax": max_ty,
            "xMax": max_tx,
            "xMin": min_tx
        }
    }

    with open(output_folder + '/info.json', 'w') as json_file:
        json.dump(data, json_file, indent=3)

