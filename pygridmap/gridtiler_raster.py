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
import numpy as np



def tiling_raster(rasters, output_folder, resolution_out, x_min, y_min, x_max, y_max, x_origin=None, y_origin=None, crs="", tile_size_cell=128, format="csv", parquet_compression="snappy", num_processors_to_use=1):
    """Tile gridded statistics from raster files.

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
        print("tile",xt,yt)

        #prepare tile cells
        cells = []


        #TODO use memory mapping ?
        #mmap = src.read(1, out_shape=(1, src.height, src.width), masked=True)
        #value = mmap[src.index(x, y)]

        #TODO: get all values of the tile in one go ?
        '''
        for key in keys:
            min_x = x_origin + xt * tile_size_geo
            max_x = x_origin + xt * tile_size_geo + (tile_size_cell-1)*resolution_out
            min_y = y_origin + yt * tile_size_geo
            max_y = y_origin + yt * tile_size_geo + (tile_size_cell-1)*resolution_out
            x_coords = np.linspace(min_x, max_x, tile_size_cell) #use range ?
            y_coords = np.linspace(min_y, max_y, tile_size_cell)
            grid_coords = [(x, y) for x in x_coords for y in y_coords]

            #convert the geographical coordinates to pixel coordinates
            src = rasters[key]["src"]
            pixel_coords = [src.index(x, y) for x, y in grid_coords]

            # Get the min and max rows and columns to determine the bounding window
            rows, cols = zip(*pixel_coords)
            min_row, max_row = min(rows), max(rows)
            min_col, max_col = min(cols), max(cols)

            # Read the window that covers all the grid points
            window = rasterio.windows.Window(min_col, min_row, max_col - min_col + 1, max_row - min_row + 1)
            data = src.read(1, window=window)
            # Extract values at the specific rows and columns
            values = [data[row - min_row, col - min_col] for row, col in pixel_coords]

            #build matrix of cell templates
            cells = []
            for xtc in range(0, tile_size_cell):
                c = []
                cells.append(c)
                for ytc in range(0, tile_size_cell):
                    c.append(build_cell())
        '''


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

