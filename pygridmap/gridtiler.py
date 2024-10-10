#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
.. _gridtiler

.. Links
.. _Eurostat: http://ec.europa.eu/eurostat/web/main
.. |Eurostat| replace:: `Eurostat <Eurostat_>`_

Tile gridded data for visualisation with GridViz javascript library.

**Dependencies**

*require*:      :mod:`os`, :mod:`csv`, :mod:`math`, :mod:`json`, :mod:`shutil`, :mod:`panda`

**Contents**
"""

# *credits*:      `jgaffuri <julien.gaffuri@ec.europa.eu>`_ 
# *since*:        May 2024

#%% Settings     

import os
import csv
from math import floor
import json
import shutil
import pandas as pd
import numbers

def grid_tiling(
    input_file,
    output_folder,
    resolution,
    tile_size_cell = 128,
    x_origin = 0.0,
    y_origin = 0.0,
    crs = "",
    clean_output_folder = False,
    input_file_delimiter = ",",
    output_file_delimiter = ",",
    format = "csv",
    parquet_compression = "snappy",
    file_extension = None
):
    """Tile an input CSV file.

    Args:
        input_file (str): The path to the input grid as CSV file. One row per grid cell. This CSV file should include a x and y columns containing the coordinates of the lower left corner of the cell. If this is not the case, use grid_transform function for example.
        output_folder (str): The path to the output folder where to store the tiles.
        resolution (float): The resolution of the input grid in the CRS UoM (usually meters).
        tile_size_cell (int, optional): The size of a tile, in number of cells. Defaults to 128.
        x_origin (float, optional): The x coordinate of the tiling scheme origin point. Defaults to 0.
        y_origin (float, optional): The y coordinate of the tiling scheme origin point. Defaults to 0.
        crs (str, optional): A text describing the grid CRS. Defaults to "".
        clean_output_folder (bool, optional): Set to True to delete the content of the output folder in the beginning: Be careful the output folder is well specified ! Defaults to False.
        input_file_delimiter (str, optional): The CSV delimiter of the input file. Defaults to ",".
        output_file_delimiter (str, optional): The CSV delimiter of the output file. Defaults to ",".
        format (str, optional): The output file encodings format, either "csv" of "parquet". Defaults to "csv".
        parquet_compression (str, optional): The parquet compression. Be aware gridviz-parquet supports only snappy encodings, currently. Defaults to "snappy".
        file_extension (str, optional): The file extension. Defaults to same as format.
    """

    #set file extension if not specified
    if file_extension == None: file_extension = format

    #compute tile size, in geo unit
    tile_size_m = resolution * tile_size_cell

    #minimum and maximum tile x,y, for info.json file
    minTX=None
    maxTX=None
    minTY=None
    maxTY=None

    #clean output folder and create it
    if clean_output_folder and os.path.exists(output_folder):
        shutil.rmtree(output_folder)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    #open file with data to tile
    with open(input_file, 'r') as csvfile:
        csvreader = csv.DictReader(csvfile, delimiter=input_file_delimiter)
        csv_header = None

        #iterate through cells from the input CSV file
        for c in csvreader:

            c["x"] = float(c["x"])
            c["y"] = float(c["y"])

            #get cell tile x,y
            xt = int(floor((c["x"] - x_origin) / tile_size_m))
            yt = int(floor((c["y"] - y_origin) / tile_size_m))

            #store extreme positions, for info.json file
            if minTX == None or xt<minTX: minTX = xt
            if maxTX == None or xt>maxTX: maxTX = xt
            if minTY == None or yt<minTY: minTY = yt
            if maxTY == None or yt>maxTY: maxTY = yt

            #compute cell position within its tile
            c["x"] = int(floor((c["x"] - x_origin) / resolution - xt * tile_size_cell))
            c["y"] = int(floor((c["y"]- y_origin) / resolution - yt * tile_size_cell))

            #check x,y values. Should be within [0,tile_size_cell-1]
            if (c["x"] < 0): print("Too low value: " + c.x + " <0")
            if (c["y"] < 0): print("Too low value: " + c.y + " <0")
            if (c["x"] > tile_size_cell - 1): print("Too high value: " + c.x + " >" + (tile_size_cell - 1))
            if (c["y"] > tile_size_cell - 1): print("Too high value: " + c.y + " >" + (tile_size_cell - 1))

            #round floats to ints, when possible, to avoid useless ".0"
            round_floats_to_ints(c)

            #create folder
            t_folder = output_folder + "/" + str(xt) + "/"
            if not os.path.exists(t_folder):
                os.makedirs(t_folder)

            #open tiled CSV file of create it if it does not exists
            file_path = t_folder + str(yt) + "." + (file_extension if format=="csv" else "csv")
            
            with open(file_path, 'a', newline='') as csvfile:

                #get header
                if csv_header==None:
                    csv_header = get_csv_header(c)

                #get writer
                writer = csv.DictWriter(csvfile, fieldnames=csv_header, delimiter=output_file_delimiter)

                #write header
                if not file_exists: writer.writeheader()

                #write cell data
                writer.writerow(c)


    #write info.json
    data = {
        "dims": [],
        "crs": crs,
        "tileSizeCell": tile_size_cell,
        "originPoint": {
            "x": x_origin,
            "y": y_origin
        },
        "resolutionGeo": resolution,
        "tilingBounds": {
            "yMin": minTY,
            "yMax": maxTY,
            "xMax": maxTX,
            "xMin": minTX
        }
    }

    with open(output_folder + '/info.json', 'w') as json_file:
        json.dump(data, json_file, indent=3)

    if format == "csv": return

    #parquet format
    csv_to_parquet(output_folder, clean=True, compression=parquet_compression, file_extension=file_extension)



def csv_to_parquet(folder_path, clean=False, compression='snappy', file_extension='parquet'):
    """Convert CSV tiled data into parquet format

    Args:
        folder_path (str): The data folder
        clean (bool, optional): Set to true to delete the initial CSV files in the end of the process. Otherwise, they will be kept. Defaults to False.
        compression (str, optional): The parquet compression. Be aware gridviz-parquet supports only snappy encodings, currently. Defaults to "snappy".
        file_extension (str, optional): The file extension. Defaults to "parquet".
    """    
    for root, _, files in os.walk(folder_path):
        for file in files:
            if not file.endswith('.csv'): continue

            csv_file_path = os.path.join(root, file)
            parquet_file_path = os.path.splitext(csv_file_path)[0] + '.' + file_extension

            #load csv file            
            df = pd.read_csv(csv_file_path)
            #save as parquet            
            df.to_parquet(parquet_file_path, engine='pyarrow', compression=compression, index=False)
            #delete csv file
            if clean: os.remove(csv_file_path)





def grid_transformation(input_file, function, output_file, input_file_delimiter = ",", output_file_delimiter = ","):
    """transform CSV data.

    Args:
        input_file (str): The path to the input grid as CSV file. One row per grid cell.
        function (dict): The function to apply to a cell row element, as a dictionnary. Properties can be added and removed within this function.
        output_file (str): The path to the output grid as CSV file. One row per grid cell.
        input_file_delimiter (str, optional): The CSV delimiter of the input file. Defaults to ",".
        output_file_delimiter (str, optional): The CSV delimiter of the output file. Defaults to ",".
    """    

    #open file to read
    with open(input_file, 'r') as infile:
        csvreader = csv.DictReader(infile, delimiter=input_file_delimiter)

        #check output file exists
        file_exists = os.path.exists(output_file)
        with open(output_file, 'w') as outfile:
            writer = None

            #iterate through cells from the input CSV file
            for c in csvreader:

                #apply function
                out = function(c)

                #filter out
                if out == False: continue

                #to save space, by removing the ".0"
                round_floats_to_ints(c)

                #create writer, if necessary, write file header
                if writer==None:
                    csv_header = get_csv_header(c)
                    writer = csv.DictWriter(outfile, fieldnames=csv_header, delimiter=output_file_delimiter)
                    if not file_exists: writer.writeheader()

                #write cell data
                writer.writerow(c)










def grid_aggregation(
    input_file,
    resolution,
    output_file,
    a,
    aggregation_fun = {},
    aggregation_rounding = 6,
    input_file_delimiter = ",",
    output_file_delimiter = ","
):
    """Aggregate grid data to a lower resolution. NB: in this version, the aggregated values are the sum of the value - other functions to be available soon.

    Args:
        input_file (str): The path to the input grid as CSV file. One row per grid cell. This CSV file should include a x and y columns containing the coordinates of the lower left corner of the cell. If this is not the case, use grid_transform function for example.
        resolution (float): The resolution of the input grid in the CRS UoM (usually meters).
        output_file (str): The path to the output grid as CSV file.
        a (int): Aggregation factor. The resolution of the output data is a multiple of the input resolution: a x resolution
        aggregation_fun (dict): Dictionnary of aggregation functions. The dictionnary keys are the cell attributes. The dictionnary values are the corresponding aggregation function to use for the column. If not specified, sum is used as aggregation function.
        aggregation_rounding (int, optional): The number of decimal places to keep for the aggregated figures. Defaults to 6.
        input_file_delimiter (str, optional): The CSV delimiter of the input file. Defaults to ",".
        output_file_delimiter (str, optional): The CSV delimiter of the output file. Defaults to ",".
    """

    # the aggregated cells, indexed by xa and then ya
    aggregation_index = {}
    target_resolution = a*resolution
    keys = None

    print("aggregation indexing...")
    with open(input_file, 'r') as infile:
        csvreader = csv.DictReader(infile, delimiter=input_file_delimiter)

        #iterate through cells from the input CSV file
        for c in csvreader:
            #get aggregated cell x,y
            xa = target_resolution * floor(float(c["x"]) / target_resolution)
            ya = target_resolution * floor(float(c["y"]) / target_resolution)

            #release memory
            del c["x"]; del c["y"]

            #store keys
            if keys==None: keys = list(c.keys())

            #add cell to its aggregation level
            try: cA_ = aggregation_index[str(xa)]
            except:
                cA_ = {}
                aggregation_index[str(xa)] = cA_
            try: cA = cA_[str(ya)]
            except:
                cA = []
                cA_[str(ya)] = cA
            cA.append(c)

    print("aggregation computation...")

    #prepare function to round aggregated figures
    tolerance = pow(10, aggregation_rounding)
    round_to_tolerance = lambda number : round(number * tolerance) / tolerance

    writer = None
    with open(output_file, 'w') as outfile:

        #aggregate cell values
        for xa, d in aggregation_index.items():
            for ya, cells in d.items():
                #print(xa,ya,len(cells))

                #make aggregated cell
                cA = { "x": xa, "y": ya }

                #compute aggregates values
                for k in keys:
                    #get list of values to aggregate
                    values = []
                    for c in cells: values.append(c[k])
                    #get aggregation function from dictionnary, or use sum by default
                    af = aggregation_fun[k] if k in aggregation_fun else aggregation_sum
                    #compute and set aggregated value
                    cA[k] = af(values, a*a)
                    if (aggregation_rounding != None and isinstance(cA[k], numbers.Number)): cA[k] = round_to_tolerance(cA[k])

                #if not, create writer and write header
                if writer == None:
                    writer = csv.DictWriter(outfile, fieldnames=get_csv_header(cA), delimiter=output_file_delimiter)
                    writer.writeheader()

                #round floats
                round_floats_to_ints(cA)

                #write aggregated cell data in output file
                writer.writerow(cA)

                #release memory immediatelly
                cells.clear()
                #del d[ya]




#some pre-defined aggregation functions

def aggregation_sum(values, _=0):
    """sum

    Args:
        values (float): The values to aggregate
        _ (int): unused

    Returns:
        float: The aggregated value
    """    
    sum = 0
    for value in values:
        if value == "" or value == None: continue
        sum += float(value)
    return sum


def aggregation_average(values, _):
    """average: the sum of values, divided by the number of values

    Args:
        values (float): The values to aggregate
        _ (int): unused

    Returns:
        float: The aggregated value
    """
    nb = len(values)
    if nb==0: return
    return aggregation_sum(values) / nb


def aggregation_average_2(values, nb):
    """average: the sum of values, divided by the total number of aggregated cells

    Args:
        values (float): The values to aggregate
        nb (int): The number of cells to aggregate (a x a)

    Returns:
        float: The aggregated value
    """
    sum = aggregation_sum(values)
    return sum/nb


def aggregation_single_value(values, _):
    """single value: return one of the cell values, the first one.

    Args:
        values (str): The values to aggregate
        _ (int): unused

    Returns:
        float: The aggregated value
    """
    return values[0]








def get_csv_header(cell):
    """Make csv header from a cell, starting with "x" and "y".

    Args:
        cell (dict): A grid cell, as a dictionnary

    Returns:
        list(str): The CSV header.
    """

    #get keys
    keys = list(cell.keys())

    #ensure x and y columns are the first ones
    if "x" in keys and "y" in keys:
        keys.remove("x")
        keys.remove("y")
        keys.insert(0, "x")
        keys.insert(1, "y")

    return keys



def round_floats_to_ints(cell):
    """ Convert rounded floats into ints so that we do not have useless ".0".
    Args:
        cell (dict): A grid cell, as a dictionnary
    """
    for key, value in cell.items():
        try:
            f = float(value)
            if f.is_integer(): cell[key] = int(f)
        except Error: pass


