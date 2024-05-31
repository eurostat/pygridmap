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


#TODO geotiff/image input

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
    parquet_compression = "snappy"
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
    """


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
            file_path = t_folder + str(yt) + ".csv"
            file_exists = os.path.exists(file_path)
            with open(file_path, 'a', newline='') as csvfile:

                #get header
                if csv_header==None:
                    csv_header = get_csv_header(c)

                #get writer
                writer = csv.DictWriter(csvfile, fieldnames=csv_header, delimiter=output_file_delimiter)

                #write header if the file was just created
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
    csv_to_parquet(output_folder, clean=True, compression=parquet_compression)



def csv_to_parquet(folder_path, clean=False, compression='snappy'):
    """Convert CSV tiled data into parquet format

    Args:
        folder_path (str): The data folder
        clean (bool, optional): Set to true to delete the initial CSV files in the end of the process. Otherwise, they will be kept. Defaults to False.
        compression (str, optional): The parquet compression. Be aware gridviz-parquet supports only snappy encodings, currently. Defaults to "snappy".
    """    
    for root, _, files in os.walk(folder_path):
        for file in files:
            if not file.endswith('.csv'): continue

            csv_file_path = os.path.join(root, file)
            parquet_file_path = os.path.splitext(csv_file_path)[0] + '.parquet'

            #load csv file            
            df = pd.read_csv(csv_file_path)
            #save as parquet            
            df.to_parquet(parquet_file_path, engine='pyarrow', compression=compression, index=False)
            #delete csv file
            if clean: os.remove(csv_file_path)





def grid_transformation(input_file, function, output_file, input_file_delimiter = ",", output_file_delimiter = ","):
    """transform CSV

    Args:
        input_file (_type_): _description_
        function (_type_): _description_
        output_file (_type_): _description_
        input_file_delimiter (str, optional): _description_. Defaults to ",".
        output_file_delimiter (str, optional): _description_. Defaults to ",".
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
                function(c)
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
    aggregation_rounding = 6,
    input_file_delimiter = ",",
    output_file_delimiter = ","
):

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

    #aggregation function
    #TODO handle other cases: average, mode, etc
    def aggregation_fun(values):
        sum = 0
        for value in values:
            if value == "" or value == None: continue
            sum += float(value)
        return sum

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
                    #compute and set aggregated value
                    cA[k] = aggregation_fun(values)
                    if (aggregation_rounding != None): cA[k] = round_to_tolerance(cA[k])

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










#make csv header from a cell, starting with "x" and "y"
def get_csv_header(cell):
    keys = list(cell.keys())
    keys.remove("x")
    keys.remove("y")
    keys.insert(0, "x")
    keys.insert(1, "y")
    return keys

#convert rounded floats into ints so that we do not have useless ".0"
def round_floats_to_ints(cell):
    for key, value in cell.items():
        try:
            f = float(value)
            if f.is_integer(): cell[key] = int(f)
        except ValueError: pass


