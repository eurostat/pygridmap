#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
.. _gridding

Grid making operations.

**Dependencies**

*require*:      :mod:`numpy`, :mod:`pandas`, :mod:`geopandas`, :mod:`shapely`

*optional*:     :mod:`multiprocessing`

*call*:         :mod:`pygridmap.base`

**Contents**
"""

# *credits*:      `gjacopo <jacopo.grazzini@ec.europa.eu>`_ 
# *since*:        Jun 2020


#%% Settings     

from six import integer_types, string_types

import math
import numpy as np

try:
    import pandas as pd
except ImportError:
    raise IOError("!!! Error importing pandas - this package is required !!!")

try:
    import geopandas as gpd
except ImportError:
    raise IOError("!!! Error importing geopandas - this package is required !!!")
        
try:
    import multiprocessing as mp
    from queue import Empty
except: 
    pass

from pygridmap.base import GridProcessor
from pygridmap.base import NPROCESSES, NCPUS

DEFPROJ             = "EPSG:4326" # "EPSG:3035"
"""Default projection for grid making.
"""


#%% Core functions/classes

#==============================================================================
# Class GridMaker
#==============================================================================

class GridMaker(GridProcessor):
    
    MODES = ['prll', 'seq', 'qtree'] # overide GridProcessor.MODES
    SORTS = ['tile', 'tilerc', 'tilecr', 'rc', 'cr']
    
    COL_X, COL_Y = '__x__', '__y__'
    COL_INTERSECTS, COL_WITHIN = '__intersects__', '__within__'
    COL_TILE = '__tile__'
    
    #/************************************************************************/
    def __init__(self, **kwargs):
        self.__xypos = None
        self.__asc = None
        self.__processor = None
        self.mode = kwargs.pop('mode', 'prll') # self.MODES[0]
        super(GridMaker,self).__init__(**kwargs)
        self.buffer = kwargs.pop('buffer', None)
        self.sorted = kwargs.pop('sorted', False)
        self.asc = kwargs.pop('asc', True)
        self.xypos = kwargs.pop('xypos ', None)
    
    #/************************************************************************/
    @property
    def mode(self):
        return self.__mode
    @mode.setter
    def mode(self, mode):
        try:
            assert (mode in self.MODES)
        except: raise TypeError("Wrong format for processing mode")
        if mode == 'prll':
            self.__processor = self.prll_process_tile
        elif mode == 'qtree':
            self.__processor = self.qtree_process_tile
        elif mode == 'seq':
            self.__processor = self.prll_process_tile
        self.__mode = mode
       
    #/************************************************************************/
    @property
    def xypos(self):
        return self.__xypos
    @xypos.setter
    def xypos(self, xypos):        
        try:
            assert (xypos is None or xypos in self.XYPOS)
        except: raise TypeError("Wrong format for (Y,Y) coordinates location in the grid cell: '%s'" % xypos)
        if xypos is None:
            xypos = self.XYPOS[0] # arbitrary default choice: 'LLc'
        self.__xypos = xypos
       
    #/************************************************************************/
    @property
    def asc(self):
        return self.__asc
    @asc.setter
    def asc(self, asc):        
        try:
            assert (isinstance(asc, bool)                                                   \
                or (isinstance(asc, (tuple,list)) and all([isinstance(a,bool) for a in asc])))
        except: raise TypeError("Wrong format for ascending parameter: '%s'" % asc)
        if asc is True and self.sorted!='tile':
            asc = [True,True] if self.sorted in ('rc','cr') else [True,True,True] 
        self.__asc = asc
     
    #/************************************************************************/
    @classmethod
    def prll_process_tile(cls, idx, gridbbox, cellsize, tilesize, xypos,
                          mask, crs, interior, sort, trim, crop, buffer):
        iy, ix = idx[:2]
        if len(idx)>2: nytiles, nxtiles = idx[2:]
        # retrieve the indexes of the grid cells within the considered tile
        tilebbox = cls.get_tile_bbox(idx, cellsize, tilesize, gridbbox, crop)
        # build the tile
        tile = gpd.GeoDataFrame(geometry = [cls.bbox_to_polygon(*tilebbox)], 
                                crs = crs)
        if buffer is not None:
            tile = tile.buffer(0)
        # test the overall tile
        ians, wans = None, 0
        if isinstance(mask, bool):
            if mask is True: ians = 1
        elif not mask is None:
            # s = list(mask.sindex.intersection(tilebounds))
            wans = tile.within(mask).tolist()[0] and 1 # note that wans = 1 => ians = 1
            ians = 1 if wans in (True,1) else tile.intersects(mask).tolist()[0]
        if ians in (True, 1):
            rows, cols = cls.get_pos_location(cellsize, tilebbox, xypos = xypos)
            tile = gpd.GeoDataFrame({'geometry': cls.build_from_pos(cellsize, rows, cols, xypos = xypos),
                                     cls.COL_X: [x for x in cols for y in rows],
                                     cls.COL_Y: [y for x in cols for y in rows]
                                    }, 
                                    # columns = ['geometry', cls.COL_X, cls.COL_Y, cls.COL_INTERSECTS, cls.COL_WITHIN]
                                    crs = crs
                                   ).set_geometry('geometry')
        else:
            tile[cls.COL_X], tile[cls.COL_Y] = tilebbox[:2] 
        # for all tiles indifferently: set the intersects/within "indicators" and the tile index
        tile[cls.COL_INTERSECTS], tile[cls.COL_WITHIN] = ians or 0, wans or 0#np.nan
        tile[cls.COL_TILE] = ix + iy*nxtiles                                          \
            if (sort is False or sort.startswith('tile') or sort.endswith == 'rc')  \
            else ix*nytiles + iy      
        if trim is False or wans in (True,1):
            return tile
        # from now on, trim is necessarily True
        if ians in (0, None,False): # i.e., ians != 1
            return None # neither within, nor intersect: trim!
        # now build the intersects/within "indicators" for all the cells in the tile
        tile[cls.COL_INTERSECTS], tile[cls.COL_WITHIN] = False, False
        #from functools import reduce
        #tile['__intersects__'] = reduce(lambda x,y : x | y, 
        #                                [tile.intersects(mask.geometry.iloc[i] for i in range(len(mask.geometry)))])
        #ibid with '__within__': we can't use lambda functions in multiprocesses (not pickable)!
        # for within/intersects operations, see https://github.com/geopandas/geopandas/issues/317
        for i in range(len(mask.geometry)):
            tile[cls.COL_WITHIN] = tile[cls.COL_WITHIN] | tile.within(mask.geometry.iloc[i])
            tile[cls.COL_INTERSECTS] = tile[cls.COL_INTERSECTS] | tile.intersects(mask.geometry.iloc[i])
            #if tile['__intersects__'].all() and tile['__within__'].all(): break
        #tile[['__intersects__','__within__']] = 1 * tile[['__intersects__','__within__']]
        if interior is True:
            return tile[tile[cls.COL_WITHIN]]
        else: 
            return tile[tile[cls.COL_WITHIN] | tile[cls.COL_INTERSECTS]]

    #/************************************************************************/
    @classmethod
    def qtree_process_tile(cls, idx, gridbbox, cellsize, tilesize, xypos,
                           mask, crs, interior, sort, trim, crop, buffer):
        iy, ix = idx[:2]
        if len(idx)>2: nytiles, nxtiles = idx[2:]
        height, width = cellsize
        # define the bounds of the tile: retrieve the grid cells within the considered tile
        tilebox = cls.get_tile_bbox(idx, cellsize, tilesize, gridbbox, crop)
        cells = []
        def qtree_recurse(xmin, ymin, xmax, ymax):  
            cheight, cwidth = ymax - ymin, xmax - xmin
            if cheight < height and cwidth < width: return # do nothing
            #elif cheight < height:                  cheight = height
            #elif cwidth < width:                    cwidth = width
            bbox = xmin, ymin, xmin+cwidth, ymin+cheight
            # build the area tile
            tile = gpd.GeoDataFrame(geometry = [cls.bbox_to_polygon(*bbox)], 
                                    crs = crs)   
            if buffer is not None:
                tile = tile.buffer(buffer)
            # test the whole tile
            ians, wans = None, 0
            if isinstance(mask, bool):
                if mask is True:
                    ians = 1
            elif not mask is None:
                for i in range(len(mask.geometry)):
                    if wans in (False,0):
                        wans = tile.within(mask.geometry.iloc[i]).tolist()[0] # note that wans = 1 => ians = 1
                    if ians in (False,None,0):
                        ians = 1 if wans in (True,1) else tile.intersects(mask.geometry.iloc[i]).tolist()[0]
                    if ians in (True,1) and wans in (True,1):
                        break
            if wans in (True,1):
                rows, cols  = cls.get_pos_location(cellsize, bbox, xypos = xypos)
                tile = gpd.GeoDataFrame({'geometry': cls.build_from_pos(cellsize, rows, cols, xypos = xypos),
                                         cls.COL_X: [x for x in cols for y in rows],
                                         cls.COL_Y: [y for x in cols for y in rows],
                                        }, 
                                        crs = crs
                                       ).set_geometry('geometry')
            else:
                tile[cls.COL_X], tile[cls.COL_Y] = tilebox[:2]
            # for all processes tiles with not distinction
            tile[cls.COL_INTERSECTS], tile[cls.COL_WITHIN] = ians or 0, wans or 0#np.nan
            tile[cls.COL_TILE] = ix + iy*nxtiles                                          \
                if (sort is False or sort.startswith('tile') or sort.endswith == 'rc')  \
                else ix*nytiles + iy      
            if wans in (True,1)                                                         \
                    or (interior is False and ians in (True,1) and cheight <= height and cwidth <= width) \
                    or trim is False:
                cells.append(tile)
                return            
            elif ians in (False,0,None) and trim is True:
                return
            # at this stage: wans in (False,0), ians in (True,1) and trim is True
            cheight, cwidth = cheight/2, cwidth/2
            qtree_recurse(xmin,        ymin,         xmin+cwidth,   ymin+cheight)
            qtree_recurse(xmin+cwidth, ymin+cheight, xmax,          ymax)
            qtree_recurse(xmin+cwidth, ymin,         xmax,          ymin+cheight)
            qtree_recurse(xmin,        ymin+cheight, xmin+cwidth,   ymax)
        qtree_recurse(*tilebox)
        if cells == []:
            return None
        return pd.concat(cells, axis=0, ignore_index=True)    
    
    #/************************************************************************/
    def __call__(self, bbox = None, mask = None, crs = None, interior = False, 
                 trim = True, crop = False, drop = False):
        # check that data are actually parsed...
        try:
            assert not (bbox is None and mask is None)
        except: raise TypeError("No input data parsed: bbox or mask need to be informed")
        # check bounding box
        try:
            assert (bbox is None or 
                    (isinstance(bbox, (tuple,list)) and all([np.isscalar(b) for b in bbox])))
        except: raise TypeError("Wrong format for grid bounds parameter")
        if isinstance(bbox, (tuple,list)):
            try:
                assert (len(bbox)==4 and bbox[0]<bbox[2] and bbox[1]<bbox[3])
            except: raise IOError("Grid bounding box parameter not recognised")
        # check mask
        try:
            assert (mask is None or isinstance(mask, (bool, pd.DataFrame, gpd.GeoSeries, gpd.GeoDataFrame)))
        except: raise TypeError("Wrong format for mask data")	
        # check crs flag
        try:
            assert (crs is None or isinstance(crs, (string_types, integer_types)))
        except: raise TypeError("Wrong format for projection")
        # update bbox/mask/rcs accordingly
        if crs is None:
            try:
                crs = mask.crs
            except: raise IOError("No information for output grid reference system")
        elif isinstance(crs, integer_types):
            crs = str(crs) 
        if not crs.startswith("EPSG"):
            crs = "EPSG:%s" % crs 
        if isinstance(mask, pd.DataFrame): 
            try: # assuming there is a geometry
                assert ("geometry" in mask.columns)
            except: raise IOError("Geometry of mask data not recognised")
            else:
                mask = gpd.GeoDataFrame(data = mask, crs = crs)
        # update bbox
        if bbox is None:
            bbox = mask.total_bounds.tolist()
        # check interior flag
        try:
            assert isinstance(interior, bool)
        except: raise TypeError("Wrong format for interior flag")
        if interior is True:
            trim = True # we force it
        # check trimming flag
        try:
            assert isinstance(trim, bool)
        except: raise TypeError("Wrong format for trimming flag")                  
        # check crop flag
        try: # hardly used...
            assert isinstance(crop, bool)
        except: raise TypeError("Wrong format for cropping flag")   
        # check drop parameter
        __drops = [self.COL_TILE, self.COL_X, self.COL_Y, 
                   self.COL_INTERSECTS, self.COL_WITHIN]
        try:
            assert (drop is None or isinstance(drop, bool)                                  \
                or (isinstance(drop,str) and drop in __drops)                               \
                or (isinstance(drop,(tuple,list)) and all([d in __drops for d in drop])))
        except: raise TypeError("Wrong format for drop parameter")
        if drop is None:              drop = False
        elif drop is True:            drop = __drops # [self.COL_TILE, self.COL_INTERSECTS, self.COL_WITHIN] #
        elif isinstance(drop,str):    drop = [drop, ]
        # settings
        cellsize, tilesize = self.cell, self.tile
        sort = self.sorted
        xypos = self.xypos
        # update tiling processing parameters
        if np.isscalar(tilesize):
            tileshape = self.set_tile_shape(tilesize)
            tilesize = self.get_tile_size(cellsize, tileshape, bbox)
        else: # if isinstance(tilesize, (tuple,list)):
            # tilesize: unchanged
            tileshape = self.get_tile_shape(cellsize, tilesize, bbox)
        nytiles, nxtiles = tileshape
        # verify processing mode
        if self.mode == 'qtree' and not all([math.log(t,2).is_integer() for t in tilesize]):
            raise IOError('Quadtree algorithm requires tile size to be a power of 2')
        elif nytiles * nxtiles > 1 and self.mode != 'seq' and self.cores > 1:
            # warnings.warn('No parallel processing available for unique tile')
            cores, processor = 1, self.prll_process_tile
        else:
            cores, processor = self.cores, self.processor
        # start processing
        pool = mp.Pool(processes = cores)
        grid_tiles = [pool.apply_async(processor,
                                       args = ([iy, ix, nytiles, nxtiles],
                                               bbox, cellsize, tilesize, xypos,
                                               mask, crs, interior,
                                               sort, trim, crop, None))
                        for iy in range(nytiles) for ix in range(nxtiles)]
        pool.close()
        pool.join()
        grid = pd.concat([t.get() for t in grid_tiles if t is not None], 
                         axis=0, ignore_index=True)
        if self.sorted != False: 
            if sort.startswith('tile'):   bycols = [cls.COL_TILE]
            else:                         bycols = []
            if sort.endswith == 'rc':     bycols.extend([cls.COL_X, cls.COL_Y])
            elif sort.endswith == 'cr':   bycols.extend([cls.COL_Y, cls.COL_X])
            grid.sort_values(by = bycols, axis = 0, ascending=asc, inplace = True, 
                             ignore_index = True)
        if drop != False: 
            grid.drop(columns = drop, axis = 1, inplace = True, errors = 'ignore')
        if False and nytiles * nxtiles == 1: # and self.COL_TILE in grid.columns:
            grid.drop(columns = self.COL_TILE, axis = 1, inplace = True, errors = 'ignore')
        return gpd.GeoDataFrame(grid, crs = crs).set_geometry('geometry')     


#%% Main for binary usage
