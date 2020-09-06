#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
.. _base

.. Links
.. _Eurostat: http://ec.europa.eu/eurostat/web/main
.. |Eurostat| replace:: `Eurostat <Eurostat_>`_

Base module enabling basic data processing and grid geoprocessing.

**Dependencies**

*require*:      :mod:`functools`, :mod:`numpy`, :mod:`pandas`, :mod:`geopandas`, :mod:`shapely`
*optional*:     :mod:`dask`, :mod:`multiprocessing`

**Contents**
"""

# *credits*:      `gjacopo <jacopo.grazzini@ec.europa.eu>`_ 
# *since*:        Jun 2020

#%% Settings     

import  functools

import json
import numpy as np

try:
    import pandas as pd
except ImportError:
    raise IOError("!!! Error importing pandas - this package is required !!!")

try:
    import dask.dataframe as dd
except:
    pass # print("!!! Error importing dask !!!")

try:
    import geopandas as gpd
except ImportError:
    raise IOError("!!! Error importing geopandas - this package is required !!!")
        
try:
    import shapely
except ImportError:
    raise IOError("!!! Error importing shapely - this package is required  !!!")
else:
    from shapely import geometry, ops, wkt

try:
    import multiprocessing as mp
    from queue import Empty
except: 
    NPROCESSES = NCPUS = 1
else:
    NPROCESSES = NCPUS = mp.cpu_count()

	
#%% Core functions/classes

#==============================================================================
# Class FrameProcessor
#==============================================================================

class FrameProcessor():

    #/************************************************************************/
    def __init__(self, cores=None):
        self.cores = cores or mp.cpu_count()
        
    #/************************************************************************/
    def __call__(self, df, func):
        pool = mp.Pool(processes=self.cores)
        df = pd.concat(
            pool.map(func, 
                     np.array_split(df, self.cores)
                    )
            )
        pool.close()
        pool.join()
        return df
    
    #/************************************************************************/
    @classmethod
    def _on_row(cls, func, df):
        #rows_iter = (row for _, row in df.iterrows())
        return df.apply(func, axis=1)

    #/************************************************************************/
    def on_row(self, df, func):
        return self.__call__(df, functools.partial(self._on_row, func))


#==============================================================================
# Class GridProcessor
#==============================================================================

class GridProcessor():
        
    MODES = ['prll', 'seq']
    SORTS = ['rc', 'cr']
    XYPOS = ['LLc', 'LRc',' URc', 'ULc', 'CC', 'centre']
    
    #/************************************************************************/
    def __init__(self, **kwargs):
        # some dumb init
        self.__mode, self.__cores = None, 1
        self.__cell, self.__tile = None, None
 	self.__xypos = None
       # self.__sorted, self.__asc = None, None
        # self.mode = kwargs.pop('mode', self.MODES[0])
        self.cores = kwargs.get('cores') or NPROCESSES
        self.cell = kwargs.get('cell')
        self.tile = kwargs.get('tile')
        self.sorted = kwargs.pop('sorted', False)
        self.asc = kwargs.pop('asc', True)
       
    #/************************************************************************/
    @property
    def cores(self):
        return self.__cores
    @cores.setter
    def cores(self, cores):
        try:
            assert isinstance(cores, int) 
        except: raise TypeError("Wrong format for cores number")
        self.__cores = cores
            
    #/************************************************************************/
    @property
    def cell(self):
        return self.__cell
    @cell.setter
    def cell(self, cell):
        try:
            assert (cell is None or isinstance(cell, int)                                  \
                or (isinstance(cell, (tuple,list)) and all([isinstance(c,int) for c in cell])))
        except: raise TypeError("Wrong format for cell size parameter")
        self.__cell = [cell, cell] if np.isscalar(cell) else cell
       
    #/************************************************************************/
    @property
    def xypos(self):
        return self.__xypos
    @xypos.setter
    def xypos(self, xypos):        
        try:
            assert (xypos is None or xypos in self.XYPOS)
        except: raise TypeError("Wrong format for (Y,Y) coordinates location in the grid cell")
        if xypos is None:
            xypos = 'LLc'
        self.__xypos = xypos
            
    #/************************************************************************/
    @property
    def tile(self):
        return self.__tile
    @tile.setter
    def tile(self, tile):
        try:
            assert (tile is None or isinstance(tile, bool) or np.isscalar(tile)         \
                or (isinstance(tile, (tuple,list)) and all([np.isscalar(t) for t in tile])))
        except: raise TypeError("Wrong format for tiling parameter")
        if isinstance(tile, bool): 
            tile = self.cores if tile is True else 1
        self.__tile = tile
    
    #/************************************************************************/
    @property
    def sorted(self):
        return self.__sorted    
    @sorted.setter
    def sorted(self, sort):        
        try:
            assert (isinstance(sort, bool) or sort in self.SORTS)
        except: raise TypeError("Wrong format for sorting parameter")
        if sort is True:              
            sort = self.SORTS[0]
        self.__sorted = sort
    
    #/************************************************************************/
    @staticmethod
    def set_ref_proj(grid, polygons):
        if grid.crs == polygons.crs:
            return polygons
        else:
            return polygons.to_crs(grid.crs)

    #/************************************************************************/
    @staticmethod
    def get_bbox(grid, polygons, gridbbox=True):
        if gridbbox is True:
            return grid.total_bounds.tolist()
        elif grid.crs == polygons.crs:
            return polygons.unary_union.union(grid.unary_union).envelope.bounds 
        else:
            return polygons.to_crs(grid.crs).unary_union.union(grid.unary_union).envelope.bounds
        
    #/************************************************************************/
    @staticmethod
    def get_grid_shape(cellsize, bbox, buffer=None):
        # Return the 'shape' of the grid covering the bounding box, i.e. the (y,x) number of unit cells 
        # in the grid
        height, width = cellsize
        buffy, buffx = buffer if buffer is not None else [0,0]
        xmin, ymin, xmax, ymax = bbox
        return [int(np.ceil((ymax-ymin+2*buffy)/height)), 
                int(np.ceil((xmax-xmin+2*buffx)/width))]  
    
    #/************************************************************************/
    @staticmethod
    def set_tile_shape(ntiles, gridshape=None):
        n = math.sqrt(ntiles)
        if gridshape is not None:
            nrows, ncols = gridshape
            ratio = max(nrows,ncols) / min(nrows, ncols)
        else: # dumb allocations
            nrows = ncols = ratio = 1
        nf, nc = int(np.floor(n/ratio)), int(np.ceil(n*ratio))
        if nf*nc>=ntiles: 
            return [nf, nc] if nrows >= ncols else [nc, nf]
        else: 
            return [nf+1, nc] if nrows >= ncols else [nc, nf+1]

    #/************************************************************************/
    @staticmethod
    def get_tile_shape(cellsize, tilesize, bbox, buffer=None):
        # Return the (x,y) shape/dimension (in # of tiles) of the tiling (made of unit cells of given size) 
        # covering the desired bounding box
        height, width = cellsize
        buffy, buffx = buffer if buffer is not None else [0,0]
        nygrid, nxgrid = tilesize
        # nrows, ncols = GridMaker.get_grid_shape(cellsize, bbox)
        xmin, ymin, xmax, ymax = bbox
        nrows, ncols = int(np.ceil((ymax-ymin+2*buffy)/height)), int(np.ceil((xmax-xmin+2*buffx)/width))  # grid shape
        return [int(np.ceil(nrows / nygrid)), 
                int(np.ceil(ncols / nxgrid))] # = #{tiles covering the grid}

    #/************************************************************************/
    @staticmethod
    def get_tile_size(cellsize, tileshape, bbox, buffer=None):
        # Return the (x,y) size (in # of unit cells of given size) of the tiling (with given shape) 
        # covering the desired bounding box
        height, width = cellsize
        buffy, buffx = buffer if buffer is not None else [0,0]
        nytile, nxtile = tileshape # tileshape = #tiles
        xmin, ymin, xmax, ymax = bbox
        nrows, ncols = int(np.ceil((ymax-ymin+2*buffy)/height)), int(np.ceil((xmax-xmin+2*buffx)/width)) # grid shape
        return [int(np.ceil(nrows / nytile)), 
                int(np.ceil(ncols / nxtile))] # = #{cells in a tile}
    
    #/************************************************************************/
    @staticmethod
    def bbox_to_polygon(west, south, east, north, density=False, buffer=False):
        # Return a polygon geometry given by the bounding box coordinates
        if density is True:     density = 10
        elif density is False:  density = 1
        elif density == 0:      density = 1
        poly = []
        poly.extend([(w, south) 
                     for w in [west + i * np.floor((east - west) / density) for i in range(density)]])
        poly.extend([(east, s)  
                     for s in [south + i * np.floor((north - south) / density) for i in range(density)]])
        poly.extend([(e, north) 
                     for e in [east - i * np.floor((east - west) / density) for i in range(density)]])
        poly.extend([(west, n)  
                     for n in [north - i * np.floor((north - south) / density) for i in range(density)]])
        poly = geometry.Polygon(poly)
        return poly if buffer is False else poly.buffer(0)
 
    #/************************************************************************/
    @classmethod
    def bbox_to_geoframe(cls, west, south, east, north, 
                         density=False, buffer=False, crs=None):
        bbox = gpd.GeoDataFrame(index=[0], 
                                geometry=[cls.bbox_to_polygon(west, south, east, north, 
                                                              density = density, buffer = buffer)
                                         ]
                               )
        if crs is not None:
            bbox.crs = crs
        return bbox
            
    #/************************************************************************/
    @staticmethod
    def get_tile_bbox(idx, cellsize, tilesize, bbox, crop, buffer=None):
        # Return the bounding box of a tile for a given index
        height, width = cellsize
        buffy, buffx = buffer if buffer is not None else [0,0]
        iy, ix = idx[:2]
        nygrid, nxgrid = tilesize
        xmin, ymin, xmax, ymax = bbox
        bxmin, bymin = xmin + ix * nxgrid * width, ymin + iy * nygrid * height
        bxmax, bymax = xmin + (ix+1) * nxgrid * width, ymin + (iy+1) * nygrid * height
        if crop is True:
            if bxmin > xmax or bymin > ymax:
                return None
            # elif bxmax > xmax or bymax > ymax:
            #    bxmax, bymax = min(bxmax, xmax), min(bymax, ymax)
            if bxmax > xmax:
                bxmax = xmin + int((xmax - xmin - GridProcessor.TOL_EPS) / width +1) * width
            if bymax > ymax:
                bymax = ymin + int((ymax - ymin - GridProcessor.TOL_EPS) / height +1) * height
        return [bxmin - buffx, bymin - buffy, bxmax + buffx, bymax + buffy]
        
    #/************************************************************************/
    @staticmethod
    def get_pos_location(cellsize, bbox, pos='LLc', buffer=None, yreverse=True):
        # Return the location of the ['LLc','LRc','URc','ULc','CC'] of the grid cells
        height, width = cellsize
        buffy, buffx = buffer if buffer is not None else [0,0]
        xstart, ystart, xend, yend = map(sum, zip(bbox, [-buffx, -buffy, buffx, buffy]))
        xsize, ysize = xend - xstart, yend - ystart
        # if loc == 'LLc':   pass # i.e.: do nothing
        if pos in ['LRc','URc']:
            xstart =  xstart+width
        if pos in ['ULc','URc']:
            ystart = ystart+height #if yreverse is False else ystart-height
        if pos in ['CC','centre']:
            xstart = xstart+width/2
            ystart = ystart+height/2 #if yreverse is False else ystart-height/2
        if True: # this way, we also deal with non integer (height,width)
            idrows = [ystart + i*height for i in range(int(np.ceil(ysize/height)))] 
            idcols = [xstart + i*width for i in range(int(np.ceil(xsize/width)))] 
        else:
            idrows = list(range(int(np.floor(ystart)), int(np.ceil(yend)), height))
            idcols = list(range(int(np.floor(xstart)), int(np.ceil(xend)), width))
        if yreverse is True:
            idrows.reverse()
        return idrows, idcols

    #/************************************************************************/
    @classmethod
    def build_from_pos(cls, cellsize, idrows, idcols, pos='LLc'):
        # Return all unit cells of a regular grid as a list of boundin box polygons
        height, width = cellsize
        # if loc == 'LLc':   pass # i.e.: do nothing
        if pos in ['LRc','URc']:
            idcols = map(lambda x: x-width, idcols)
        if pos in ['ULc','URc']:
            idrows = map(lambda y: y-height, idrows)
        if pos in ['CC','centre']:
            idcols = map(lambda x: x-width/2, idcols)
            idrows = map(lambda y: y-height/2, idrows)
        polygrid = []
        [polygrid.append(cls.bbox_to_polygon(x, y, x+width, y+height)) 
         for x in idcols for y in idrows] # note the order: cols then rows
        return polygrid
    
    #/************************************************************************/
    @staticmethod
    def align_pos_location(cellsize, bbox, loc, pos='LLc', maxsize=None):
        # Return a bounding box for a regular grid that will encompass the provided 
        # bounding box and whose cells with given resolution will pass through the
        # given locations 
        height, width = cellsize
        if not isinstance(bbox, (tuple,list)): # bbox is gpd.GeoSeries
            bbox = bbox.total_bounds.tolist()
        if maxsize is None:
            maxsize = max(height, width)
        ceildist = lambda b, p, size: size * np.ceil(abs(b - p)/size)
        if len(loc) == 2: # only one location (x,y) is parsed
            loc = [loc[0], loc[1], # in the following, no impact with putting the same
                   loc[0], # + ceildist(bbox[2], loc[0], width) - width,
                   loc[1], # + ceildist(bbox[3], loc[1], height) - height
                  ]
        # if loc == 'LLc':   pass # i.e.: do nothing
        if pos in ['LRc','URc']:
            loc[0], loc[2] = loc[0] - width, loc[2] - width
        if pos in ['ULc','URc']:
            loc[1], loc[3] = loc[1] - height, loc[3] - height
        if pos in ['CC','centre']:
            loc[0], loc[2] = loc[0] - width/2, loc[2] - width/2
            loc[1], loc[3] = loc[1] - height/2, loc[3] - height/2
        xmax, ymax = max(bbox[2], loc[2]+maxsize), max(bbox[3], loc[3]+maxsize)
        loc_new = [loc[0]         if loc[0]<=bbox[0]  else loc[0]-ceildist(loc[0], bbox[0], width),
                   loc[1]         if loc[1]<=bbox[1]  else loc[1]-ceildist(loc[1], bbox[1], height),
                   loc[2]+maxsize if loc[2]>=bbox[2]  else loc[2]+ceildist(loc[2], xmax, width),
                   loc[3]+maxsize if loc[3]>=bbox[3]  else loc[3]+ceildist(loc[3], ymax, height)
                  ]
        # if loc == 'LLc':   pass # i.e.: do nothing
        if pos in ['LRc','URc']:
            loc_new[0], loc_new[2] = loc_new[0] + width, loc_new[2] + width
        if pos in ['ULc','URc']:
            loc_new[1], loc_new[3] = loc_new[1] + height, loc_new[3] + height
        if pos in ['CC','centre']:
            loc_new[0], loc_new[2] = loc_new[0] + width/2, loc_new[2] + width/2
            loc_new[1], loc_new[3] = loc_new[1] + height/2, loc_new[3] + height/2
        return loc_new
    
    #/************************************************************************/
    @classmethod
    def intersection(cls, *arg, **kwargs):
        # Return the intersection of an input geometry with an arbitrary number of 
        # other geometries. The type of the returned object depends on the type of  
        # the input geomety (bounding box or polygon)
        crs = kwargs.pop('crs', None)
        rtree = kwargs.pop('rtree', False)
        def _intersection_bbox(bbox, other):
            if isinstance(other, (gpd.GeoDataFrame,gpd.GeoSeries)):
                other = other.total_bounds.tolist()
            elif isinstance(other, (geometry.Polygon,geometry.MultiPolygon)):
                other = other.bounds
            return [max(bbox[0], other[0]), max(bbox[1], other[1]),
                    min(bbox[2], other[2]), min(bbox[3], other[3])]
        def _intersection_polygon(polygon, other):
            if rtree is True:
                if isinstance(other, (gpd.GeoDataFrame,gpd.GeoSeries)):
                    sindex = other.sindex
                elif isinstance(other, (geometry.Polygon,geometry.MultiPolygon)):
                    sindex = strtree.STRtree(other)
                else:
                    sindex = None
            if isinstance(other, (tuple,list)):
                other = cls.bbox_to_polygon(*other) # cls.bbox_to_geoframe(*other, crs=crs)
            elif isinstance(other, (gpd.GeoDataFrame,gpd.GeoSeries)):
                other = other.geometry
            #elif isinstance(other, (geometry.Polygon,geometry.MultiPolygon)):
            #    other = gpd.GeoDataFrame(index=[0], geometry=[geom])
            #    if crs is not None:
            #        other.crs = crs
            if rtree is True and sindex is not None:
                try:
                    sidx = list(sindex.intersection(other))
                except:
                    sidx = sindex.query(other)
                if sidx == []:
                    return None
                try:
                    polygon = polygon.iloc[sidx].copy()
                except: 
                    polygon = polygon[sidx]
            if True:
                return polygon.intersection(other)
            else:
                return gpd.overlay(polygon, other, how='intersection')
        if geom is None and len(args)>=1: # support iterative use
            geom, args = args[0], args[1:]
        if len(args)==0:
            return geom
        return functools.reduce(_intersection_bbox if isinstance(geom, (tuple,list)) else _intersection_polygon, 
                                [geom, *args])   
    
    #/************************************************************************/
    @classmethod
    def unary_union(cls, geom, *args, **kwargs):
        # Perform repeated unary_union operations of a given geometry with an arbitrary 
        # number of other geometries. The type of the returned object depends on the type 
        # of the input geomety (bounding box or polygon)
        crs = kwargs.pop('crs', None)
        def _union_bbox(bbox, other):
            if isinstance(other, (gpd.GeoDataFrame,gpd.GeoSeries)):
                other = other.total_bounds.tolist()
            elif isinstance(other, (geometry.Polygon,geometry.MultiPolygon)):
                other = other.bounds
            return [min(bbox[0], other[0]), min(bbox[1], other[1]),
                    max(bbox[2], other[2]), max(bbox[3], other[3])]
        def _union_polygon(polygon, other):
            if isinstance(other, (tuple,list)):
                other = cls.bbox_to_polygon(*other) # cls.bbox_to_geoframe(*other, crs=crs)
            elif isinstance(other, (gpd.GeoDataFrame,gpd.GeoSeries)):
                other = other.unary_union
            #elif isinstance(other, (geometry.Polygon,geometry.MultiPolygon)):
            #    other = gpd.GeoDataFrame(index=[0], geometry=[geom])
            #    if crs is not None:
            #        other.crs = crs
            if True:
                return polygon.union(other)
            else:
                return gpd.overlay(polygon, other, how='union')
        if geom is None and len(args)>=1: # support iterative use
            geom, args = args[0], args[1:]
        if len(args)==0:
            if isinstance(geom, (gpd.GeoDataFrame,gpd.GeoSeries)):
                return geom.unary_union
            elif isinstance(geom, (geometry.Polygon,geometry.MultiPolygon)):
                return ops.unary_union(geom)
            else:
                return geom
        return functools.reduce(_union_bbox if isinstance(geom, (tuple,list)) else _union_polygon, 
                                [geom, *args])    


#%% Main for binary usage
