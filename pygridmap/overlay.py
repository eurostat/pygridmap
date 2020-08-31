#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
.. _gridding

Grid making operations.

**Dependencies**

*require*:      :mod:`numpy`, :mod:`pandas`, :mod:`geopandas`, :mod:`shapely`

*optional*:     :mod:`multiprocessing`

*call*:         :mod:`pygridmap.base`, :mod:`pygridmap.gridding`

**Contents**
"""

# *credits*:      `gjacopo <jacopo.grazzini@ec.europa.eu>`_ 
# *since*:        Mon 6 2020


#%% Settings     

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

from pygridmap.base import FrameProcessor, GridProcessor#analysis:ignore
from pygridmap.gridding import GridMaker
from pygridmap.base import NPROCESSES, NCPUS

#%% Core functions/classes

#==============================================================================
# Class GridOverlay
#==============================================================================

class GridOverlay(GridProcessor):

    HOWS = ['intersection', 'union']
    RULES = ['sum', 'max', 'min', 'list'] # 'gsum', 'hsum'
    SORTS = ['rc', 'cr']
    
    try:
        COL_TILE = GridMaker.COL_TILE
        COL_INTERSECTS, COL_WITHIN = GridMaker.COL_INTERSECTS, GridMaker.COL_WITHIN
        COL_X, COL_Y = GridMaker.COL_X, GridMaker.COL_Y
    except:
        COL_TILE = '__tile__'
        COL_INTERSECTS, COL_WITHIN = '__intersects__', '__within__'
        COL_X, COL_Y = '__x__', '__y__'
        
    COL_GRIDX, COL_POLIDX = '__gridx__', '__polidx__'
    COL_AREA, COL_PCT = '__area__', '__area_pct__'
    
	NAN = np.nan if True else 0 #?
 
    #/************************************************************************/
    def __init__(self, **kwargs):
        self.__mode, self.__how = None, None
        self.__keep_geom_type, self.__preserve_polygon = True, False
        super(GridOverlay,self).__init__(**kwargs)
        self.mode = 'prll' # kwargs.pop('mode', 'prll') # ignored
        self.how = kwargs.pop('how', 'intersection') # self.HOWS[0]
        self.keep_geom_type = kwargs.pop('keep_geom_type', True)
        self.preserve_polygon = kwargs.pop('preserve_polygon', False)
        
    #/************************************************************************/
    @property
    def mode(self):
        return self.__mode
    @mode.setter
    def mode(self, mode):
        try:
            assert isinstance(mode, str)
        except: raise TypeError("Wrong format for overlay processing mode")
        try:
            assert (mode in self.MODES)
        except: raise IOError("Wrong value for overlay processing mode")
        self.__mode = mode

    #/************************************************************************/
    @property
    def how(self):
        return self.__how
    @how.setter
    def how(self, how):
        try:
            assert isinstance(how, str)
        except: raise TypeError("Wrong format for overlay geometrical operation")
        try:
            assert (how in self.HOWS)
        except: raise IOError("Wrong value for overlay geometrical operation")
        self.__how = how
        
    #/************************************************************************/
    @property
    def preserve_polygon(self):
        return self.__preserve_polygon
    @preserve_polygon.setter
    def preserve_polygon(self, pres):        
        try:
            assert isinstance(pres, bool)
        except: raise TypeError("Wrong format for preserve_polygon parameter")
        self.__preserve_polygon = pres
        
    #/************************************************************************/
    @property
    def keep_geom_type(self):
        return self.__keep_geom_type
    @keep_geom_type.setter
    def keep_geom_type(self, keep):        
        try:
            assert isinstance(keep, bool)
        except: raise TypeError("Wrong format for keep_geom_type parameter")
        self.__keep_geom_type = keep
        
    #/************************************************************************/
    @property
    def asc(self):
        return self.__asc
    @asc.setter
    def asc(self, asc):        
        try:
            assert (isinstance(asc, bool)                                                   \
                or (isinstance(asc, (tuple,list)) and all([isinstance(a,bool) for a in asc])))
        except: raise TypeError("Wrong format for ascending parameter")
        self.__asc = [True,True] if asc is True else asc
    
    #/************************************************************************/
    @classmethod
    def trim_grid(cls, nxtiles, iy, ix, grid, polyarea):
        tile = grid[grid[cls.COL_TILE] == iy*nxtiles+ix].copy()
        if tile.iloc[0, tile.columns.get_loc(cls.COL_INTERSECTS)] == 0:
            return None
        elif tile.iloc[0, tile.columns.get_loc(cls.COL_WITHIN)] == 1:
            return tile
        # note on the use of spatial indexing in sjoin: "if both dataframes have a spatial index, 
        # the right_df’s index will be used preferentially"
        poly = gpd.sjoin(tile, polyarea, op = 'intersects', how = 'left')
        return poly.drop(columns=list(set(tile.columns).difference(set(poly.columns))), errors = 'ignore')
    
    #/************************************************************************/
    @classmethod
    def crop_grid_from_tile(cls, idx, grid, gridbbox, cellsize, tilesize,
                            sorted_, preserve_tile=False):
        iy, ix = idx[:2]
        tilebbox = None
        # create a single tile object for clipping
        if tilesize is True: 
            tilesize = cls.COL_TILE
        if isinstance(tilesize,str):
            if tilesize not in grid.columns:
                raise IOError('tile column name not found')
            index = grid[grid[tilesize] == iy].index
        elif cellsize is not None:
            xmin, ymin, xmax, ymax = cls.get_tile_bbox(idx, cellsize, tilesize, gridbbox, False) # no crop
            tilebbox = cls.bbox_to_polygon (xmin, ymin, xmax, ymax)
            # retrieve the indexes of the grid cells within the considered tile
            if cls.COL_X in grid.columns and cls.COL_Y in grid.columns:
                index = grid[(grid[cls.COL_X] >= xmin) & (grid[cls.COL_X]  <= xmax) \
                            & (grid[cls.COL_Y] >= ymin)  & (grid[cls.COL_Y] <= ymax)].index
            elif sorted_ not in (False, None):
                nygrid, nxgrid = tilesize
                nrows, ncols = cls.get_grid_shape(cellsize, gridbbox)
                index = [(iy*nygrid + j)*ncols + (ix*nxgrid + i) if sorted_ == 'rc'    \
                         else (iy*nygrid + j) + (ix*nxgrid + i)*nrows               \
                          for j in range(nygrid) for i in range(nxgrid)]
                # index = [(ix*nxgrid + i)*nrows + (iy*nygrid + j) for j in range(nygrid) for i in range(nxgrid)]
                index = [ind for ind in index if ind<=nrows*ncols - 1]
            elif True:
                bbox = tilebbox.buffer(-min(cellsize)*1e-2).bounds
                index = list(grid.sindex.intersection(bbox))
            else:
                geometry = tilebbox.buffer(0).geometry
                clip_bbox = gpd.GeoDataFrame(geometry = [geometry], 
                                             crs = grid.crs)
                # buffer to solve inconsistent geometries
                return gpd.clip(grid, clip_bbox, keep_geom_type=True).buffer(0) 
        if index is None or index is []:
            return tilebbox if preserve_tile is True else None
        crop_grid = (grid
                     .iloc[index] # .copy()
                     .reset_index() # index reset at the level of the processing tile
                     .rename(columns={'index': cls.COL_GRIDX})
                    )
        return crop_grid
    
    #/************************************************************************/
    @classmethod
    def geometry_overlay(cls, gridarea1, polygon2, how_overlay,
                         keep_geom_type = True, preserve_polygon = False): 
        # Perform set/geometric intersection/union operations (how_overlay) of a grid (gridarea1)
        # and a model/polygon layer (polygon2) at tile level
        # 1. build the dissolved version of the polyarea2 (e.g. "merge" the grid geometries into a subgrid)
        # gridarea1.reset_index(inplace=True)
        union1 = gridarea1.geometry.unary_union # it is a tile in the case a subgrid is in the interior
        # or (faster?):
        # union1 = GeoProcesor.bbox_to_geoframe(*gridarea1.total_bounds, crs = gridarea1.crs)
        # 2. create a bounding box for the initial intersection
        bbox1 = union1.bounds # one/unary geometry only, same as polyarea1.total_bounds
        # 3. retrieve the spatial index - note tht it can be created outside when calling it the first time
        sindex2 = polygon2.sindex
        # 4. get a list of id's for all polygons in polygon2 that overlap/intersect the polyarea1 bounding box 
        sidx2 = list(sindex2.intersection(bbox1))
        # checking whether the geometries actually overlap/intersect, otherwise return the polyarea1
        if sidx2 == []:
            return polyarea1 if preserve_polygon is True else None
        # 5. clip / subset the polygons in polygon2 that intersect
        if cls.COL_POLIDX not in polygon2.columns:
            clipped2 = (polygon2.iloc[sidx2]
                        .reset_index() # index reset at the level of the processing tile
                        .rename(columns={'index': cls.COL_POLIDX})
                       ) 
        else:
            clipped2 = polygon2.iloc[sidx2].copy() # the presence of 'copy' will prevent a SettingWithCopyWarning
        if cls.COL_AREA not in polygon2.columns:
            clipped2[cls.COL_AREA] = clipped2.area
        if True:
            clipped2['geometry'] = clipped2.intersection(union1)
        else: # geopandas >0 .8
            clipped2 = gpd.clip(clipped2, 
                                union1,
                                keep_geom_type = keep_geom_type
                               )
        # filter out clipped linestrings and points with null area
        clipped2 = clipped2[clipped2.area>0] # note: clipped2[cls.COL_AREA] != clipped2.area
        # clean by getting rid of empty geometries, possibly returning the subgrid when empty
        clipped2 = clipped2[~clipped2.is_empty] 
        if clipped2.is_empty.all():
            return gridarea1 if preserve_polygon is True else None
        # 6. when not empty, retrieve the (geometrical) union/intersection of the subgrid with the 
        # clipped polygons
        # see overlay implementation at https://github.com/geopandas/geopandas/blob/master/geopandas/tools/overlay.py
        return gpd.overlay(gridarea1, clipped2, how_overlay, 
                           keep_geom_type = keep_geom_type, # True,
                           make_valid = True)
        # return gpd.overlay(clipped2, gridarea1, how_overlay, keep_geom_type = True, make_valid = True)
    
    #/************************************************************************/
    @classmethod
    def cover_overlay(cls, overlay, pct_thres = None):  
        lookup = (overlay[~np.isnan(overlay[cls.COL_PCT])]
                 .groupby(cls.COL_GRIDX, as_index = False)
                 .aggregate({cls.COL_POLIDX: list})
                 .set_index(cls.COL_GRIDX, drop = True)
                )
        overlay[cls.COL_POLIDX] = (overlay[cls.COL_GRIDX]
                                   .map(lookup.to_dict()[cls.COL_POLIDX])
                                   .fillna(overlay[cls.COL_POLIDX])
                                  )
        return overlay
    
    #/************************************************************************/
    @classmethod
    def area_overlay(cls, overlay, pct_thres = None):  
        # the area information in overlay comes from the polygons, we store it in a variable
        poly_area = overlay[cls.COL_AREA].copy() if cls.COL_AREA in overlay.columns else 0
        # null area, will ensure a division /0 and Inf values
        # we recompute the area of the overlayed geometries
        overlay[cls.COL_AREA] = overlay.area # we udpate the info
        overlay[cls.COL_PCT] = overlay[cls.COL_AREA].div(poly_area)
        overlay.loc[~np.isfinite(overlay[cls.COL_PCT]), cls.COL_PCT] = cls.NAN
        return overlay
    
    #/************************************************************************/
    @classmethod
    def attribute_overlay(cls, polyarea1, overlay2, rule, columns):  
        if rule in ['min', 'max']:
            overlay2 = (overlay2
                        .sort_values(cls.COL_PCT, ascending = True if rule == 'min' else False)
                      )
        else:
            overlay2.update(
                overlay2
                .iloc[:, [overlay2.columns.get_loc(c) for c in columns]]
                .mul(overlay2[cls.COL_PCT], axis = 0)
            )
            overlay2.update(
                overlay2
                .groupby(cls.COL_GRIDX)[columns]
                .transform(rule) # rule = 'sum'
            ) 
        # in all cases, keep one representative grid cell only, the first:
        overlay2 = (overlay2 
                    .groupby(cls.COL_GRIDX)
                    .first() 
                    .reset_index()
                  )
        # note: another approach for ['min', 'max']
        #    index = (overlay2
        #             .sort_values(cls.COL_PCT, ascending = True if rule=='min' else False)
        #             .drop_duplicates(cls.COL_GRIDX)
        #             .index
        #            )
        #    overlay2 = overlay2.iloc[index]
        # update the resulting selected geometries
        return (polyarea1
                .merge(overlay2
                       .drop(columns = [cls.COL_PCT, cls.COL_AREA, 'geometry'], errors = 'ignore'), 
                       on = cls.COL_GRIDX
                      ) 
               )

    #/************************************************************************/
    @classmethod
    def process_split_tile(cls, polygons, subgrid, 
                           area, cover, rule, columns, how_overlay,
                           keep_geom_type, preserve_polygon):
        # define the subgrid over the specific tile considered for processing
        if columns not in (None,[]): 
            # and (intercols:=set(columns).intersection(set(subgrid.columns))) != set()):
            subgrid.drop(columns = columns, inplace = True, errors = 'ignore')
        # define the geometric overlay (e.g., how_overlay=intersection or union) of both the subgrid 
        # and the polygon representation
        overlay = cls.geometry_overlay(subgrid, polygons, how_overlay, 
                                       keep_geom_type, preserve_polygon)
        if overlay is None or overlay.is_empty.all():
            return None
        # udpate the attributes available in the overlay geometries to take into account the proportional 
        # area coverage
        if area is True or rule is not None:
            overlay = cls.area_overlay(overlay, pct_thres = None)
        if cover is True:
            overlay = cls.cover_overlay(overlay, pct_thres = None)
        if rule is None or columns in (None,[]):
            return overlay
        overlay.drop(columns = [cls.COL_INTERSECTS, cls.COL_WITHIN, cls.COL_TILE, cls.COL_X, cls.COL_Y], 
                     inplace = True, errors = 'ignore')
        overlay = cls.attribute_overlay(subgrid, 
                                        overlay, 
                                        rule,
                                        columns)
        return overlay
    
    #/************************************************************************/
    @classmethod
    def process_tile(cls, idx, polygons, grid, gridbbox, cellsize, tilesize,
                     area, cover, rule, columns, how_overlay, 
                     keep_geom_type, sorted_, preserve_polygon):
        # define the subgrid over the specific tile considered for processing
        subgrid = cls.crop_grid_from_tile(idx, grid, gridbbox,
                                          cellsize, tilesize,  
                                          sorted_) 
        if subgrid is None or subgrid.is_empty.all():
            return None        
        # define the overlay (e.g., how_overlay=intersection or union) geometries of both the polygons
        # and the subgrid
        return cls.process_split_tile(polygons, subgrid, 
                                      area, cover, columns, rule, how_overlay, 
                                      keep_geom_type, preserve_polygon)
      
    #/************************************************************************/
    def __call__(self, polygons, grid, 
                 area = True, cover = False, rule = None, columns = None, 
                 memory_split = True, drop = None):  
        # check area flag
        try:
            assert isinstance(area, bool)
        except: raise TypeError("Wrong format for area flag")
        # check cover flag
        try:
            assert isinstance(cover, bool)
        except: raise TypeError("Wrong format for cover flag")
        # check rule
        try:
            assert (rule is None or isinstance(rule, str))
        except: raise TypeError("Wrong format for feature assignment rule")
        try:
            assert (rule is None or rule in self.RULES)
        except: raise IOError("Wrong value for feature assignment rule")
        # check columns
        try:
            assert (columns is None or isinstance(columns, str)
                    or (isinstance(columns, (tuple,list)) 
                        and all([isinstance(a, str) for a in columns])))
        except: raise TypeError("Wrong format for columns variable")   
        if isinstance(columns,str): 
            columns = [columns,]
        try:
            assert (columns is None or 
                    set(columns).difference(set(polygons.columns)) == set())
        except: raise IOError("Wrong values for columns variable") 
        # check memory_split flag
        try:
            assert isinstance(memory_split, bool)
        except: raise TypeError("Wrong format for memory_split flag")
        # check drop parameter
        __drops = [self.COL_TILE, self.COL_X, self.COL_Y, 
                   self.COL_INTERSECTS, self.COL_WITHIN,
                   self.COL_AREA, self.COL_PCT,
                   self.COL_POLIDX, self.COL_GRIDX
                  ]
        try:
            assert (drop is None or isinstance(drop, bool)                                  \
                or (isinstance(drop,str) and drop in __drops)                               \
                or (isinstance(drop,(tuple,list)) and all([d in __drops for d in drop])))
        except: raise TypeError("Wrong format for drop parameter")
        if drop is None:              drop = False
        elif isinstance(drop, str):   drop = [drop,]
        # settings
        cellsize, tilesize = self.cell, self.tile or 1
        # defining the grid dimensions
        gridbbox = self.get_bbox(grid, polygons, gridbbox=True)
        if isinstance(tilesize, str):
            iytiles, ixtiles = list(set(grid[tilesize].values)), [0]
            nytiles, nxtiles = len(iytiles), 1
        elif np.isscalar(tilesize): # or isinstance(tilesize, (tuple,list)) and len(tilesize)==1
            tileshape = self.set_tile_shape(tilesize)
            if cellsize is not None:
                tilesize = self.get_tile_size(cellsize, tileshape, gridbbox)
            else:
                tilesize = None
        elif isinstance(tilesize, (tuple,list)):
            if cellsize is not None:
                tileshape = self.get_tile_shape(cellsize, tilesize, gridbbox)
                # tilesize: unchanged
            else:
                tileshape, tilesize = tilesize, None
        if not isinstance(tilesize, str):
            nytiles, nxtiles = tileshape
            iytiles, ixtiles = list(range(nytiles)), list(range(nxtiles))
        # start processing        
        pool = mp.Pool(processes = self.cores)
        if cover is True:
            # polygons = (polygons
            #             .reset_index()
            #             .set_index('index', drop = False)
            #             .rename(columns={'index': self.COL_POLIDX})
            #            ) 
            polygons[self.COL_POLIDX] = polygons.index.copy()
        if memory_split is True:
            overlay_tile = []
            for ix in ixtiles:
                for iy in iytiles:
                    subgrid = self.crop_grid_from_tile([iy, ix, nytiles, nxtiles], 
                                                       grid, gridbbox, cellsize, tilesize,  
                                                       self.sorted) 
                    if subgrid is None or subgrid.is_empty.all():
                        continue        
                    overlay_tile.append(pool.apply_async(self.process_split_tile,
                                                         args = (polygons, subgrid, 
                                                                 area, cover, rule, columns, 
                                                                 self.how, self.keep_geom_type, 
                                                                 self.preserve_polygon
                                                                ))
                                       )
        else:
            overlay_tile = [pool.apply_async(self.process_tile,
                                             args = ([iy, ix, nytiles, nxtiles], 
                                                     polygons, grid, gridbbox, 
                                                     cellsize, tilesize, 
                                                     area, cover, rule, columns,
                                                     self.how, self.keep_geom_type,
                                                     self.sorted, 
                                                     self.preserve_polygon
                                                    ))    
                            for iy in iytiles for ix in ixtiles]       
        pool.close()
        pool.join()
        overlayed = pd.concat([t.get() for t in overlay_tile if t is not None], 
                         axis=0, ignore_index=True)
        if drop != False: 
            if drop is True:
                drop = [self.COL_TILE, self.COL_INTERSECTS, self.COL_WITHIN, self.COL_POLIDX] #__drops
                drop.extend(polygons.columns.tolist())
            if area is True: # we keep the cover indexing column
                drop = set(drop).difference({self.COL_AREA, self.COL_PCT})
            if cover is True: # we keep the cover indexing column
                drop = set(drop).difference({self.COL_POLIDX})
            if rule is not None: # we keep still keep those, otherwise there was no point!
                drop = set(drop).difference(set(columns))
            drop = set(drop).difference({self.COL_GRIDX, 'geometry'}) # don't drop those two!!!
            overlayed.drop(columns = list(drop), axis = 1, inplace = True, errors = 'ignore')
        return gpd.GeoDataFrame(overlayed, crs = grid.crs)


#%% Main for binary usage		