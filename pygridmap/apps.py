#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
.. _apps 

Main applications: grid maker and areal interpolation.

**Dependencies**

*call*:         :mod:`pygridmap.base`, :mod:`pygridmap.gridding`, :mod:`pygridmap.overlay`

**Contents**
"""

# *credits*:      `gjacopo <jacopo.grazzini@ec.europa.eu>`_ 
# *since*:        Jun 2020


#%% Settings     

from pygridmap.base import GridProcessor, NPROCESSES, NCPUS
from pygridmap.gridding import GridMaker, DEFPROJ
from pygridmap.overlay import GridOverlay


#%% Classes / functions     

#==============================================================================
# Method grid_maker
#==============================================================================

def grid_maker(source, cell, bbox = None, crs = None, tile = None, cores = None, interior = False):
    """Gridding of a vector layer onto a regular (square) grid
    
        >>> grid = grid_maker(source, cell, bounding_box = None, crs = None, 
                              tile = None, cores = None, interior = False)
    
    Description
    -----------
    
    Given a vector layer `source` and a `cell` resolution, a regular grid with (identical) unit cells of
    dimension `cell` is created so as to cover the spatial coverage of the source data. 
    
    Because the gridding is not unique, the domain is adjusted throught the introduction of the bounding 
    box `bbox` which can be used to parse the start location (LLc) of the grid. Otherwise, the bounding box
    of the `source` data is used.
    """
    cores = cores or NCPUS
    mode = 'prll' # default for sequential as well
    tile = tile or NPROCESSES
    trim, buffer, crop = True, False, False
    crs = source.crs if crs in (True,None) else (DEFPROJ if crs is False else crs)
    bbox = bbox or source.total_bounds.tolist()

    proc = GridMaker(cores = cores, mode = mode, cell = cell, tile = tile, buffer = buffer)
    return proc(bbox, mask = source, crs = crs, interior = interior, 
                trim = trim, crop = crop)


#==============================================================================
# Method area_interpolate
#==============================================================================

def area_interpolate(source, target, extensive_variables, 
                     cell = None, tile = None, cores = None, memory_split = False):
    """Areal interpolation of extensive variable.
    
        >>> inter = area_interpolate(source, target, attribute, cell = None, 
                                    tile = None, cores = None, memory_split = False)
    
    Description
    -----------
    
    Given a vector layer `source` with extensive variable(s)/field(s) `attribute` (*e.g.*, a population
    variable represented at a given scale), running:
    
        >>> from pygridmap import apps
        >>> estimate = gridding.area_interpolate(source, target, attribute) 
        
    is then equivalent to:
    
        >>> from tobler import area_weighted
        >>> estimate = area_weighted.area_interpolate(source, target, extensive_variables = attribute)
        
    where the target layer `target` is typically a regular grid used to "rasterise" the variable(s) 
    `attribute` from a coarser polygonal resolution (vector representation available in `source`) to a 
    finer square resolution (unit cells of `target` grid). 
    
    Note
    ----
    
    The function actually accepts any polygonal representation as a target layer, *i.e.* `target` does not 
    have to be necessarly a grid.
    """
    cores = cores or NCPUS
    tile = tile or NPROCESSES
    rule = 'sum'
    area, cover, drop = True, True, True
    
    proc = GridOverlay(mode = mode, how = 'intersection', cell = cell, cores = cores, tile = tile)
    try:
        assert memory_split is None
        proc.memory_split = False
        return proc(source, target, rule = rule, columns = extensive_variables, 
                    area = area, cover = cover, drop = drop)
    except: # at this stage: either memory_split was None and the proc crashed, or it was not None
        proc.memory_split = memory_split or True # if it is None, then try with True this time
        return proc(source, target, rule = rule, columns = extensive_variables, 
                    area = area, cover = cover, drop = drop)
    