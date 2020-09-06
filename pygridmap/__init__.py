#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
.. __init__

Initialisation module of `pygridmap` package.

**Contents**
"""

# *credits*:      `gjacopo <jacopo.grazzini@ec.europa.eu>`_ 
# *since*:        June 2020

# 88""Yb Yb  dP  dP""b8 88""Yb 88 8888b.  8b    d8    db    88""Yb 
# 88__dP  YbdP  dP   `" 88__dP 88  8I  Yb 88b  d88   dPYb   88__dP 
# 88"""    8P   Yb  "88 88"Yb  88  8I  dY 88YbdP88  dP__Yb  88"""  
# 88      dP     YboodP 88  Yb 88 8888Y"  88 YY 88 dP""""Yb 88     
                                                                
#%% Settings

from os import path as osp

#%% Global vars

#__THISFILE          = __file__ # useles...
#__THISDIR           = osp.dirname(__THISFILE)

PACKNAME            = 'pygridmap' # this package...
"""Name of this package.
"""

#PACKPATH            = getattr(sys.modules[PACKNAME], '__path__')[0]
PACKPATH            = osp.dirname(__file__) 
"""Path to this package.
"""

__all__             = [ ] # 'base', 'gridding', 'overlay']#analysis:ignore
__all__.extend(['__version__']) # , '__start__'
