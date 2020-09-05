#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
.. __init__

Initialisation module of `pygridmap` package.

**Contents**
"""

# *credits*:      `gjacopo <jacopo.grazzini@ec.europa.eu>`_ 
# *since*:        Jun 2020

__THISFILE          = __file__ # useles...
__THISDIR           = osp.dirname(__THISFILE)

__all__             = ['base', 'gridding', 'overlay']#analysis:ignore

#PROJS               = ["EPSG:4326", "EPSG:3035"]
#"""Available projections.
#""" See pyproj
DEFPROJ             = "EPSG:4326" # "EPSG:3035"
"""Default projection in calculations.
"""
