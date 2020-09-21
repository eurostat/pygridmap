pygridmap
=========

Geoprocessing tools for efficient geometric mapping and set operations over geospatial grids.
---

The `pygridmap` package enable you to perform some basic geometric/set operations over vector datasets, including regular grid maker (_i.e._, rasterisation) and grid overlay (_i.e._ intersection and union). The code takes advantage of multiprocessor capabilities for efficient tile-based processing. 

**Description**

The package `pygridmap` supports the following methods/algorithms: 
* basic [**vector data handling**](https://saylordotorg.github.io/text_essentials-of-geographic-information-systems/s11-geospatial-analysis-i-vector-o.html) over geospatial grid (_e.g._, bounding box manipulations, locations definition, _etc_...) with the [`GridProcessor`](pygridmap/base.py) class,
<!-- ![bounding boxes](docs/bbox_manipulation.png)-->
<table align="center"> 
        <header> <td align="centre">Regular gridding with different reference locations</td></header> 
        <tr> <td align="center" width="800px"> <img src="docs/bbox_manipulation.png"></img></td></tr> 
</table>

* user-defined [**rasterisation**](https://en.wikipedia.org/wiki/Rasterisation) of polygon regions (_e.g._, vector boundaries) into regular grids of various resolutions (_i.e._, dimension of one-size square grid cells) with the [`GridMaker`](pygridmap/gridding.py) class (see also `Java` implementation [GridMaker](https://github.com/eurostat/GridMaker)),
<!-- ![bounding boxes](docs/BE_interior_gridding.png)-->
<table align="center">
        <header> <td align="centre" colspan=2>Regular gridding striclty covering a spatial domain (left) or its interior only (right)</td></header> 
        <tr> <td align="center" width="300px"> <img src="docs/BE_gridding.png"></img></td>
        <td align="center" width="300px"> <img src="docs/BE_interior_gridding.png"></img></td>
        </tr>
</table>

* customised [**vector overlay operations**](https://docs.qgis.org/3.10/en/docs/user_manual/processing_algs/qgis/vectoroverlay.html) (_e.g._, intersection, union, overlay/overlap and merging) between any polygon layer and a geospatial grid ("regular" or not) with the [`Gridoverlay`](pygridmap/overlay.py) class, including weighted areal interpolation.
<!-- ![bounding boxes](docs/BE_overlay.png)-->
<table align="center">
        <header> <td align="centre">Intersection of a regular grid with an arbitrary vector layer</td></header> 
        <tr> <td align="center" width="320px"> <img src="docs/BE_overlay.png"></img></td></tr> 
</table>

**Quick launch**

* **Run (and update to your convenience) the test/example notebooks from the [tests/](tests) folder in [`binder`](https://mybinder.org/)**. We provide the interactive environments with already installed packages to query and access _Eurostat_ database for notebook resources below (current build with commit [???](https://github.com/eurostat/pygridmap/commit/???)): <!-- generate new badges: https://mybinder.readthedocs.io/en/latest/howto/badges.html -->
[![badge](https://img.shields.io/badge/Eurostat%20code-binder-579ACA.svg?logo=data:image/png;base64,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)](https://mybinder.org/v2/gh/eurostat/pygridmap/) <!-- [![Binder](https://mybinder.org/badge_logo.svg)](http://mybinder.org/v2/gh/eurostat/statistics-coded/master)--> <!-- [![Binder](https://mybinder.org/badge_logo.svg)](http://mybinder.org/v2/gh/eurostat/statistics-coded/master?filepath=index.ipynb)-->

* **Run the notebooks in [`Google colab`](https://colab.research.google.com/)** separately, namely:
    * [base.ipynb](https://nbviewer.jupyter.org/github/eurostat/pygridmap/blob/master/tests/base.ipynb) for basic data handling operations [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/eurostat/pygridmap/blob/master/tests/base.ipynb),
    * [gridding.ipynb](https://nbviewer.jupyter.org/github/eurostat/pygridmap/blob/master/tests/gridding.ipynb) for regular grid making [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/eurostat/pygridmap/blob/master/tests/gridding.ipynb),
    * [overlay.ipynb](https://nbviewer.jupyter.org/github/eurostat/pygridmap/blob/master/tests/overlay.ipynb) for overlay/overlap, merging, intersection and union operations [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/eurostat/pygridmap/blob/master/tests/overlay.ipynb).

**Quick install**

The `pygridmap` package can be installed using `pip`:
```python
>>> pip install pygridmap
```
If you wish to use the latest available version from _github_:
```python
>>> pip install git+https://github.com/eurostat/pygridmap.git
```

**Usage**

*In your `Python` script*

The following classes are available: 
```python
>>> from pygridmap.base import FrameProcessor, GridProcessor
>>> from pygridmap.gridding import GridMaker
>>> from pygridmap.overlay import GridOverlay
```

*Through the `bash` command*

**Notes**

*Processing* 

The implementation of the methods above adopts (customised, using the [`multiprocessing`](https://docs.python.org/3/library/multiprocessing.html) module) multitprocessor tiling processing together with (native) vector processing whenever possible in order to take advantage of multiprocessor compute capabilities. Note  that the [`FrameProcessor`](pygridmap/base.py) class can be used to run embarassing parallel calculations on dataframe rows.
<!-- ![tile processing](docs/BE_tile_processing.png)-->
<table align="center">
        <header> <td align="centre">Example of block-processing tiles</td></header> 
        <tr> <td align="center" width="320px"> <img src="docs/BE_tile_processing.png"></img></td> </tr> 
</table>

Geometric operations are making an extensive use of the geometric/geospatial `Python` module [`geopandas`](https://geopandas.org/) since this module supports essential features (most of them derived from the [`shapely`](https://shapely.readthedocs.io/en/latest/manual.html) [`fiona`](https://fiona.readthedocs.io/en/latest/manual.html) libraries), like:
  * basic [set](https://geopandas.org/set_operations.html) and [geometric](https://geopandas.org/geometric_manipulations.html) operations, and
  * [spatial indexing](https://geopandas.org/mergingdata.html?highlight=spatial%20index) through the [R-Tree algorithm](https://automating-gis-processes.github.io/site/notebooks/L3/spatial_index.html).
  
These implementations are probably not optimal since they are wrappers to the [`GEOS`](https://trac.osgeo.org/geos/) library, itself a port of the Java Topology Suite ([`JTS`](https://projects.eclipse.org/projects/locationtech.jts)), similarly to what is used in QGis.

*Algorithm*

On areal interpolation, the package `tobler` that is distributed with the [`pysal`](https://github.com/pysal) library enables to perform more generic operations (*i.e.*, not only considering grid as a target layer). Whenever one needs to  project/interpolate an extensive feature `attribute` from a `source` vector layer onto a `target` vector layer (*e.g.*, a regular grid), the following command:
```python
>>> from pygridmap import apps
>>> estimate = apps.area_interpolate(source, target, attribute) 
```
is equivalent to running:
```python
>>> from tobler import area_weighted
 >>> estimate = area_weighted.area_interpolate(source, target, extensive_variables = attribute)
```
<table align="center">
        <header> <td align="centre">Areal interpolation of an extensive variable (left) from coarse to fine using <a href="https://github.com/pysal/tobler/blob/master/notebooks/02_areal_interpolation_example.ipynb"><code>tobler</code> code</a> (left) and <code>pygridmap</code> algorithm (right)</td></header> 
        <tr> <td align="center" width="900px"> <img src="docs/overlay_tobler.png"></img></td> </tr> 
</table>

With respect to the latter implementation, the current algorithm supports a tile-based multicore approach for the interpolation (through the setting of the `tile` parameter). The `tobler` algorithm will however help you "project" intensive variables as well.

**About**

<table align="center">
    <tr> <td align="left"><i>documentation</i></td> <td align="left"><b>in construction</b></td>  </tr> 
    <tr> <td align="left"><i>status</i></td> <td align="left">since 2020 &ndash; <b>ongoing</b></td></tr> 
    <tr> <td align="left"><i>contributors</i></td> 
    <td align="left" valign="middle">
<a href="https://github.com/gjacopo"><img src="https://github.com/gjacopo.png" width="40"></a>
</td> </tr> 
    <tr> <td align="left"><i>license</i></td> <td align="left"><a href="https://joinup.ec.europa.eu/sites/default/files/eupl1.1.-licence-en_0.pdfEUPL">EUPL</a> </td> </tr> 
</table>

**<a name="Requirements"></a>Requirements**

* Data handling: [`numpy`](https://numpy.org/) and [`pandas`](http://pandas.pydata.org).
* Geospatial data processing: [`geopandas`](http://geopandas.org), [`fiona`](https://fiona.readthedocs.io/en/latest/manual.html), [shapely](https://pypi.org/project/Shapely/), [`RTree`](https://toblerity.org/rtree/), [`geopy`](https://github.com/geopy/geopy) and [`pyproj`](https://github.com/pyproj4/pyproj).
* Map visualisations: [`ipyleaflet`](https://github.com/jupyter-widgets/ipyleaflet) or [`folium`](https://github.com/python-visualization/folium).

**<a name="Resources"></a>Other resources**

* Project [python-geospatial](https://github.com/giswqs/python-geospatial), a collection of `Python` packages for geospatial analysis, and earthdatascience, a repository of `Python` [tutorials](https://www.earthdatascience.org/tutorials/python/) and [workshop](https://www.earthdatascience.org/workshops/gis-open-source-python/).
* Library [`pysal`](https://github.com/pysal/pysal) for geospatial data science with an emphasis on geospatial vector data.
* Package [`tobler`](https://github.com/pysal/tobler) for an alternative implementation of areal interpolation (module [`area_weighted`](https://github.com/pysal/tobler/tree/master/tobler/area_weighted)), as well as dasymetric mapping, and change of support.
* Software [`GridMaker`](https://github.com/eurostat/GridMaker) for alternative grid maker. 

**<a name="References"></a>References**

* Gimond, M. (2019): [**Intro to GIS and Spatial Analysis**](https://mgimond.github.io/Spatial/index.html).
<!-- * Lovelace R., Nowosad J. and Muenchow J. (2019): [**Geocomputation with R**](https://geocompr.robinlovelace.net/), _Chapman & Hall/CRC_. -->
* de Smith M.J., Goodchild M.F. and Longley P.A. (2018): [**Geospatial Analysis: A Comprehensive Guide to Principles, Techniques and Software Tools**](https://www.spatialanalysisonline.com/HTML/index.html), _The Winchelsea Press_. 
* Garrard C. (2016): **Geoprocessing with Python**, _Manning Publications_.
* Westra E. (2016): [**Python Geospatial Development**](https://www.programmer-books.com/wp-content/uploads/2019/04/Python-Geospatial-Development-3rd-Edition.pdf), _PACKT Publishing_.
* Bahgat K. (2016): **Python Geospatial Development Essentials**, _PACKT Publishing_.
* Lawed J. (2015): [**QGIS Python Programming CookBook**](https://www.programmer-books.com/wp-content/uploads/2019/05/QGIS-Python-Programming-Cookbook.pdf), _PACKT Publishing_.
* Kresse W. and Danko D.M., eds. (2012): **Handbook of Geographic Information**, _Springer_, doi: [10.1007/978-3-540-72680-7](https://doi.org/10.1007/978-3-540-72680-7).
*  Campbell J. and Shin M. (2011): [**Essentials of Geographic Information Systems**](https://saylordotorg.github.io/text_essentials-of-geographic-information-systems/index.html) ([pdf](https://resources.saylor.org/wwwresources/archived/site/textbooks/Essentials%20of%20Geographic%20Information%20Systems.pdf)), _Saylor Academy Open Textbooks_.
