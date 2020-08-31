pygridmap
=========

Geometric grid and set mapping operations.
---


**About**

The `pygridmap` package enable you to perform some basic geometric/set operations over vector datasets, including regular grid maker (_i.e._, rasterisation) and grid overlay (_i.e._ intersection and union). The code takes advantage of multiprocessor capabilities for efficient tile-based processing. 

<table align="center">
    <tr> <td align="left"><i>documentation</i></td> <td align="left"><b>in construction</b></td>  </tr> 
    <tr> <td align="left"><i>status</i></td> <td align="left">since 2020 &ndash; <b>ongoing</b></td></tr> 
    <tr> <td align="left"><i>contributors</i></td> 
    <td align="left" valign="middle">
<a href="https://github.com/gjacopo"><img src="https://github.com/gjacopo.png" width="40"></a>
</td> </tr> 
    <tr> <td align="left"><i>license</i></td> <td align="left"><a href="https://joinup.ec.europa.eu/sites/default/files/eupl1.1.-licence-en_0.pdfEUPL">EUPL</a> </td> </tr> 
</table>

**Description**

The package `pygridmap` supports the following methods/algorithms: 
* basic geometric/set operations and geoprocessing over regular grid with the [`GridProcessor`](pygridmap/base.py) class,
* user-defined [rasterisation](https://en.wikipedia.org/wiki/Rasterisation) of polygon regions (_e.g._, vector boundaries) into regular grids of various resolutions (_i.e._, dimension of square grid cells) with the [`GridMaker`](pygridmap/gridding.py) class (see also [GridMaker](https://github.com/eurostat/GridMaker)),
* customised [vector overlay operations](https://docs.qgis.org/3.10/en/docs/user_manual/processing_algs/qgis/vectoroverlay.html) (_e.g._, intersection, union, overlay/overlap and merging) between any polygon layer and a regular grid with the [`Gridoverlay`](pygridmap/overlay.py) class.

**Quick install and start**

**Usage**

*In your `Python` script*

*Through the `bash` command*

**Note**

The operations above are making an extensive use of the geometric/geospatial `Python` module [`geopandas`](https://geopandas.org/) since it supports essential features (most of them derived from the [`shapely`](https://shapely.readthedocs.io/en/latest/manual.html) [`fiona`](https://fiona.readthedocs.io/en/latest/manual.html) libraries), like:
  * basic [set](https://geopandas.org/set_operations.html) and [geometric](https://geopandas.org/geometric_manipulations.html) operations, and
  * [spatial indexing](https://geopandas.org/mergingdata.html?highlight=spatial%20index) through the [R-Tree algorithm](https://automating-gis-processes.github.io/site/notebooks/L3/spatial_index.html).
  
These implementations are probably not optimal since they are wrappers to the [`GEOS`](https://trac.osgeo.org/geos/) library, itself a port of the Java Topology Suite ([`JTS`](https://projects.eclipse.org/projects/locationtech.jts)), similarly to what is used in QGis.

They also adopt (customised, using the [`multiprocessing`](https://docs.python.org/3/library/multiprocessing.html) module) multitprocessor tiling processing together with (native) vector processing whenever possible (the regular grid often makes the computation embarassingly parallel) in order to take advantage of multiprocessor compute capabilities.

**<a name="Software"></a>Software resources/dependencies**

* Data handling: [`numpy`](https://numpy.org/) and [`pandas`](http://pandas.pydata.org).
* Geospatial data processing: [`geopandas`](http://geopandas.org),[`fiona`](https://fiona.readthedocs.io/en/latest/manual.html) libraries) and [shapely](https://pypi.org/project/Shapely/).
* Map visualisations: [`ipyleaflet`](https://github.com/jupyter-widgets/ipyleaflet) or [`folium`](https://github.com/python-visualization/folium).

**<a name="References"></a>References**
