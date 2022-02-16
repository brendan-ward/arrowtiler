# arrowtiler

`arrowtiler` creates Mapbox Vector Tiles from Feather (Apache Arrow IPC) format
files created according to the
[geo-arrow-spec](https://github.com/geopandas/geo-arrow-spec).

Feather files are a compact binary format that supports very fast I/O, and
GeoArrow files are Feather files that contain geometry objects in addition to
other data types.

This package aims to leverage the fast I/O and controlled data types of GeoArrow
combined with better control over simplification and rendering of geometries into
tiles.

WARNING: this package is under heavy development! APIs, functionality, etc are
all subject to change. We are still validating output tiles against those
created using tippecanoe.

Requires:

-   Go >= 1.17
-   GEOS >= 3.11.0

## Why arrowtiler?

The primary tool for creating vector tiles from file sources is
[tippecanoe](https://github.com/mapbox/tippecanoe).

Tippecanoe is great; we've used it in a bunch of projects.

However, we've run into 2 issues with tippecanoe:

-   serializing / deserializing GeoJSON is super slow, involves large
    intermediate files, and does not do a great job of preserving data types
-   it has several configuration parameters that control things like point
    density and simplification, though it is not clear just which parameters can
    be combined, and these have not worked well in our use cases

We have started using GeoArrow files as intermediates in our data processing
pipelines and really appreciate the fast I/O and preservation of data types.
While it is possible to build in support for geospatial formats using
[GDAL](https://gdal.org/), we've found that Feather files are quite a bit faster
for I/O.

At present, geometries are encoded as WKB bytes in a GeoArrow file.

Input geometries must be in geographic coordinates (WGS 1984 / "EPSG:4326").

Geometries must be of a uniform type: all points / multipoints, linestrings /
multilinestrings, or polygons / multipolygons. GeometryCollections are not
supported.

## How it works

This package uses the [GEOS](http://libgeos.org/) C API to provide the
underlying geometry objects and operations. Unlike other implementations of
GEOS in Go, this uses a vectorized approach to limit the number of calls to CGO.
Instead, GEOS geometries are handled as slices of pointers in Go, and the entire
slice is passed to CGO to perform operations in C using GEOS, which returns
either a new array of geometries, or other array-type return values (e.g.,
integer indexes from STRtree query operation).

Input geometries must be in geographic (latitude / longitude) coordinates.
These are clipped to world bounds supported by spherical Mercator (-180, -85,
180, 85) and projected to spherical Mercator.

We then create an STRtree (spatial index) of those geometries for very fast
querying by the bounding box of each tile.

All non-geometry attributes are encoded into the MVT value format up front
rather than each time an individual tile is rendered.

For each tile, geometries are projected to tile coordinates and converted into
MVT geometry format. At present, only minimal simplification is performed, and
no limit of tile number of features or size is enforced.

## Credits

Inspired by:

-   [tippecanoe](https://github.com/mapbox/tippecanoe)
-   [Golang MVT writer](https://github.com/tidwall/mvt)
-   Some Geographic / Spherical Mercator calculations derived from [mercantile](https://github.com/mapbox/mercantile)
