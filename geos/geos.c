#include "float.h"
#include "math.h"
#include "string.h"
#include <geos_c.h>
#include <stdio.h>
#include <stdlib.h>

#include "geos.h"
#include "kvec.h"

// MVT encoding parameters
const uint32_t MOVETO = 1;
const uint32_t LINETO = 2;
const uint32_t CLOSEPATH = 7;

const unsigned char MVT_GEOM_POINT = 1;
const unsigned char MVT_GEOM_LINE = 2;
const unsigned char MVT_GEOM_POLYGON = 3;

// projection constants
const double MERCATOR_ORIGIN = 20037508.342789;
const double CE = 2 * 6378137.0 * M_PI;
const double RE = 6378137.0;
const double DEG2RAD = M_PI / 180;
// calculate value to multiply by longitude to get Mercator
const double LON_FAC = RE * DEG2RAD;

// outermost bounds that can be safely represented in spherical Mercator
const double WORLD_BOUNDS[4] = {-180.0, -85.051129, 180.0, 85.051129};

// callback from GEOS to Go
void geos_error_handler(const char *message, void *userdata) {
  // FIXME: temp logging
  printf("GEOS ERROR: %s\n", message);
  void geos_errorMessageHandlerCallback(const char *, void *);
  geos_errorMessageHandlerCallback(message, userdata);
}

/* Frees each nested subarray in values that is non-NULL, then frees the
 * containing array of values and sizes.
 */
void free_uint32_subarrays(uint32_t **values, size_t *sizes, size_t size) {
  if (values == NULL || size == 0) {
    return;
  }

  for (int i = 0; i < size; i++) {
    if (values[i] != NULL) {
      free(values[i]);
    }
  }

  free(values);
  free(sizes);
}

/* Free all GEOS geometries in array.  Caller is responsible for freeing
 * containing array; it may be a Go slice.
 *
 */
void free_geometry_array(GEOSGeometry **geoms, size_t count) {
  if (geoms == NULL) {
    return;
  }

  GEOS_INIT;
  for (int i = 0; i < count; i++) {
    if (geoms[i] != NULL) {
      GEOSGeom_destroy_r(ctx, geoms[i]);
    }
  }
  GEOS_FINISH;
}

/* Read GEOS geometries from array of WKT strings.
 * Polygon and multi geometries will be normalized.
 * Any geometries that fail to parse properly are returned as NULL.
 *
 * Returns NULL on error.
 *
 * Parameters
 * ----------
 * wkts: pointer to array of WKT strings
 * count: number of WKT geometries
 *
 * Returns
 * -------
 * pointer to array of GEOS Geometry objects (same length as count)
 */
GEOSGeometry **from_wkt(char **wkts, size_t count) {
  GEOSGeometry **geoms = NULL;
  GEOSWKTReader *reader;
  int i, j;

  GEOS_INIT;

  reader = GEOSWKTReader_create_r(ctx);
  geoms = malloc(sizeof(GEOSGeometry *) * count);

  for (i = 0; i < count; i++) {
    geoms[i] = NULL; // default to NULL

    if (wkts[i] == NULL || strlen(wkts[i]) == 0) {
      continue;
    }

    geoms[i] = GEOSWKTReader_read_r(ctx, reader, wkts[i]);
    if (geoms[i] == NULL) {
      // since this is a parsing error, we don't want to fail the entire
      // operation; output for this geom will be NULL
      printf("WARNING: could not read WKT geometry from at index %i\n", i);
    }

    // normalize coordinate order and subgeometry order for multi types
    if (GEOSGeomTypeId_r(ctx, geoms[i]) != GEOS_POINT) {
      if (GEOSNormalize_r(ctx, geoms[i]) == -1) {
        printf("ERROR: could not normalize polygon\n");
        errstate = GEOSERROR;
        goto finish;
      }
    }
  }

// cleanup on error
finish:
  if (errstate == GEOSERROR) {
    if (geoms != NULL) {
      for (j = 0; j < i; j++) {
        if (geoms[j] != NULL) {
          GEOSGeom_destroy_r(ctx, geoms[j]);
        }
      }
      free(geoms);
      geoms = NULL;
    }
  }

  GEOS_FINISH;

  return geoms;
}

/**
 * Write the GEOS geometries to WKT.
 * Any NULL geometry is returned as NULL.
 * Caller is responsible for freeing the strings and containing array.
 *
 * Parameters
 * ----------
 * geometries: array of GEOS geometry objects
 * count: number of geometries in array
 * precision: decimal precision for formatting; -1 indicates no precision set
 *
 * Returns
 * -------
 * pointer to array of WKT strings.
 * Returns NULL on error.
 */
char **to_wkt(GEOSGeometry **geometries, size_t count, int precision) {
  GEOSWKTWriter *writer;
  char **wkts = NULL;
  int i, j;

  GEOS_INIT;
  writer = GEOSWKTWriter_create_r(ctx);
  if (precision != -1) {
    GEOSWKTWriter_setRoundingPrecision_r(ctx, writer, precision);
  }

  wkts = malloc(sizeof(char *) * count);

  for (i = 0; i < count; i++) {
    wkts[i] = NULL; // default to NULL
    if (geometries[i] == NULL) {
      continue;
    }
    wkts[i] = GEOSWKTWriter_write_r(ctx, writer, geometries[i]);
    if (wkts[i] == NULL) {
      printf("ERROR: could not write geometry to WKT at index %i\n", i);
      errstate = GEOSERROR;
      goto finish;
    }
  }

// cleanup on error
finish:
  if (errstate == GEOSERROR) {
    if (wkts != NULL) {
      for (j = 0; j < i; j++) {
        if (wkts[j] != NULL) {
          free(wkts[i]);
        }
      }
      free(wkts);
      wkts = NULL;
    }
  }

  GEOS_FINISH;

  return wkts;
}

/* Read GEOS geometries from array of WKB bytes.
 * Polygon and multi geometries will be normalized.
 * Any geometries that fail to parse properly are returned as NULL.
 *
 * Returns NULL on error.
 *
 * Parameters
 * ----------
 * wkbs: pointer to array of WKB byte arrays
 * count: number of WKB geometries
 *
 * Returns
 * -------
 * pointer to array of GEOS Geometry objects (same length as count)
 */
GEOSGeometry **from_wkb(unsigned char **wkbs, size_t *counts, size_t count) {
  GEOSGeometry **geoms = NULL;
  GEOSWKBReader *reader;

  int i, j;

  GEOS_INIT;
  reader = GEOSWKBReader_create_r(ctx);
  geoms = malloc(sizeof(GEOSGeometry *) * count);

  for (i = 0; i < count; i++) {
    geoms[i] = NULL; // default to NULL

    if (wkbs[i] == NULL) {
      continue;
    }

    geoms[i] = GEOSWKBReader_read_r(ctx, reader, (const unsigned char *)wkbs[i],
                                    counts[i]);
    if (geoms[i] == NULL) {
      // since this is a parsing error, we don't want to fail the entire
      // operation; output for this geom will be NULL
      printf("WARNING: could not read WKB geometry from at index %i\n", i);
    }

    // normalize coordinate order and subgeometry order for multi types
    if (GEOSGeomTypeId_r(ctx, geoms[i]) != GEOS_POINT) {
      if (GEOSNormalize_r(ctx, geoms[i]) == -1) {
        printf("ERROR: could not normalize polygon\n");
        errstate = GEOSERROR;
        goto finish;
      }
    }
  }

// cleanup on error
finish:
  if (errstate == GEOSERROR) {
    if (geoms != NULL) {
      for (j = 0; j < i; j++) {
        if (geoms[j] != NULL) {
          GEOSGeom_destroy_r(ctx, geoms[j]);
        }
      }
      free(geoms);
      geoms = NULL;
    }
  }

  GEOS_FINISH;

  return geoms;
}

/* Create a Polygon from bounding coordinates.
 *
 * Parameters
 * ----------
 * ctx: GEOS context ctx
 * xmin: minimum X value
 * ymin: minimum Y value
 * xmax: maximum X value
 * ymax: maximum Y value
 *
 * Returns
 * -------
 * GEOSGeometry* on success (owned by caller) or
 * NULL on failure
 */
GEOSGeometry *create_box(GEOSContextHandle_t ctx, double xmin, double ymin,
                         double xmax, double ymax) {

  GEOSCoordSequence *coords = NULL;
  GEOSGeometry *geom = NULL;
  GEOSGeometry *ring = NULL;

  // Construct coordinate sequence and set vertices
  coords = GEOSCoordSeq_create_r(ctx, 5, 2);
  if (coords == NULL) {
    return NULL;
  }

  // TODO: setxy
  // CCW from bottom right (xmax, ymin) to match shapely
  if (!(GEOSCoordSeq_setX_r(ctx, coords, 0, xmax) &&
        GEOSCoordSeq_setY_r(ctx, coords, 0, ymin) &&
        GEOSCoordSeq_setX_r(ctx, coords, 1, xmax) &&
        GEOSCoordSeq_setY_r(ctx, coords, 1, ymax) &&
        GEOSCoordSeq_setX_r(ctx, coords, 2, xmin) &&
        GEOSCoordSeq_setY_r(ctx, coords, 2, ymax) &&
        GEOSCoordSeq_setX_r(ctx, coords, 3, xmin) &&
        GEOSCoordSeq_setY_r(ctx, coords, 3, ymin) &&
        GEOSCoordSeq_setX_r(ctx, coords, 4, xmax) &&
        GEOSCoordSeq_setY_r(ctx, coords, 4, ymin))) {
    if (coords != NULL) {
      GEOSCoordSeq_destroy_r(ctx, coords);
    }

    return NULL;
  }

  // Construct linear ring then use to construct polygon
  // Note: coords are owned by ring; if ring fails to construct, it will
  // automatically clean up the coords
  ring = GEOSGeom_createLinearRing_r(ctx, coords);
  if (ring == NULL) {
    return NULL;
  }

  // Note: ring is owned by polygon; if polygon fails to construct, it will
  // automatically clean up the ring
  geom = GEOSGeom_createPolygon_r(ctx, ring, NULL, 0);
  if (geom == NULL) {
    return NULL;
  }

  return geom;
}

/* Get the bounds of a geometry
 *
 * Parameters
 * ----------
 * ctx: GEOS context ctx
 * geom: GEOS geometry object
 * xmin: double pointer to receive xmin value
 * ymin: double pointer to receive ymin value
 * xmax: double pointer to receive xmax value
 * ymax: double pointer to receive ymax value
 *
 * * Returns 0 on failure
 */
int get_geometry_bounds(GEOSContextHandle_t ctx, GEOSGeometry *geom,
                        double *xmin, double *ymin, double *xmax,
                        double *ymax) {

  // can't get bounds of empty geometries; fail
  if (geom == NULL || GEOSisEmpty_r(ctx, geom)) {
    return 0;
  }

  if (!(GEOSGeom_getXMin_r(ctx, geom, xmin) &&
        GEOSGeom_getYMin_r(ctx, geom, ymin) &&
        GEOSGeom_getXMax_r(ctx, geom, xmax) &&
        GEOSGeom_getYMax_r(ctx, geom, ymax))) {
    return 0;
  }

  return 1;
}

/* Get the total bounds of an array of geometries
 *
 * Parameters
 * ----------
 * ctx: GEOS context ctx
 * geometries: pointer to array fo GEOS geometry objects
 * count: number of geometries in array
 * xmin: double pointer to receive xmin value
 * ymin: double pointer to receive ymin value
 * xmax: double pointer to receive xmax value
 * ymax: double pointer to receive ymax value
 *
 *
 * Returns 0 on failure
 */
int get_total_bounds(GEOSGeometry **geometries, size_t count, double *xmin,
                     double *ymin, double *xmax, double *ymax) {

  double total_xmin = DBL_MAX;
  double total_ymin = DBL_MAX;
  double total_xmax = -DBL_MAX;
  double total_ymax = -DBL_MAX;
  double geom_xmin, geom_ymin, geom_xmax, geom_ymax;
  int i;
  GEOSGeometry *geom;

  GEOS_INIT;

  for (i = 0; i < count; i++) {
    geom = geometries[i];

    // skip empty / null geometries
    if (geom == NULL || GEOSisEmpty_r(ctx, geom)) {
      continue;
    }

    if (!get_geometry_bounds(ctx, geom, &geom_xmin, &geom_ymin, &geom_xmax,
                             &geom_ymax)) {
      return 0;
    }
    if (geom_xmin < total_xmin) {
      total_xmin = geom_xmin;
    }
    if (geom_ymin < total_ymin) {
      total_ymin = geom_ymin;
    }
    if (geom_xmax > total_xmax) {
      total_xmax = geom_xmax;
    }
    if (geom_ymax > total_ymax) {
      total_ymax = geom_ymax;
    }
  }
  *xmin = total_xmin;
  *ymin = total_ymin;
  *xmax = total_xmax;
  *ymax = total_ymax;

  GEOS_FINISH;

  return 1;
}

/* Callback called by strtree_query with item data of each intersecting geometry
 * and a dynamic vector to push that item onto.
 *
 * Parameters
 * ----------
 * item: index of intersected geometry in the tree
 * user_data: pointer to dynamic vector
 */
void query_callback(void *item, void *user_data) {
  kv_push(uint32_t, *(uint32_vec_t *)user_data, *(uint32_t *)item);
}

/* dummy query callback used during tree construction */
void dummy_query_callback(void *item, void *user_data) {}

/* Create an STRtree of geometries
 *
 * Parameters
 * ----------
 * geometries: array of GEOSGeometry pointers
 * count: size of the array of geometries
 *
 * Returns
 * -------
 * STRtree wrapper around GEOSStrtree
 */
STRtree *create_tree(GEOSGeometry **geometries, size_t count) {
  GEOSGeometry *geom;
  GEOSGeometry *pt;
  STRtree *tree;
  GEOSSTRtree *geos_tree;
  uint32_t *indexes;
  uint32_t i;
  size_t tree_size = 0;
  size_t node_capacity = 10;

  GEOS_INIT;

  geos_tree = GEOSSTRtree_create_r(ctx, node_capacity);
  if (geos_tree == NULL) {
    printf("ERROR: could not construct STRtree");
    errstate = GEOSERROR;
    goto finish;
  }

  indexes = (uint32_t *)malloc(count * sizeof(uint32_t));

  for (i = 0; i < count; i++) {
    indexes[i] = i;
    geom = geometries[i];
    if (!(geom == NULL || GEOSisEmpty_r(ctx, geom))) {
      GEOSSTRtree_insert_r(ctx, geos_tree, geom, &(indexes[i]));
      tree_size++;
    }
  }

  // run a dummy query to force tree to construct
  pt = GEOSGeom_createPointFromXY_r(ctx, 0, 0);
  if (pt == NULL) {
    printf("ERROR: could not init STRtree");
    errstate = GEOSERROR;
    goto finish;
  }
  GEOSSTRtree_query_r(ctx, geos_tree, pt, dummy_query_callback, NULL);

  tree = (STRtree *)malloc(sizeof(STRtree *));
  tree->geos_tree = geos_tree;
  tree->indexes = indexes;

finish:
  if (errstate == GEOSERROR) {
    if (indexes != NULL) {
      free(indexes);
    }
    if (geos_tree != NULL) {
      GEOSSTRtree_destroy_r(ctx, geos_tree);
    }
    if (tree != NULL) {
      free(tree);
    }
  }

  GEOS_FINISH;

  return tree;
}

void destroy_tree(STRtree *tree) {
  if (tree == NULL) {
    return;
  }

  GEOS_INIT;
  if (tree != NULL) {
    GEOSSTRtree_destroy_r(ctx, tree->geos_tree);
    free(tree->indexes);
    free(tree);
  }
  GEOS_FINISH;
}

/* Query the tree for geometries within xmin, ymin, xmax, ymax and return
 * an array of integer indexes into tree geometries of the hits and a hit_count.
 *
 * Returns
 * -------
 * 0 on failure, 1 on success
 *
 * Parameters
 * ----------
 * tree: STRtree pointer (wrapper around GEOSStrtree)
 * xmin...ymax: bounding coordinates to query tree
 * hits: pointer to array of uint32s that will receive the integer indexes
 *  of intersected geometries in the tree.  NOTE: the caller must free this.
 * hit_count: number of hits returned in array
 */
int query_tree(STRtree *tree, double xmin, double ymin, double xmax,
               double ymax, uint32_t **hits, size_t *hit_count) {

  GEOS_INIT;

  GEOSGeometry *envelope = create_box(ctx, xmin, ymin, xmax, ymax);
  if (envelope == NULL) {
    printf("could not create box to query STRtree");
    errstate = GEOSERROR;
    goto finish;
  }

  uint32_vec_t indexes;

  kv_init(indexes);
  GEOSSTRtree_query_r(ctx, tree->geos_tree, envelope, query_callback, &indexes);

  // Caller must free hits!
  *hits = indexes.a;
  *hit_count = kv_size(indexes);

finish:
  if (errstate == GEOSERROR) {
    kv_destroy(indexes);
  }

  if (envelope != NULL) {
    GEOSGeom_destroy_r(ctx, envelope);
  }

  GEOS_FINISH;

  return 1;
}

/**
 * Project geographic coordinates to spherical Mercator.
 * WARNING: assumes that geometry has already been clipped to
 * Mercator world bounds.
 *
 * Parameters
 * ----------
 * x: input x
 * y: input y
 * x_out: pointer to receive projected x
 * y_out: pointer to receive projected y
 * userdata: not used for this function
 */
int geo_to_mercator_coords(double x, double y, double *x_out, double *y_out,
                           void *userdata) {
  *x_out = LON_FAC * x;
  *y_out = RE * log(tan((M_PI * 0.25) + (0.5 * DEG2RAD * y)));

  // if y is very small, make it 0
  if (*y_out > -1e-6 && *y_out < 1e-6) {
    *y_out = 0;
  }

  return 1;
}

/**
 * Project Mercator coordinates to tile coordinates.
 *
 * Parameters
 * ----------
 * x: input x
 * y: input y
 * x_out: pointer to receive projected x
 * y_out: pointer to receive projected y
 * userdata: struct containing xmin, ymin, xmax, ymax for scaling
 */
int mercator_to_tile_coords(double x, double y, double *x_out, double *y_out,
                            void *userdata) {

  ScaleParams *params = (ScaleParams *)(userdata);
  *x_out = x * params->xscale + params->xoffset;
  *y_out = y * params->yscale + params->yoffset;

  return 1;
}

/*
 * Applies coord_func to the coordinates in geom for simple geometry types.
 * Geometry must be non-empty.
 *
 * Parameters
 * ----------
 * ctx: GEOS context handle
 * geom: GEOS Geometry object (Point, LineString, LinearRing only)
 * coord_func: function applied to coordinates to perform projection (see
 * coord_project_func type) Returns
 * --------
 * pointer to GEOS Geometry object
 */
GEOSGeometry *project_geometry(GEOSContextHandle_t ctx,
                               const GEOSGeometry *geom,
                               coord_project_func *coord_func, void *userdata) {
  GEOSGeometry *out = NULL;
  const GEOSCoordSequence *seq;
  GEOSCoordSequence *out_seq;
  unsigned int size;
  int geom_type;
  double lon, lat;
  double x, y;

  seq = GEOSGeom_getCoordSeq_r(ctx, geom);
  if (seq == NULL) {
    printf("ERROR: could not get coordinate sequence for geometry\n");
    return NULL;
  }

  if (!GEOSCoordSeq_getSize_r(ctx, seq, &size)) {
    printf("ERROR: could not get coordinates for coordinate sequence\n");
    return NULL;
  }

  out_seq = GEOSCoordSeq_create_r(ctx, size, 2);
  if (out_seq == NULL) {
    printf("ERROR: could not create new coordinate sequence\n");
    return NULL;
  }

  for (int i = 0; i < size; i++) {
    if (!GEOSCoordSeq_getXY_r(ctx, seq, i, &lon, &lat)) {
      printf("ERROR: could not get coordinates for coordinate sequence\n");
      goto fail;
    }

    coord_func(lon, lat, &x, &y, userdata);

    if (!GEOSCoordSeq_setXY_r(ctx, out_seq, i, x, y)) {
      printf("ERROR: could not set coordinates into coordinate sequence\n");
      goto fail;
    }
  }

  // Note: the coordinate sequence becomes owned by the new geometry
  switch (GEOSGeomTypeId_r(ctx, geom)) {
  case GEOS_POINT: {
    out = GEOSGeom_createPoint_r(ctx, out_seq);
    break;
  }
  case GEOS_LINESTRING: {
    out = GEOSGeom_createLineString_r(ctx, out_seq);
    break;
  }
  case GEOS_LINEARRING: {
    out = GEOSGeom_createLinearRing_r(ctx, out_seq);
    break;
  }
  }

  if (out == NULL) {
    goto fail;
  }

  return out;

fail:
  GEOSCoordSeq_destroy_r(ctx, out_seq);
  return NULL;
}

/*
 * Applies coord_func to the coordinates in geom for singular polygon type
 * Geometry must be non-empty.
 *
 * Parameters
 * ----------
 * ctx: GEOS context handle
 * geom: GEOS Geometry object (Polygon only)
 * coord_func: function applied to coordinates to perform projection (see
 * coord_project_func type) Returns
 * --------
 * pointer to GEOS Geometry object
 */
GEOSGeometry *project_polygon(GEOSContextHandle_t ctx, const GEOSGeometry *geom,
                              coord_project_func *coord_func, void *userdata) {

  const GEOSGeometry *ring;
  GEOSGeometry *out, *out_ring;
  GEOSGeometry **out_rings;
  int num_rings, i, j;
  char errstate = SUCCESS;

  ring = GEOSGetExteriorRing_r(ctx, geom);
  if (ring == NULL) {
    printf("ERROR: could not get exterior ring for polygon\n");
    return NULL;
  }

  out_ring = project_geometry(ctx, ring, coord_func, userdata);
  if (out_ring == NULL) {
    // error should already be set
    return NULL;
  }

  num_rings = GEOSGetNumInteriorRings_r(ctx, geom);
  out_rings = malloc(sizeof(GEOSGeometry *) * num_rings);

  for (i = 0; i < num_rings; i++) {
    ring = GEOSGetInteriorRingN_r(ctx, geom, i);
    if (ring == NULL) {
      printf("ERROR: could not get interior ring %i for polygon\n", i);
      errstate = GEOSERROR;
      goto finish;
    }

    out_rings[i] = project_geometry(ctx, ring, coord_func, userdata);
    if (out_rings[i] == NULL) {
      // error should already be set
      errstate = GEOSERROR;
      goto finish;
    }
  }

  // new rings become owned by out geometry
  out = GEOSGeom_createPolygon_r(ctx, out_ring, out_rings, num_rings);
  if (out == NULL) {
    printf("ERROR: could not create new polygon\n");
    errstate = GEOSERROR;
  }

finish:
  if (errstate == GEOSERROR) {
    if (out_ring != NULL) {
      GEOSGeom_destroy_r(ctx, out_ring);
    }
    if (out_rings != NULL) {
      for (j = 0; j < i; j++) {
        if (out_rings[i] != NULL) {
          GEOSGeom_destroy_r(ctx, out_rings[i]);
        }
      }
    }
  }

  if (out_rings != NULL) {
    free(out_rings);
  }

  return out;
}

/*
 * Applies coord_func to the coordinates in geom for singular polygon type
 * Geometry must be non-empty.
 *
 * Parameters
 * ----------
 * ctx: GEOS context handle
 * geom: GEOS Geometry object (Polygon only)
 * coord_func: function applied to coordinates to perform projection (see
 * coord_project_func type) Returns
 * --------
 * pointer to GEOS Geometry object
 */
GEOSGeometry *project_multi_geometry(GEOSContextHandle_t ctx,
                                     const GEOSGeometry *geom,
                                     coord_project_func *coord_func,
                                     void *userdata) {

  const GEOSGeometry *subgeom;
  GEOSGeometry *out;
  GEOSGeometry **out_sub_geoms;
  int i, j, num_geoms, geom_type;
  char errstate = SUCCESS;
  geometry_project_func *proj_func = NULL;

  geom_type = GEOSGeomTypeId_r(ctx, geom);
  num_geoms = GEOSGetNumGeometries_r(ctx, geom);
  out_sub_geoms = malloc(sizeof(GEOSGeometry *) * num_geoms);

  if (geom_type == GEOS_MULTIPOLYGON) {
    proj_func = project_polygon;
  } else {
    proj_func = project_geometry;
  }

  for (i = 0; i < num_geoms; i++) {
    subgeom = GEOSGetGeometryN_r(ctx, geom, i);
    if (subgeom == NULL) {
      printf("ERROR: Could not get geometry %i\n", i);
      errstate = GEOSERROR;
      goto finish;
    }

    out_sub_geoms[i] = proj_func(ctx, subgeom, coord_func, userdata);

    if (out_sub_geoms[i] == NULL) {
      // error should already be reported
      errstate = GEOSERROR;
      goto finish;
    }
  }

  out = GEOSGeom_createCollection_r(ctx, geom_type, out_sub_geoms, num_geoms);
  if (out == NULL) {
    printf("ERROR: could not create collection from new geometries\n");
    errstate = GEOSERROR;
    goto finish;
  }

finish:
  if (errstate == GEOSERROR) {
    if (out_sub_geoms != NULL) {
      for (j = 0; j < i; j++) {
        if (out_sub_geoms[i] != NULL) {
          GEOSGeom_destroy_r(ctx, out_sub_geoms[i]);
        }
      }
    }
  }

  if (out_sub_geoms != NULL) {
    free(out_sub_geoms);
  }

  return out;
}

/**
 * Project geometries from geographic coordinates to spherical Mercator.
 * Geometries are clipped to Mercator world bounds if necessary.
 *
 * GeometryCollections are not supported.
 *
 * Parameters
 * ----------
 * geometries: array of GEOS Geometry objects
 * count: length of array
 *
 * Returns
 * -------
 * pointer to array of GEOS Geometry objects
 * Returns NULL on error.
 */
GEOSGeometry **project_to_mercator(GEOSGeometry **geometries, size_t count) {

  GEOSGeometry *geom = NULL;
  GEOSGeometry *new_geom = NULL;
  const GEOSGeometry *subgeom;
  const GEOSCoordSequence *seq;
  GEOSGeometry **out_geoms = NULL;
  int i, j, geom_type;
  double geom_xmin, geom_ymin, geom_xmax, geom_ymax;
  geometry_project_func *proj_func = NULL;

  GEOS_INIT;

  out_geoms = malloc(sizeof(GEOSGeometry *) * count);

  for (i = 0; i < count; i++) {
    out_geoms[i] = NULL; // default everything to NULL
    new_geom = NULL;

    geom = geometries[i];
    if (geom == NULL || GEOSisEmpty_r(ctx, geom)) {
      continue;
    }

    geom_type = GEOSGeomTypeId_r(ctx, geom);

    if (geom_type == GEOS_GEOMETRYCOLLECTION) {
      printf("ERROR: GeometryCollections are not supported\n");
      errstate = GEOSERROR;
      goto finish;
    }

    if (!get_geometry_bounds(ctx, geom, &geom_xmin, &geom_ymin, &geom_xmax,
                             &geom_ymax)) {
      printf("ERROR: could not calculate bounds of geometry at index %i\n", i);
      errstate = GEOSERROR;
      goto finish;
    }

    if (!(geom_xmin > WORLD_BOUNDS[0] && geom_ymin > WORLD_BOUNDS[1] &&
          geom_xmax < WORLD_BOUNDS[2] && geom_ymax < WORLD_BOUNDS[3])) {

      // printf("INFO: Clipping geometry to world bounds\n");

      // TODO: wrap X if needed

      out_geoms[i] =
          GEOSClipByRect_r(ctx, geom, WORLD_BOUNDS[0], WORLD_BOUNDS[1],
                           WORLD_BOUNDS[2], WORLD_BOUNDS[3]);

      if (out_geoms[i] == NULL) {
        printf("ERROR: could not clip geometry to tile bounds plus buffer\n");
        errstate = GEOSERROR;
        goto finish;
      }

      if (GEOSisEmpty_r(ctx, out_geoms[i])) {
        // outside world bounds, set it to NULL
        GEOSGeom_destroy_r(ctx, out_geoms[i]);
        out_geoms[i] = NULL;
        continue;
      }

      // get type again in case changed by clip
      geom_type = GEOSGeomTypeId_r(ctx, out_geoms[i]);
      geom = out_geoms[i];
    }

    if (geom_type == GEOS_POINT || geom_type == GEOS_LINESTRING) {
      proj_func = project_geometry;
    } else if (geom_type == GEOS_POLYGON) {
      proj_func = project_polygon;
    } else {
      proj_func = project_multi_geometry;
    }

    new_geom = proj_func(ctx, geom, geo_to_mercator_coords, NULL);
    if (new_geom == NULL) {
      // error should already be reported
      errstate = GEOSERROR;
      goto finish;
    }

    // delete prior, if exists, before setting new values
    if (out_geoms[i] != NULL) {
      GEOSGeom_destroy_r(ctx, out_geoms[i]);
    }
    out_geoms[i] = new_geom;
  }

finish:
  if (errstate == GEOSERROR) {
    if (out_geoms != NULL) {
      // free any geometries created here
      for (j = 0; j < i; j++) {
        if (out_geoms[j] != NULL) {
          GEOSGeom_destroy_r(ctx, out_geoms[j]);
        }
      }
      free(out_geoms);
      out_geoms = NULL;
    }
  }

  GEOS_FINISH;

  return out_geoms;
}

/**
 * Clip and project geometries to tile.
 * Geometries must already be in Mercator coordinates and clipped to world
 * bounds. Geometries beyond edge of tile will be clipped to a 256px buffer
 * around tile. Geometries are scaled to tile coordinates.
 *
 * Parameters
 * ----------
 * geometries: pointer to array of GEOS geometry objects
 * count: size of geometry array
 * xmin...ymax: Web Mercator bounds of tile
 */
GEOSGeometry **clip_project_to_tile(GEOSGeometry **geometries, size_t count,
                                    double xmin, double ymin, double xmax,
                                    double ymax, uint16_t extent) {

  double geom_xmin, geom_ymin, geom_xmax, geom_ymax;
  double x, y;
  int i, j, geom_type, new_geom_type, num_coords, num_sub_geoms;
  char is_same_type;
  GEOSGeometry *geom, *new_geom;
  geometry_project_func *proj_func = NULL;

  GEOS_INIT;

  GEOSGeometry **out_geoms = malloc(sizeof(GEOSGeometry *) * count);

  double xrange = xmax - xmin;
  double yrange = ymax - ymin;
  double min_width = (xrange / extent) / 2.0;
  double min_height = (yrange / extent) / 2.0;

  // buffer the tile by 256 units (same as PostGIS default) out of the extent
  double xbuffer = xrange * (256.0 / (double)extent);
  double ybuffer = yrange * (256.0 / (double)extent);

  // values for scaling coordinates
  double xscale = extent / xrange;
  double yscale = -(extent / yrange);
  double xoffset = -xmin * xscale;
  double yoffset = -ymax * yscale;

  const ScaleParams scale_params = {xscale, xoffset, yscale, yoffset};

  for (i = 0; i < count; i++) {
    geom = geometries[i];
    out_geoms[i] = NULL; // default everything to NULL

    if (geom == NULL || GEOSisEmpty_r(ctx, geom)) {
      continue;
    }

    geom_type = GEOSGeomTypeId_r(ctx, geom);

    // check bounds: if too small to represent, then return null (similar to
    // what PostGIS does)
    if (geom_type == GEOS_POLYGON || geom_type == GEOS_MULTIPOLYGON ||
        geom_type == GEOS_LINESTRING || geom_type == GEOS_MULTILINESTRING) {
      if (!get_geometry_bounds(ctx, geom, &geom_xmin, &geom_ymin, &geom_xmax,
                               &geom_ymax)) {
        printf("ERROR: could not calculate bounds of geometry\n");
        errstate = GEOSERROR;
        goto finish;
      }

      if ((geom_xmax - geom_xmin) < min_width &&
          (geom_ymax - geom_ymin) < min_height) {
        // too small: skip
        // printf("INFO: geometry is too small to represent in tile\n");
        continue;
      }
    }

    // if geometry is a point, we already ensured it is contained in tile
    // using the tree query; otherwise if its bounds are entirely within the
    // tile, it also doesn't need to be clipped
    if (!(geom_type == GEOS_POINT || (geom_xmin > xmin && geom_xmax < xmax &&
                                      geom_ymin > ymin && geom_ymax < ymax))) {
      // printf("INFO: Clipping geometry by tile bounds plus buffer\n");

      out_geoms[i] = GEOSClipByRect_r(ctx, geom, xmin - xbuffer, ymin - ybuffer,
                                      xmax + xbuffer, ymax + ybuffer);
      if (out_geoms[i] == NULL) {
        printf("ERROR: could not clip geometry to tile bounds plus buffer\n");
        errstate = GEOSERROR;
        goto finish;
      }

      // in case type changed by clip
      geom_type = GEOSGeomTypeId_r(ctx, out_geoms[i]);
      geom = out_geoms[i];
    }

    if (geom_type == GEOS_POINT || geom_type == GEOS_LINESTRING) {
      proj_func = project_geometry;
    } else if (geom_type == GEOS_POLYGON) {
      proj_func = project_polygon;
    } else {
      proj_func = project_multi_geometry;
    }

    new_geom =
        proj_func(ctx, geom, mercator_to_tile_coords, (void *)(&scale_params));
    if (new_geom == NULL) {
      // error should already be reported
      errstate = GEOSERROR;
      goto finish;
    }

    // delete prior, if exists, before setting new values
    if (out_geoms[i] != NULL) {
      GEOSGeom_destroy_r(ctx, out_geoms[i]);
    }
    out_geoms[i] = new_geom;

    // Reduce precision to 1 pixel
    // use gridSize = 1 (1 pixel)
    // use default flags (preserve topology, drop collapsed geometries)
    geom = GEOSGeom_setPrecision_r(ctx, out_geoms[i], 1, 0);
    if (geom == NULL) {
      printf("ERROR: could not reduce precision for geometry; likely due to "
             "validity error\n");
      errstate = GEOSERROR;
      goto finish;
    }

    // destroy prev geom before overwriting
    GEOSGeom_destroy_r(ctx, out_geoms[i]);
    out_geoms[i] = geom;

    if (!(geom_type == GEOS_POINT || geom_type == GEOS_MULTIPOINT)) {
      // reduce unnecessary vertices on straight lines
      // TODO: verify this is giving us correct level of simplification here
      geom = GEOSSimplify_r(ctx, out_geoms[i], 1);
      if (geom == NULL) {
        printf("ERROR: could not simplify geometry\n");
        return NULL;
      }

      // destroy prev geom before overwriting
      GEOSGeom_destroy_r(ctx, out_geoms[i]);
      out_geoms[i] = geom;
    }

    // attempt to collapse multi type to singular
    // WARNING: this may be a big performance hit for polygons
    // TODO: consider geom_type == GEOS_MULTILINESTRING || geom_type ==
    // GEOS_MULTIPOLYGON
    if (geom_type == GEOS_MULTIPOINT) {
      geom = GEOSUnaryUnion_r(ctx, out_geoms[i]);
      if (geom == NULL) {
        printf("ERROR: could not simplify multi geometry to singular\n");
        return NULL;
      }

      // destroy prev geom before overwriting
      GEOSGeom_destroy_r(ctx, out_geoms[i]);
      out_geoms[i] = geom;
    }

    // geometry was collapsed or converted to different type
    // conversions from multi to singular are OK
    new_geom_type = GEOSGeomTypeId_r(ctx, out_geoms[i]);
    is_same_type =
        (geom_type == new_geom_type) ||
        (geom_type == GEOS_MULTIPOINT && new_geom_type == GEOS_POINT) ||
        (geom_type == GEOS_POINT && new_geom_type == GEOS_MULTIPOINT) ||
        (geom_type == GEOS_MULTILINESTRING &&
         new_geom_type == GEOS_LINESTRING) ||
        (geom_type == GEOS_LINESTRING &&
         new_geom_type == GEOS_MULTILINESTRING) ||
        (geom_type == GEOS_MULTIPOLYGON && new_geom_type == GEOS_POLYGON) ||
        (geom_type == GEOS_POLYGON && new_geom_type == GEOS_MULTIPOLYGON);
    if (GEOSisEmpty_r(ctx, geom) || !is_same_type) {
      // printf(
      //     "INFO: Geometry (type %i) was collapsed or converted to a different
      //     " "dimensionality (type %i, is empty? %i); it will not be included
      //     in " "tile\n", geom_type, new_geom_type, GEOSisEmpty_r(ctx, geom));
      out_geoms[i] = NULL;
      continue;
    }

    // check for and make valid (in place)
    if (geom_type != GEOS_POINT) {
      if (!GEOSisValid_r(ctx, geom)) {
        if (!GEOSMakeValid_r(ctx, geom)) {
          printf("ERROR: could not make geometry valid\n");
          return NULL;
        }
      }
    }
  }

finish:
  if (errstate == GEOSERROR) {
    if (out_geoms != NULL) {
      for (j = 0; j < i; j++) {
        if (out_geoms[j] != NULL) {
          GEOSGeom_destroy_r(ctx, out_geoms[j]);
        }
      }
      free(out_geoms);
      out_geoms = NULL;
    }
  }

  GEOS_FINISH;

  return out_geoms;
}

int encode_geometry(GEOSContextHandle_t ctx, const GEOSGeometry *geom,
                    uint32_t **values, int *pos, double *prevX, double *prevY) {
  unsigned int size;
  double x, y;
  uint32_t deltaX, deltaY;

  const GEOSCoordSequence *seq = GEOSGeom_getCoordSeq_r(ctx, geom);
  if (seq == NULL) {
    printf("ERROR: could not get coordinate sequence for geometry\n");
    return 0;
  }

  if (!GEOSCoordSeq_getSize_r(ctx, seq, &size)) {
    printf("ERROR: could not get coordinates for coordinate sequence\n");
    return 0;
  }

  for (int i = 0; i < size; i++) {
    if (!GEOSCoordSeq_getXY_r(ctx, seq, i, &x, &y)) {
      printf("ERROR: could not get coordinates for coordinate sequence\n");
      return 0;
    }

    deltaX = (uint32_t)(x - (*prevX));
    deltaY = (uint32_t)(y - (*prevY));

    (*values)[(*pos)] = (uint32_t)((deltaX << 1) ^ -(deltaX >> 31));
    (*values)[(*pos) + 1] = (uint32_t)((deltaY << 1) ^ -(deltaY >> 31));
    (*pos) += 2;
    (*prevX) = x;
    (*prevY) = y;

    if (i == 0 && size > 1) {
      // command/count for continuing the line
      (*values)[(*pos)] = (LINETO & 0x7) | ((size - 1) << 3);
      (*pos)++;
    }
  }

  return 1;
}

int encode_polygon(GEOSContextHandle_t ctx, const GEOSGeometry *geom,
                   uint32_t **values, int *pos, double *prevX, double *prevY) {

  const GEOSGeometry *ring = GEOSGetExteriorRing_r(ctx, geom);
  if (ring == NULL) {
    printf("ERROR: could not get exterior ring for polygon\n");
    return 0;
  }

  // encode outer ring
  (*values)[(*pos)] = (MOVETO & 0x7) | (1 << 3);
  (*pos)++;
  if (!encode_geometry(ctx, ring, values, pos, prevX, prevY)) {
    return 0;
  }
  // close ring
  (*values)[(*pos)] = CLOSEPATH;
  (*pos)++;

  int num_rings = GEOSGetNumInteriorRings_r(ctx, geom);
  for (int i = 0; i < num_rings; i++) {
    (*values)[(*pos)] = (MOVETO & 0x7) | (1 << 3);
    (*pos)++;

    ring = GEOSGetInteriorRingN_r(ctx, geom, i);
    if (ring == NULL) {
      printf("ERROR: could not get interior ring %i for polygon\n", i);
      return 0;
    }
    if (!encode_geometry(ctx, ring, values, pos, prevX, prevY)) {
      return 0;
    }
    // close ring
    (*values)[(*pos)] = CLOSEPATH;
    (*pos)++;
  }

  return 1;
}

/**
 * Encode array of GEOS geometries into MVT format.
 *
 * Parameters
 * ----------
 * geometries: pointer to array of GEOS geometry objects
 * count: size of geometry array
 * types: pointer to array of unsigned chars (byte) to receive MVT (not GEOS)
 * geometry types (will be 0 for NULL geoms) values: pointer to array of uint32
 * arrays to receive encoded values sizes: pointer to array that will hold sizes
 * of each uint32 array
 */
int encode_geometries(GEOSGeometry **geometries, size_t count,
                      unsigned char **types, uint32_t ***values,
                      size_t **sizes) {

  int i, j, num_sub_geoms, num_rings, num_coords, pos;
  GEOSGeometry *geom;
  const GEOSGeometry *subgeom;
  int geom_type;
  size_t length;
  double prevX, prevY;

  GEOS_INIT;

  unsigned char *geom_types = malloc(sizeof(int *) * count);
  uint32_t **buffers = malloc(sizeof(uint32_t *) * count);
  size_t *lengths = malloc(sizeof(uint32_t *) * count);

  // each encoded geometry is geometry type plus command / parameter
  // combinations
  for (i = 0; i < count; i++) {
    buffers[i] = NULL; // default to NULL
    geom_types[i] = 0;
    lengths[i] = 0;

    geom = geometries[i];

    if (geom == NULL || GEOSisEmpty_r(ctx, geom)) {
      continue;
    }

    geom_type = GEOSGeomTypeId_r(ctx, geom);
    num_sub_geoms = GEOSGetNumGeometries_r(ctx, geom);
    num_coords = GEOSGetNumCoordinates_r(ctx, geom);

    prevX = 0;
    prevY = 0;

    // TODO: sort sub geoms on Hilbert or Z order curve to minimize deltas
    // NOTE: oriented from upper left

    switch (geom_type) {
    case GEOS_POINT:
    case GEOS_MULTIPOINT: {
      geom_types[i] = MVT_GEOM_POINT;

      // command/count, <x, y>...
      length = num_sub_geoms + num_coords * 2;
      buffers[i] = malloc(sizeof(uint32_t) * length);
      lengths[i] = length;

      // command is repeated num_sub_geoms types
      buffers[i][0] = (MOVETO & 0x7) | (num_sub_geoms << 3);
      pos = 1;

      for (j = 0; j < num_sub_geoms; j++) {
        subgeom = GEOSGetGeometryN_r(ctx, geom, j);
        if (subgeom == NULL) {
          printf("Could not get geometry %i for line / multiline\n", j);
          return 0;
        }

        if (!encode_geometry(ctx, subgeom, &buffers[i], &pos, &prevX, &prevY)) {
          printf("Could not get encode coordinate sequence\n");
          return 0;
        }
      }

      break;
    }
    case GEOS_LINESTRING:
    case GEOS_MULTILINESTRING: {
      geom_types[i] = MVT_GEOM_LINE;

      // 2 commands per line
      length = (num_sub_geoms * 2) + (num_coords * 2);
      buffers[i] = malloc(sizeof(uint32_t) * length);
      lengths[i] = length;

      pos = 0;

      for (j = 0; j < num_sub_geoms; j++) {
        // MOVETO for start of each line
        buffers[i][pos] = (MOVETO & 0x7) | (1 << 3);
        pos++;

        subgeom = GEOSGetGeometryN_r(ctx, geom, j);
        if (subgeom == NULL) {
          printf("Could not get geometry %i for line / multiline\n", j);
          return 0;
        }

        if (!encode_geometry(ctx, subgeom, &buffers[i], &pos, &prevX, &prevY)) {
          printf("Could not get encode coordinate sequence\n");
          return 0;
        }
      }

      break;
    }
    case GEOS_POLYGON:
    case GEOS_MULTIPOLYGON: {
      geom_types[i] = MVT_GEOM_POLYGON;

      // count up rings first
      num_rings = num_sub_geoms;
      for (j = 0; j < num_sub_geoms; j++) {
        subgeom = GEOSGetGeometryN_r(ctx, geom, j);
        if (subgeom == NULL) {
          printf("Could not get geometry %i for polygon\n", j);
          return 0;
        }

        num_rings += GEOSGetNumInteriorRings_r(ctx, subgeom);
      }

      // 3 commands per ring: MoveTo(1), LineTo(n-1), ClosePath()
      length = (num_rings * 3) + (num_coords * 2);
      buffers[i] = malloc(sizeof(uint32_t) * length);
      lengths[i] = length;

      pos = 0;

      for (j = 0; j < num_sub_geoms; j++) {
        subgeom = GEOSGetGeometryN_r(ctx, geom, j);

        if (!encode_polygon(ctx, subgeom, &buffers[i], &pos, &prevX, &prevY)) {
          printf("Could not get encode coordinate polygon\n");
          return 0;
        }
      }

      break;
    }
    }
  }

  if (errstate == GEOSERROR) {
    if (buffers != NULL) {
      for (j = 0; j < i; j++) {
        if (buffers[j] != NULL) {
          free(buffers[j]);
        }
      }
      free(buffers);
    }

    buffers = NULL;
    if (lengths != NULL) {
      free(lengths);
      lengths = NULL;
    }

    if (geom_types != NULL) {
      free(geom_types);
      geom_types = NULL;
    }
  }

  *types = geom_types;
  *values = buffers;
  *sizes = lengths;

  GEOS_FINISH;

  return 1;
}
