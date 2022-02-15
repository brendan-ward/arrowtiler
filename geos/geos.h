
#ifndef GEOS_H
#define GEOS_H

#define GEOS_USE_ONLY_R_API
#include <geos_c.h>
#include <stdint.h>
#include <stdlib.h>

typedef struct {
  size_t n, m;
  uint32_t *a;
} uint32_vec_t;

/* struct that wraps GEOSSTRtree to keep track of input geometry indices */
struct STRtree_t {
  GEOSSTRtree *geos_tree;
  uint32_t *indexes;
};

typedef struct STRtree_t STRtree;

/* struct that returns an array of size_t values; used to pass back array to Go
 * as a return value;
 */
struct Uint32Array_t {
  uint32_t *values;
  size_t size;
};

typedef struct Uint32Array_t Uint32Array;

struct ScaleParams_t {
  double xscale;
  double xoffset;
  double yscale;
  double yoffset;
};

typedef struct ScaleParams_t ScaleParams;

// error handling
enum { GEOSERROR, ERROR, WARNING, SUCCESS };

// error callback from GEOS to go
void geos_error_handler(const char *message, void *userdata);

// Define initialization / finalization macros that wrap calls to GEOS
#define _GEOS_INIT_DEF                                                         \
  char errstate = SUCCESS;                                                     \
  char last_error[1024] = "";                                                  \
  char last_warning[1024] = "";                                                \
  GEOSContextHandle_t ctx

#define _GEOS_INIT                                                             \
  ctx = GEOS_init_r();                                                         \
  GEOSContext_setErrorMessageHandler_r(ctx, geos_error_handler, last_error)

// TODO:
// GEOSContext_setNoticeMessageHandler_r(ctx,
// geos_notice_handler, last_warning)

#define GEOS_INIT                                                              \
  _GEOS_INIT_DEF;                                                              \
  _GEOS_INIT

#define GEOS_FINISH GEOS_finish_r(ctx)

// TODO:
// GEOS_HANDLE_ERR

void free_geometry_array(GEOSGeometry **geoms, size_t count);

GEOSGeometry **from_wkb(unsigned char **wkbs, size_t *counts, size_t count);

GEOSGeometry **from_wkt(char **wkts, size_t count);

char **to_wkt(GEOSGeometry **geometries, size_t count, int precision);

int get_total_bounds(GEOSGeometry **geometries, size_t num_geometries,
                     double *xmin, double *ymin, double *xmax, double *ymax);

void free_uint32_subarrays(uint32_t **values, size_t *sizes, size_t size);

/* STRtree functions */
STRtree *create_tree(GEOSGeometry **geometries, size_t count);

int query_tree(STRtree *tree, double xmin, double ymin, double xmax,
               double ymax, uint32_t **hits, size_t *hit_count);

void destroy_tree(STRtree *tree);

/* Project and tile encoding functions */

GEOSGeometry **project_to_mercator(GEOSGeometry **geometries, size_t count);

GEOSGeometry **clip_project_to_tile(GEOSGeometry **geometries,
                                    size_t num_geometries, double xmin,
                                    double ymin, double xmax, double ymax,
                                    uint16_t extent);

int encode_geometries(GEOSGeometry **geometries, size_t num_geometries,
                      unsigned char **types, uint32_t ***buffer,
                      size_t **sizes);

#endif // GEOS_H