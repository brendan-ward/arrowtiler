package geos

// #cgo LDFLAGS: -lgeos_c
// #include "geos.h"
import "C"
import (
	"fmt"
	"sort"
	"strings"
	"unsafe"

	"github.com/brendan-ward/arrowtiler/tiles"
)

// GeometryArray holds GEOS Geometry pointers (to CGO objects) and a tree (STRtree)
// that is created as part of Query() calls.
// GeometryArray must be manually freed using Release()
type GeometryArray struct {
	geometries []GEOSGeometry
	tree       STRtree
}

// Release GEOS Geometry objects
func (g *GeometryArray) Release() {
	freeGeometries(g.geometries)
	if g.tree != nil {
		C.destroy_tree(g.tree)
	}

	// clear out previous references
	*g = GeometryArray{}
}

func freeGeometries(geometries []GEOSGeometry) {
	if len(geometries) > 0 {
		C.free_geometry_array((**C.GEOSGeometry)(&(geometries[0])), (C.size_t)(len(geometries)))
	}
}

// newGeometryArrayFromGEOS creates a new GeometryArray from a pointer to a C
// array of GEOS geometries and a size.
// Caller must free the C array.
func newGeometryArrayFromGEOS(ptr *GEOSGeometry, size int) *GeometryArray {
	// copy from C array to Go slice (C array must be freed by caller)
	cArr := unsafe.Slice((**C.GEOSGeometry)(ptr), size)
	geometries := make([]GEOSGeometry, size)
	for i := 0; i < size; i++ {
		if cArr[i] != nil {
			geometries[i] = (GEOSGeometry)(unsafe.Pointer(cArr[i]))
		}
	}

	g := &GeometryArray{
		geometries: geometries,
	}

	return g
}

// Create a new GeometryArray from a slice of Geometry Well-Known Text strings.
// The GeometryArray must be freed manually be calling Release().
func NewGeometryArrayFromWKT(wkts []string) (*GeometryArray, error) {
	size := len(wkts)
	// copy from Go strings into C char arrays
	buffer := make([](*C.char), size)
	for i := 0; i < size; i++ {
		buffer[i] = C.CString(wkts[i])
	}

	// char arrays must be deallocated
	defer func() {
		for i := 0; i < size; i++ {
			C.free(unsafe.Pointer(buffer[i]))
		}
	}()

	var ptr *GEOSGeometry = (*GEOSGeometry)(C.from_wkt((**C.char)(&(buffer[0])), (C.size_t)(size)))
	if ptr == nil {
		// TODO: check GEOS error
		return nil, fmt.Errorf("could not parse WKTs")
	}
	defer C.free(unsafe.Pointer(ptr))

	return newGeometryArrayFromGEOS(ptr, size), nil
}

// Create a new GeometryArray from a slice of Geometry Well-Known Binary byte slices.
// The GeometryArray must be freed manually be calling Release().
func NewGeometryArrayFromWKB(wkbs [][]byte) (*GeometryArray, error) {
	size := len(wkbs)
	counts := make([](C.size_t), size)
	// copy from Go strings into C uchar arrays
	buffer := make([](*C.uchar), size)
	for i := 0; i < size; i++ {
		buffer[i] = (*C.uchar)(C.CBytes(wkbs[i]))
		counts[i] = C.size_t(len(wkbs[i]))
	}

	// char arrays must be deallocated
	defer func() {
		for i := 0; i < size; i++ {
			C.free(unsafe.Pointer(buffer[i]))
		}
	}()

	var ptr *GEOSGeometry = (*GEOSGeometry)(C.from_wkb((**C.uchar)(&(buffer[0])), (*C.size_t)(&(counts[0])), C.size_t(len(wkbs))))
	if ptr == nil {
		// TODO: check GEOS error
		return nil, fmt.Errorf("could not parse WKTs")
	}
	defer C.free(unsafe.Pointer(ptr))

	return newGeometryArrayFromGEOS(ptr, len(wkbs)), nil
}

// ToWKT writes the GEOS Geometries using Well-Known Text, according to the
// specified decimal precision.
func (g *GeometryArray) ToWKT(precision int) ([]string, error) {
	size := len(g.geometries)
	if size == 0 {
		return nil, nil
	}

	ptr := C.to_wkt((**C.GEOSGeometry)(&(g.geometries[0])), C.size_t(size), (C.int)(precision))
	if ptr == nil {
		// TODO: check GEOS error
		return nil, fmt.Errorf("could not write to WKT")
	}

	cArr := unsafe.Slice((**C.char)(ptr), size)
	defer func() {
		for i := 0; i < size; i++ {
			if cArr[i] != nil {
				C.free(unsafe.Pointer(cArr[i]))
			}
		}
	}()

	out := make([]string, size)
	for i := 0; i < size; i++ {
		out[i] = C.GoString(cArr[i])
	}

	return out, nil
}

func (g *GeometryArray) Size() int {
	if g == nil {
		return 0
	}

	return len(g.geometries)
}

func (g *GeometryArray) String() string {
	if len(g.geometries) == 0 {
		return ""
	}

	truncate := 60

	wkts, err := g.ToWKT(2)
	if err != nil {
		panic(err)
	}
	var b strings.Builder

	b.WriteString("[")

	for i := 0; i < len(wkts); i++ {
		b.WriteString("<")
		if len(wkts[i]) > truncate {
			b.WriteString(wkts[i][:truncate-3] + "...")
		} else {
			b.WriteString(wkts[i])
		}
		b.WriteString(">")
		if i < len(wkts) {
			b.WriteString(", ")
		}
	}
	b.WriteString("]")
	return b.String()
}

func (g *GeometryArray) TotalBounds() ([4]float64, error) {
	if g == nil {
		panic("GeometryArray not initialized")
	}

	bounds := [4]float64{}

	if C.get_total_bounds((**C.GEOSGeometry)(&(g.geometries[0])),
		(C.size_t)(len(g.geometries)), (*C.double)(&bounds[0]), (*C.double)(&bounds[1]), (*C.double)(&bounds[2]), (*C.double)(&bounds[3])) == 0 {
		return bounds, fmt.Errorf("could not calculate outer bounds of GeometryArray")
	}

	return bounds, nil
}

func (g *GeometryArray) createTree() {
	g.tree = C.create_tree((**C.GEOSGeometry)(&(g.geometries[0])), C.size_t(len(g.geometries)))
	if g.tree == nil {
		panic("could not create tree for GeometryArray")
	}
}

// Query returns a slice of integer indexes into GeometryArray that overlap with
// the bounds defined by xmin, ymin, xmax, ymax.
// Will return nil if there are no results.
func (g *GeometryArray) Query(xmin, ymin, xmax, ymax float64) ([]int, error) {
	if g == nil {
		panic("GeometryArray not initialized")
	}

	if g.tree == nil {
		g.createTree()
	}

	var cArr *C.uint32_t
	var cSize C.size_t
	ret := int(C.query_tree(g.tree, C.double(xmin), C.double(ymin), C.double(xmax), C.double(ymax), (**C.uint32_t)(&cArr), (*C.size_t)(&cSize)))
	if ret != 1 {
		return nil, fmt.Errorf("failed during query of tree")
	}

	defer C.free(unsafe.Pointer(cArr))

	size := int(cSize)
	values := unsafe.Slice((*C.uint32_t)(cArr), size)

	// copy values from uint32_t to int
	indexes := make([]int, size)
	for i := 0; i < size; i++ {
		indexes[i] = (int)(values[i])
	}

	// results are in tree-traversal order; put them into incremental order
	sort.Ints(indexes)

	return indexes, nil
}

// Return a new GeometryArray with coordinates projected to Mercator.
// GeometryCollections are not supported.
// Geometries may be null if outside Mercator world bounds.
func (g *GeometryArray) ToMercator() (*GeometryArray, error) {
	var ptr *GEOSGeometry = (*GEOSGeometry)(C.project_to_mercator((**C.GEOSGeometry)(&(g.geometries[0])), C.size_t(len(g.geometries))))
	if ptr == nil {
		// TODO: check GEOS error
		return nil, fmt.Errorf("could not project to Mercator")
	}
	defer C.free(unsafe.Pointer(ptr))

	return newGeometryArrayFromGEOS(ptr, len(g.geometries)), nil
}

// Return a new GeometryArray with coordinates projected to Mercator, and
// release the previous array.
func (g *GeometryArray) ToMercatorInPlace() error {
	newGeoms, err := g.ToMercator()
	if err != nil {
		return err
	}

	g.Release()

	*g = *newGeoms

	return err
}

// Take creates a new GeometryArray by taking geometries from the GeometryArray
// specified by integer indexes.  Out of bounds indexes will cause a panic.
// The new GeometryArray points to the same underlying geometries.
// TODO: this may fail badly if the master array from which these are taken is released first!
// Should we clone geometries instead!
func (g *GeometryArray) Take(indexes []int) *GeometryArray {
	geometries := make([]GEOSGeometry, len(indexes))
	for i, index := range indexes {
		geometries[i] = g.geometries[index]
	}

	return &GeometryArray{
		geometries: geometries,
		tree:       nil,
	}
}

// Returns integer of indexes into original array for geometries encoded into the tile,
// MVT geometry type, and geometries encoded to MVT uint32 commands and coordinates.
// Input GeometryArray must already be in Mercator coordinates.
func (g *GeometryArray) ToTile(t *tiles.TileID) ([]int, []byte, [][]uint32, error) {
	if g == nil {
		panic("GeometryArray not initialized")
	}

	if len(g.geometries) == 0 {
		return nil, nil, nil, nil
	}

	extent := 4096 // hardcode default for now

	// first figure out if there are any geometries in tile
	if g.tree == nil {
		g.createTree()
	}

	xmin, ymin, xmax, ymax := t.MercatorBounds()
	hits, err := g.Query(xmin, ymin, xmax, ymax)
	if err != nil {
		return nil, nil, nil, err
	}

	size := len(hits)

	if size == 0 {
		return nil, nil, nil, nil
	}

	inGeoms := make([]GEOSGeometry, size)
	for i := 0; i < size; i++ {
		inGeoms[i] = g.geometries[hits[i]]
	}

	var ptr *GEOSGeometry = (*GEOSGeometry)(C.clip_project_to_tile((**C.GEOSGeometry)(&(inGeoms[0])), C.size_t(size), C.double(xmin), C.double(ymin), C.double(xmax), C.double(ymax), C.uint16_t(extent)))

	if ptr == nil {
		return nil, nil, nil, fmt.Errorf("could not extract GeometryArray in tile %v", t)
	}
	defer C.free(unsafe.Pointer(ptr))

	// copy from C array to Go slice
	indexes := make([]int, 0, size+1)
	geomsToEncode := make([]GEOSGeometry, size)

	cArr := unsafe.Slice((**C.GEOSGeometry)(ptr), size)
	for i := 0; i < size; i++ {
		if cArr[i] != nil {
			indexes = append(indexes, hits[i])
			geomsToEncode[i] = (GEOSGeometry)(unsafe.Pointer(cArr[i]))
		}
	}
	// make sure that geometries are released
	defer freeGeometries(geomsToEncode)

	// TODO: additional simplification steps would occur here

	// encode to MVT geometries
	types, encodedGeoms, err := encodeMVTGeometries(geomsToEncode)
	if err != nil {
		return nil, nil, nil, err
	}

	return indexes, types, encodedGeoms, nil
}

func encodeMVTGeometries(geometries []GEOSGeometry) ([]byte, [][]uint32, error) {
	if len(geometries) == 0 {
		return nil, nil, fmt.Errorf("cannot encode empty array of geometries")
	}

	size := len(geometries)

	var cTypesArrPtr *C.uchar
	var cArrPtr **C.uint32_t
	var cSizes *C.size_t

	ret := (int)(C.encode_geometries((**C.GEOSGeometry)(&(geometries[0])), (C.size_t)(size), (**C.uchar)(&cTypesArrPtr), (***C.uint32_t)(&cArrPtr), (**C.size_t)(&cSizes)))
	if ret != 1 {
		return nil, nil, fmt.Errorf("encode geometries failed")
	}

	// free subarrays and containing arrays
	defer C.free_uint32_subarrays((**C.uint32_t)(cArrPtr), (*C.size_t)(cSizes), (C.size_t)(size))
	defer C.free(unsafe.Pointer(cTypesArrPtr))

	cTypes := unsafe.Slice((*C.uchar)(cTypesArrPtr), size)
	cArrs := unsafe.Slice((**C.uint32_t)(cArrPtr), size)
	sizes := unsafe.Slice((*C.size_t)(cSizes), size)

	types := make([]byte, size)
	buffers := make([][]uint32, size)

	for i := 0; i < size; i++ {
		if cArrs[i] == nil {
			buffers[i] = nil
			continue
		}

		types[i] = byte(cTypes[i])
		buffers[i] = make([]uint32, int(sizes[i]))

		cValues := unsafe.Slice((*C.uint32_t)(cArrs[i]), (int)(sizes[i]))
		for j, v := range cValues {
			buffers[i][j] = uint32(v)
		}
	}

	return types, buffers, nil
}
