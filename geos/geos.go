package geos

// #cgo LDFLAGS: -lgeos_c
// #include "geos.h"
import "C"
import "unsafe"

type GEOSError string
type GEOSGeometry *C.GEOSGeometry
type STRtree *C.STRtree

// Note: can't use components because GEOS_VERSION_PATCH may be int or string-like
const GEOSVersion string = C.GEOS_VERSION

func (e GEOSError) Error() string {
	return string(e)
}

// Export callback to be able to call from C.
//export geos_errorMessageHandlerCallback
func geos_errorMessageHandlerCallback(message *C.char, userdata unsafe.Pointer) {
	errP := (*error)(userdata)
	*errP = GEOSError(C.GoString(message))
}
