package tiles

import (
	"fmt"
	"math"
)

var RE float64 = 6378137.0
var ORIGIN = RE * math.Pi
var CE float64 = 2.0 * ORIGIN
var DEG2RAD float64 = math.Pi / 180.0

// WebMercator tile, numbered starting from upper left
type TileID struct {
	Zoom uint8
	X    uint32
	Y    uint32
}

func NewTileID(zoom uint8, x uint32, y uint32) *TileID {
	return &TileID{zoom, x, y}
}

func GeoToMercator(lon float64, lat float64) (x float64, y float64) {
	// truncate incoming values to world bounds
	lon = math.Min(math.Max(lon, -180), 180)
	lat = math.Min(math.Max(lat, -85.051129), 85.051129)

	x = lon * ORIGIN / 180.0
	y = RE * math.Log(math.Tan((math.Pi*0.25)+(0.5*DEG2RAD*lat)))
	return
}

// GeoToTile calculates the tile x,y at zoom that contains longitude, latitude
func GeoToTile(zoom uint8, x float64, y float64) *TileID {
	z2 := 1 << zoom
	zoomFactor := float64(z2)
	eps := 1e-14

	// truncate incoming values to world bounds
	x = math.Min(math.Max(x, -180), 180)
	y = math.Min(math.Max(y, -85.051129), 85.051129)

	var tileX uint32
	var tileY uint32

	x = math.Max(x/360.0+0.5, 0.0)
	if x >= 1 {
		tileX = uint32(z2 - 1)
	} else {
		tileX = uint32(math.Floor((x + eps) * zoomFactor))
	}

	y = math.Sin(y * math.Pi / 180)
	y = 0.5 - 0.25*math.Log((1.0+y)/(1.0-y))/math.Pi
	if y >= 1 {
		tileY = uint32(z2 - 1)
	} else {
		tileY = uint32((y + eps) * zoomFactor)
	}

	return &TileID{
		Zoom: zoom,
		X:    tileX,
		Y:    tileY,
	}
}

// TileRange calculates the min tile x, min tile y, max tile x, max tile y tile
// range for Mercator coordinates xmin, ymin, xmax, ymax at a given zoom level.
// Assumes bounds have already been clipped to Mercator world bounds.
func TileRange(zoom uint8, bounds [4]float64) (*TileID, *TileID) {
	z2 := 1 << zoom
	zoomFactor := float64(z2)
	origin := -ORIGIN
	eps := 1.0e-11

	xmin := math.Min(math.Max(math.Floor(((bounds[0]-origin)/CE)*zoomFactor), 0), zoomFactor-1)
	// ymin isn't right yet, spilling over
	ymin := math.Min(math.Max(math.Floor(((1.0-(((bounds[1]-origin)/CE)+eps))*zoomFactor)), 0), zoomFactor-1)
	xmax := math.Min(math.Max(math.Floor((((bounds[2]-origin)/CE)-eps)*zoomFactor), 0), zoomFactor-1)
	ymax := math.Min(math.Max(math.Floor((1.0-((bounds[3]-origin)/CE))*zoomFactor), 0), zoomFactor-1)

	// tiles start in upper left, flip y values
	minTile := &TileID{Zoom: zoom, X: uint32(xmin), Y: uint32(ymax)}
	maxTile := &TileID{Zoom: zoom, X: uint32(xmax), Y: uint32(ymin)}

	return minTile, maxTile
}

func (t *TileID) String() string {
	return fmt.Sprintf("Tile(zoom: %v, x: %v, y:%v)", t.Zoom, t.X, t.Y)
}

func (t *TileID) GeoBounds() (float64, float64, float64, float64) {
	z2 := 1 << t.Zoom
	zoomFactor := (float64)(z2)
	x := (float64)(t.X)
	y := (float64)(t.Y)

	xmin := x/zoomFactor*360.0 - 180.0
	ymin := math.Atan(math.Sinh(math.Pi*(1.0-2.0*((y+1.0)/zoomFactor)))) * (180.0 / math.Pi)
	xmax := (x+1.0)/zoomFactor*360.0 - 180.0
	ymax := math.Atan(math.Sinh(math.Pi*(1.0-2.0*y/zoomFactor))) * (180.0 / math.Pi)

	return xmin, ymin, xmax, ymax
}

func (t *TileID) MercatorBounds() (float64, float64, float64, float64) {
	z2 := 1 << t.Zoom
	tileSize := CE / (float64)(z2)
	xmin := (float64)(t.X)*tileSize - CE/2.0
	xmax := xmin + tileSize
	ymax := CE/2 - (float64)(t.Y)*tileSize
	ymin := ymax - tileSize

	return xmin, ymin, xmax, ymax
}
