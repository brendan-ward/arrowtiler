package tiles

import (
	"fmt"
	"math"
)

var CE float64 = 2 * 6378137.0 * math.Pi

// WebMercator tile, numbered starting from upper left
type TileID struct {
	Zoom uint32
	X    uint32
	Y    uint32
}

func NewTileID(zoom, x, y uint32) *TileID {
	return &TileID{zoom, x, y}
}

// GeoToTile calculates the tile x,y at zoom that contains longitude, latitude
func GeoToTile(zoom uint32, x float64, y float64) *TileID {
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
// range for geographic coordinates xmin, ymin, xmax, ymax at a given zoom level
// TODO: convert to use Mercator coordinates instead
func TileRange(zoom uint32, bounds [4]float64) (*TileID, *TileID) {
	eps := 1.0e-11

	minTile := GeoToTile(zoom, bounds[0], bounds[3])
	maxTile := GeoToTile(zoom, bounds[2]-eps, bounds[1]+eps)

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
