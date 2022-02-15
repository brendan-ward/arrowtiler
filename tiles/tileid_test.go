package tiles

import (
	"math"
	"testing"
)

func closeEnough(a, b, tolerance float64) bool {
	return math.Abs(a-b) < tolerance
}

func Test_GeoToTile(t *testing.T) {
	tests := []struct {
		zoom uint16
		lon  float64
		lat  float64
		x    uint32
		y    uint32
	}{
		{zoom: 0, lon: 0, lat: 0, x: 0, y: 0},
		{zoom: 1, lon: 0, lat: 0, x: 1, y: 1},
		{zoom: 1, lon: 0, lat: -0.1, x: 1, y: 1},
		{zoom: 1, lon: 0, lat: 0.1, x: 1, y: 0},
		{zoom: 1, lon: -180, lat: -90, x: 0, y: 1},
		{zoom: 1, lon: -180, lat: 90, x: 0, y: 0},
		{zoom: 1, lon: -20, lat: 20, x: 0, y: 0},
		{zoom: 1, lon: 0, lat: 0, x: 1, y: 1},
		{zoom: 4, lon: -20, lat: 20, x: 7, y: 7},
	}

	for _, tc := range tests {
		tile := GeoToTile(tc.zoom, tc.lon, tc.lat)
		if tile.Zoom != tc.zoom {
			t.Errorf("zoom: %v (%f, %f) => %v | zoom: %v not expected value\n", tc.zoom, tc.lon, tc.lat, tile, tile.Zoom)
		}
		if tile.X != tc.x {
			t.Errorf("zoom: %v (%f, %f) => %v | tile x not expected value: %v\n", tc.zoom, tc.lon, tc.lat, tile, tc.x)
		}
		if tile.Y != tc.y {
			t.Errorf("zoom: %v (%f, %f) => %v | tile y not expected value: %v\n", tc.zoom, tc.lon, tc.lat, tile, tc.y)
		}
	}
}

func Test_TileRange(t *testing.T) {
	tests := []struct {
		zoom    uint16
		bounds  [4]float64
		minTile *TileID
		maxTile *TileID
	}{
		{zoom: 0, bounds: [4]float64{-180, -90, 180, 90}, minTile: &TileID{0, 0, 0}, maxTile: &TileID{0, 0, 0}},
		{zoom: 1, bounds: [4]float64{-180, -90, 180, 90}, minTile: &TileID{1, 0, 0}, maxTile: &TileID{1, 1, 1}},
		{zoom: 1, bounds: [4]float64{-180, -90, 0, 90}, minTile: &TileID{1, 0, 0}, maxTile: &TileID{1, 0, 1}},
		{zoom: 1, bounds: [4]float64{-180, -90, -170, -80}, minTile: &TileID{1, 0, 1}, maxTile: &TileID{1, 0, 1}},
		{zoom: 4, bounds: [4]float64{-100, -20, -20, 20}, minTile: &TileID{4, 3, 7}, maxTile: &TileID{4, 7, 8}},
		{zoom: 4, bounds: [4]float64{-118.22826385, 17.93906593, -65.33223724, 48.99859619}, minTile: &TileID{4, 2, 5}, maxTile: &TileID{4, 5, 7}},
		{zoom: 4, bounds: [4]float64{-1e-6, -1e-6, 1e-6, 1e-6}, minTile: &TileID{4, 7, 7}, maxTile: &TileID{4, 8, 8}},
	}

	for _, tc := range tests {
		// convert to Mercator
		xmin, ymin := GeoToMercator(tc.bounds[0], tc.bounds[1])
		xmax, ymax := GeoToMercator(tc.bounds[2], tc.bounds[3])

		minTile, maxTile := TileRange(tc.zoom, [4]float64{xmin, ymin, xmax, ymax})

		if minTile.Zoom != tc.zoom {
			t.Errorf("minTile %v | zoom not expected value: %v\n", minTile, tc.zoom)
		}
		if minTile.X != tc.minTile.X {
			t.Errorf("minTile %v | x not expected value: %v\n", minTile, tc.minTile)
		}
		if minTile.Y != tc.minTile.Y {
			t.Errorf("minTile %v | y not expected value: %v\n", minTile, tc.minTile)
		}
		if maxTile.Zoom != tc.zoom {
			t.Errorf("maxTile %v | zoom not expected value: %v\n", maxTile, tc.zoom)
		}
		if maxTile.X != tc.maxTile.X {
			t.Errorf("maxTile %v | x not expected value: %v\n", maxTile, tc.maxTile)
		}
		if maxTile.Y != tc.maxTile.Y {
			t.Errorf("maxTile %v | y not expected value: %v\n", maxTile, tc.maxTile)
		}
	}
}

func Test_GeoBounds(t *testing.T) {
	tests := []struct {
		zoom uint16
		x    uint32
		y    uint32
		xmin float64
		ymin float64
		xmax float64
		ymax float64
	}{
		{zoom: 0, x: 0, y: 0, xmin: -180.0, ymin: -85.051129, xmax: 180.0, ymax: 85.051129},
		{zoom: 1, x: 1, y: 1, xmin: 0, ymin: -85.051129, xmax: 180, ymax: 0},
		{zoom: 10, x: 20, y: 30, xmin: -172.968750, ymin: 84.016022, xmax: -172.617188, ymax: 84.052561},
	}

	tol := 1e-6

	for _, tc := range tests {
		tile := NewTileID(tc.zoom, tc.x, tc.y)
		xmin, ymin, xmax, ymax := tile.GeoBounds()
		if !closeEnough(xmin, tc.xmin, tol) {
			t.Errorf("%v | xmin: %f does not match expected value: %f\n", tile, xmin, tc.xmin)
		}
		if !closeEnough(ymin, tc.ymin, tol) {
			t.Errorf("%q | ymin: %f does not match expected value: %f\n", tile, ymin, tc.ymin)
		}
		if !closeEnough(xmax, tc.xmax, tol) {
			t.Errorf("%q | xmax: %f does not match expected value: %f\n", tile, xmax, tc.xmax)
		}
		if !closeEnough(ymax, tc.ymax, tol) {
			t.Errorf("%q | ymax: %f does not match expected value: %f\n", tile, ymax, tc.ymax)
		}
	}
}

func Test_MercatorBounds(t *testing.T) {
	tests := []struct {
		zoom uint16
		x    uint32
		y    uint32
		xmin float64
		ymin float64
		xmax float64
		ymax float64
	}{
		{zoom: 0, x: 0, y: 0, xmin: -20037508.342789, ymin: -20037508.342789, xmax: 20037508.342789, ymax: 20037508.342789},
		{zoom: 1, x: 1, y: 1, xmin: 0, ymin: -20037508.342789, xmax: 20037508.342789, ymax: 0},
		{zoom: 10, x: 20, y: 30, xmin: -19254793.173149, ymin: 18824299.829847, xmax: -19215657.414667, ymax: 18863435.588329},
	}

	tol := 1e-6

	for _, tc := range tests {
		tile := NewTileID(tc.zoom, tc.x, tc.y)
		xmin, ymin, xmax, ymax := tile.MercatorBounds()
		if !closeEnough(xmin, tc.xmin, tol) {
			t.Errorf("tile: %v | xmin: %f does not match expected value: %f\n", tile, xmin, tc.xmin)
		}
		if !closeEnough(ymin, tc.ymin, tol) {
			t.Errorf("tile: %q | ymin: %f does not match expected value: %f\n", tile, ymin, tc.ymin)
		}
		if !closeEnough(xmax, tc.xmax, tol) {
			t.Errorf("tile: %q | xmax: %f does not match expected value: %f\n", tile, xmax, tc.xmax)
		}
		if !closeEnough(ymax, tc.ymax, tol) {
			t.Errorf("tile: %q | ymax: %f does not match expected value: %f\n", tile, ymax, tc.ymax)
		}
	}
}
