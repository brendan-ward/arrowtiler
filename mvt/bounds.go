package mvt

import "math"

var RAD2DEG float64 = 180 / math.Pi
var EARTH_RADIUS float64 = 6378137.0

func MercatorBoundsToGeoBounds(b [4]float64) [4]float64 {
	return [4]float64{
		b[0] * RAD2DEG / EARTH_RADIUS,
		((math.Pi * 0.5) - 2.0*math.Atan(math.Exp(-b[1]/EARTH_RADIUS))) * RAD2DEG,
		b[2] * RAD2DEG / EARTH_RADIUS,
		((math.Pi * 0.5) - 2.0*math.Atan(math.Exp(-b[3]/EARTH_RADIUS))) * RAD2DEG,
	}
}
