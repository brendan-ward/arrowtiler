package tiles

// EncodingParams holds parameters used for encoding tiles
type EncodingConfig struct {
	Extent         uint16
	Buffer         uint16
	Precision      uint8
	Simplification uint8
}

func NewDefaultEncodingConfig() *EncodingConfig {
	return &EncodingConfig{
		Extent:         4096, // pixels, vector tile default
		Buffer:         256,  // pixels, PostGIS default
		Precision:      1,    // pixels
		Simplification: 1,    // pixels
	}
}

func NewEncodingConfig(extent uint16, buffer uint16, precision uint8, simplification uint8) *EncodingConfig {
	return &EncodingConfig{
		Extent:         extent,
		Buffer:         buffer,
		Precision:      precision,
		Simplification: simplification,
	}
}
