package mvt

import (
	"encoding/binary"
	"math"
	"unsafe"
)

// Byte identifiers of fields in MVT protobuf
// These are encoded as: (field ID << 3) | field type
// see https://developers.google.com/protocol-buffers/docs/encoding
const (
	// Feature fields
	FEATURE_ID_FIELD        byte = 8  // (1 << 3) | 0
	FEATURE_TAGS_FIELD      byte = 18 // (2 << 3) | 2
	FEATURE_GEOM_TYPE_FIELD byte = 24 // (3 << 3) | 0
	FEATURE_GEOMETRY_FIELD  byte = 34 // (4 << 3) | 2

	// Layer fields
	LAYER_VERSION_FIELD  byte = 120 // (15 << 3) | 0
	LAYER_NAME_FIELD     byte = 10  // (1 << 3) | 2
	LAYER_FEATURES_FIELD byte = 18  // (2 << 3) | 2
	LAYER_KEY_FIELD      byte = 26  // (3 << 3) | 2
	LAYER_VALUE_FIELD    byte = 34  // (4 << 3) | 2
	LAYER_EXTENT_FIELD   byte = 40  // (5 << 3) | 0

	// Value fields
	VALUE_STRING_FIELD    byte = 10 // (1 << 3) | 2
	VALUE_FLOAT32_FIELD   byte = 21 // (2 << 3) | 5
	VALUE_FLOAT64_FIELD   byte = 24 // (3 << 3) | 1
	VALUE_INT64_FIELD     byte = 32 // (4 << 3) | 0
	VALUE_UVARINT64_FIELD byte = 40 // (5 << 3) | 0
	VALUE_VARINT64_FIELD  byte = 48 // (6 << 3) | 0
	VALUE_BOOL_FIELD      byte = 56 // (6 << 3) | 0

	// Tile fields
	TILE_LAYERS_FIELD byte = 26 // (3 << 3) | 2
)

// EncodeString returns a byte buffer with the a Uvarint of the length of the string
// followed by the string value
func EncodeString(v string) []byte {
	return append(EncodeUvarint(uint64(len(v))), v...)
}

// EncodeBytes returns a byte buffer with the a Uvarint of the length of the byte slice
// followed by the bytes
func EncodeBytes(v []byte) []byte {
	return append(EncodeUvarint(uint64(len(v))), v...)
}

// EncodeUvarint returns a byte buffer with the bytes representing the Uvarint value
// of v
func EncodeUvarint(v uint64) []byte {
	buffer := make([]byte, binary.MaxVarintLen64)
	size := binary.PutUvarint(buffer, v)
	return buffer[:size]
}

// EncodeUvarint returns a byte buffer with the bytes representing the Varint value
// of v
func EncodeVarint(v int64) []byte {
	buffer := make([]byte, binary.MaxVarintLen64)
	size := binary.PutVarint(buffer, v)
	return buffer[:size]
}

func EncodeUint32(v uint32) []byte {
	// assumes only running on little endian systems
	return (*[4]byte)(unsafe.Pointer(&v))[:]
}

// wrapValue wraps the byte slice as a protobuf Value
func wrapValue(b []byte, specificType byte) []byte {
	// add specific type before value bytes
	value := append([]byte{specificType}, b...)
	valueLength := EncodeUvarint(uint64(len(b) + 1))

	// reserve enough space for length and type
	buffer := make([]byte, 1, len(b)+len(valueLength)+2)
	buffer[0] = LAYER_VALUE_FIELD
	buffer = append(buffer, valueLength...)
	return append(buffer, value...)
}

// EncodeByteValue encodes a byte slice as a Value for the MVT protobuf
func EncodeByteValue(v []byte) []byte {
	return wrapValue(EncodeBytes(v), VALUE_STRING_FIELD)
}

// EncodeByteValue encodes a string as a Value for the MVT protobuf
func EncodeStringValue(v string) []byte {
	return wrapValue(EncodeString(v), VALUE_STRING_FIELD)
}

// EncodeByteValue encodes a uint64 as a Value (of type Uvarint) for the MVT protobuf
func EncodeUint64Value(v uint64) []byte {
	return wrapValue(EncodeUvarint(v), VALUE_UVARINT64_FIELD)
}

// EncodeByteValue encodes a int64 as a Value (of type Varint) for the MVT protobuf
func EncodeInt64Value(v int64) []byte {
	return wrapValue(EncodeVarint(v), VALUE_VARINT64_FIELD)
}

// EncodeByteValue encodes a float32 as a Value for the MVT protobuf
func EncodeFloat32Value(v float32) []byte {
	buffer := make([]byte, 4)
	binary.LittleEndian.PutUint32(buffer, math.Float32bits(v))
	return wrapValue(buffer, VALUE_FLOAT32_FIELD)
}

// EncodeByteValue encodes a float36 as a Value for the MVT protobuf
func EncodeFloat64Value(v float64) []byte {
	buffer := make([]byte, 8)
	binary.LittleEndian.PutUint64(buffer, math.Float64bits(v))
	return wrapValue(buffer, VALUE_FLOAT64_FIELD)
}

// EncodeByteValue encodes a bool as a Value for the MVT protobuf
func EncodeBoolValue(v bool) []byte {
	buffer := make([]byte, 1)
	if v {
		buffer[0] = 1
	} else {
		buffer[0] = 0
	}
	return wrapValue(buffer, VALUE_BOOL_FIELD)
}

// EncodeKey returns a buffer representing the key type, followed by the length
// of the string, followed by the bytes representing the string
func EncodeKey(key string) []byte {
	s := EncodeString(key)
	buffer := make([]byte, 0, len(s)+2)
	buffer = append(buffer, LAYER_KEY_FIELD)
	return append(buffer, EncodeString(key)...)
}
