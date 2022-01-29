package mvt

import (
	"fmt"
	"os"

	"github.com/apache/arrow/go/arrow"
	"github.com/apache/arrow/go/arrow/array"
	"github.com/apache/arrow/go/arrow/ipc"
	"github.com/brendan-ward/arrowtiler/geos"
	"github.com/brendan-ward/arrowtiler/tiles"
)

var EXTENT uint16 = 4096 // vector tile default extent

type LayerInfo struct {
	Name        string            `json:"id"`
	Description string            `json:"description"`
	Minzoom     uint32            `json:"minzoom"`
	Maxzoom     uint32            `json:"maxzoom"`
	Fields      map[string]string `json:"fields"`
}

// FeatureTable is a data structure for holding feature information as a
// GeometryArray and array of encoded attribute data
type FeatureTable struct {
	ids        []uint64
	geometries *geos.GeometryArray
	columns    []*ByteColumn
}

func ReadFeather(path string, geomColName string, idColName string) (*FeatureTable, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	r, err := ipc.NewFileReader(f)
	if err != nil {
		return nil, err
	}
	defer r.Close()

	// read records into a simple Arrow table
	records := make([]array.Record, r.NumRecords())
	for i := 0; i < r.NumRecords(); i++ {
		records[i], err = r.RecordAt(i)
		if err != nil {
			return nil, err
		}
	}
	t := array.NewTableFromRecords(r.Schema(), records)

	i := 0
	var col *array.Chunked

	// Extract id field (optional)
	var ids []uint64
	idColIdx := -1
	if idColName != "" {
		if idColName == "__auto__" {
			fmt.Println("Autogenerating id field")

			ids = make([]uint64, t.NumRows())
			for i = 0; i < len(ids); i++ {
				ids[i] = uint64(i)
			}
		} else {
			fmt.Printf("Using '%v' column for id\n", idColName)

			idColIdxs := r.Schema().FieldIndices(idColName)
			if idColIdxs == nil {
				return nil, fmt.Errorf("'%s' column must be present to use as id", idColName)
			}
			idColIdx = idColIdxs[0]
			col = t.Column(idColIdx).Data()
			ids = make([]uint64, t.NumRows())

			i = 0
			for _, chunk := range col.Chunks() {
				for j := 0; j < chunk.Len(); j++ {
					id, err := getFeatureID(chunk, i)
					if err != nil {
						return nil, err
					}
					ids[i] = id
					i++
				}
			}
		}
	}

	// Extract geometries
	geomColIdxs := r.Schema().FieldIndices(geomColName)
	if geomColIdxs == nil {
		return nil, fmt.Errorf("'%v' column must be present", geomColName)
	}
	geomColIdx := geomColIdxs[0]

	col = t.Column(geomColIdx).Data()

	wkbs := make([][]byte, t.NumRows())
	i = 0
	for _, chunk := range col.Chunks() {
		c := chunk.(*array.Binary)
		for j := 0; j < chunk.Len(); j++ {
			wkbs[i] = c.Value(j)
			i++
		}
	}

	geometries, err := geos.NewGeometryArrayFromWKB(wkbs)
	if err != nil {
		return nil, err
	}

	// project to Mercator
	err = geometries.ToMercatorInPlace()
	if err != nil {
		return nil, err
	}

	// encode non-geometry, non-id columns
	var columns []*ByteColumn

	fields := r.Schema().Fields()
	// columns := make([]string, len(fields))
	for fieldIdx, field := range fields {
		// skip geometry and id
		if (fieldIdx == geomColIdx) || idColIdx != -1 && fieldIdx == idColIdx {
			continue
		}

		colType, err := getColumnType(field)
		if err != nil {
			return nil, err
		}

		column := NewEmptyByteColumn(field.Name, colType, int(t.NumRows()))

		i = 0
		for _, chunk := range t.Column(fieldIdx).Data().Chunks() {
			if err != nil {
				return nil, err
			}

			for j := 0; j < chunk.Len(); j++ {
				encodedValue, err := getEncodedValue(chunk, j)
				if err != nil {
					return nil, err
				}
				column.SetValue(i, encodedValue)

				i++
			}
		}
		if columns == nil {
			columns = make([]*ByteColumn, 1)
			columns[0] = column
		} else {
			columns = append(columns, column)
		}
	}

	return &FeatureTable{
		ids:        ids,
		geometries: geometries,
		columns:    columns,
	}, nil
}

func (t *FeatureTable) Geometry() *geos.GeometryArray {
	return t.geometries
}

func (t *FeatureTable) Column(i int) *ByteColumn {
	return t.columns[i]
}

func (t *FeatureTable) NumCols() int {
	return len(t.columns)
}

func (t *FeatureTable) Size() int {
	if t == nil {
		return 0
	}
	return t.geometries.Size()
}

func getColumnType(field arrow.Field) (string, error) {
	switch field.Type.ID() {
	case arrow.STRING, arrow.BINARY:
		return "String", nil
	case arrow.INT8, arrow.INT16, arrow.INT32, arrow.INT64, arrow.UINT8, arrow.UINT16, arrow.UINT32, arrow.UINT64:
		return "Number", nil
	case arrow.BOOL:
		return "Boolean", nil
	case arrow.FLOAT32, arrow.FLOAT64:
		return "Float", nil
	default:
		return "", fmt.Errorf("type not supported: %v", field.Type.ID())
	}
}

func getFeatureID(chunk array.Interface, i int) (uint64, error) {
	// TODO: raise error on null value
	switch chunk.DataType().ID() {
	case arrow.INT8:
		id := chunk.(*array.Int8).Value(i)
		if id < 0 {
			return 0, fmt.Errorf("cannot use column with negative value for id")
		}
		return uint64(id), nil
	case arrow.INT16:
		id := chunk.(*array.Int16).Value(i)
		if id < 0 {
			return 0, fmt.Errorf("cannot use column with negative value for id")
		}
		return uint64(id), nil
	case arrow.INT32:
		id := chunk.(*array.Int32).Value(i)
		if id < 0 {
			return 0, fmt.Errorf("cannot use column with negative value for id")
		}
		return uint64(id), nil
	case arrow.INT64:
		id := chunk.(*array.Int64).Value(i)
		if id < 0 {
			return 0, fmt.Errorf("cannot use column with negative value for id")
		}
		return uint64(id), nil
	case arrow.UINT8:
		return uint64(chunk.(*array.Uint8).Value(i)), nil
	case arrow.UINT16:
		return uint64(chunk.(*array.Uint16).Value(i)), nil
	case arrow.UINT32:
		return uint64(chunk.(*array.Uint32).Value(i)), nil
	case arrow.UINT64:
		return chunk.(*array.Uint64).Value(i), nil
	default:
		return 0, fmt.Errorf("cannot use non-integer column for id")
	}
}

func getEncodedValue(chunk array.Interface, i int) ([]byte, error) {
	// TODO: check for null
	switch chunk.DataType().ID() {
	case arrow.BINARY:
		return EncodeByteValue(chunk.(*array.Binary).Value(i)), nil
	case arrow.STRING:
		return EncodeStringValue(chunk.(*array.String).Value(i)), nil
	case arrow.INT8:
		return EncodeInt64Value(int64(chunk.(*array.Int8).Value(i))), nil
	case arrow.INT16:
		return EncodeInt64Value(int64(chunk.(*array.Int16).Value(i))), nil
	case arrow.INT32:
		return EncodeInt64Value(int64(chunk.(*array.Int32).Value(i))), nil
	case arrow.INT64:
		return EncodeInt64Value(chunk.(*array.Int64).Value(i)), nil
	case arrow.UINT8:
		return EncodeUint64Value(uint64(chunk.(*array.Uint8).Value(i))), nil
	case arrow.UINT16:
		return EncodeUint64Value(uint64(chunk.(*array.Uint16).Value(i))), nil
	case arrow.UINT32:
		return EncodeUint64Value(uint64(chunk.(*array.Uint32).Value(i))), nil
	case arrow.UINT64:
		return EncodeUint64Value(chunk.(*array.Uint64).Value(i)), nil
	case arrow.FLOAT32:
		return EncodeFloat32Value(chunk.(*array.Float32).Value(i)), nil
	case arrow.FLOAT64:
		return EncodeFloat64Value(chunk.(*array.Float64).Value(i)), nil
	case arrow.BOOL:
		return EncodeBoolValue(chunk.(*array.Boolean).Value(i)), nil

	default:
		return nil, fmt.Errorf("no available encoder for %v", chunk.DataType().ID())
	}
}

func (t *FeatureTable) Take(indexes []int) *FeatureTable {
	hasId := t.ids != nil
	var ids []uint64
	if hasId {
		ids = make([]uint64, len(indexes))
		for i := 0; i < len(indexes); i++ {
			index := indexes[i]
			if hasId {
				ids[i] = t.ids[index]
			}
		}
	}

	var columns []*ByteColumn
	if len(t.columns) > 0 {
		columns = make([]*ByteColumn, len(t.columns))
		for i, col := range t.columns {
			columns[i] = col.Take(indexes)
		}
	}

	return &FeatureTable{
		ids:        ids,
		geometries: t.geometries.Take(indexes),
		columns:    columns,
	}
}

func (t *FeatureTable) EncodeToLayer(name string, tile *tiles.TileID) ([]byte, error) {
	// WARNING: assumes that all columns are unique names
	// and that all values are at least empty for their type

	hasId := t.ids != nil

	// project geometries to tile and encode to MVT format
	index, geomTypes, encodedGeoms, err := t.geometries.ToTile(tile)

	if err != nil {
		return nil, err
	}
	if len(index) == 0 {
		return nil, nil
	}

	// track counters used for feature-level tags
	keyIndex := make(map[int]uint32)
	var keyCounter uint32 = 0
	valueIndex := make(map[string]uint32)
	var valueCounter uint32 = 0
	var keyBuffer []byte
	var valueBuffer []byte

	// encode features (tags + geometries per feature)
	featuresBuffer := make([]byte, 0)

	// collect tags and convert geometry uint32 => uvarint
	for i, rowIdx := range index {
		featureBuffer := make([]byte, 0)
		geomBuffer := make([]byte, 0)
		for _, v := range encodedGeoms[i] {
			// convert uint32 => MVT varint
			geomBuffer = append(geomBuffer, EncodeUvarint(uint64(v))...)
		}

		tagsBuffer := make([]byte, 0)
		for colIdx, col := range t.columns {
			value := col.GetValue(rowIdx)

			// TODO: verify correct handling of nil values
			if value != nil {
				if keyIdx, ok := keyIndex[colIdx]; !ok {
					keyBuffer = append(keyBuffer, EncodeKey(col.Name)...)
					tagsBuffer = append(tagsBuffer, EncodeUvarint(uint64(keyCounter))...)
					keyIndex[colIdx] = keyCounter
					keyCounter++
				} else {
					tagsBuffer = append(tagsBuffer, EncodeUvarint(uint64(keyIdx))...)
				}

				// this should leverage fast string key for map without actually
				// converting bytes to string
				if valIdx, ok := valueIndex[string(value)]; !ok {
					valueBuffer = append(valueBuffer, value...)
					tagsBuffer = append(tagsBuffer, EncodeUvarint(uint64(valueCounter))...)
					valueIndex[string(value)] = valueCounter
					valueCounter++
				} else {
					tagsBuffer = append(tagsBuffer, EncodeUvarint(uint64(valIdx))...)
				}
			}
		}

		// feature ID
		if hasId {
			featureBuffer = append(featureBuffer, FEATURE_ID_FIELD)
			featureBuffer = append(featureBuffer, EncodeUvarint(t.ids[rowIdx])...)
		}

		// feature geometry type
		// tippecanoe encodes this before tags
		featureBuffer = append(featureBuffer, FEATURE_GEOM_TYPE_FIELD)
		featureBuffer = append(featureBuffer, geomTypes[i])

		// feature tags
		featureBuffer = append(featureBuffer, FEATURE_TAGS_FIELD)
		featureBuffer = append(featureBuffer, EncodeUvarint(uint64(len(tagsBuffer)))...)
		featureBuffer = append(featureBuffer, tagsBuffer...)

		// feature geometry
		featureBuffer = append(featureBuffer, FEATURE_GEOMETRY_FIELD)
		featureBuffer = append(featureBuffer, EncodeUvarint(uint64(len(geomBuffer)))...)
		featureBuffer = append(featureBuffer, geomBuffer...)

		featuresBuffer = append(featuresBuffer, LAYER_FEATURES_FIELD)
		featuresBuffer = append(featuresBuffer, EncodeUvarint(uint64(len(featureBuffer)))...)
		featuresBuffer = append(featuresBuffer, featureBuffer...)
	}

	// Write the final layer buffer
	buffer := make([]byte, 0)

	// version
	buffer = append(buffer, 120, 2)

	// layer name
	buffer = append(buffer, LAYER_NAME_FIELD)
	buffer = append(buffer, EncodeString(name)...)

	// layer extent
	buffer = append(buffer, LAYER_EXTENT_FIELD)
	buffer = append(buffer, EncodeUvarint(uint64(EXTENT))...)

	// layer tag keys
	buffer = append(buffer, keyBuffer...)

	// layer tag values
	buffer = append(buffer, valueBuffer...)

	// features
	buffer = append(buffer, featuresBuffer...)

	// add total size of layer to front of buffer
	var tmp []byte
	tmp = append(tmp, TILE_LAYERS_FIELD)
	tmp = append(tmp, EncodeUvarint(uint64(len(buffer)))...)
	buffer = append(tmp, buffer...)

	return buffer, nil
}

func (t *FeatureTable) GetLayerInfo(name string, description string, minZoom uint32, maxZoom uint32) (*LayerInfo, error) {
	fields := make(map[string]string)
	for _, col := range t.columns {
		fields[col.Name] = col.Type
	}

	info := &LayerInfo{
		Name:        name,
		Description: description,
		Minzoom:     minZoom,
		Maxzoom:     maxZoom,
		Fields:      fields,
	}

	return info, nil
}
