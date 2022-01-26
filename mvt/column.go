package mvt

type colValue []byte

// ByteColumn represents a ByteColumn of MVT encoded value data
type ByteColumn struct {
	Name   string
	Type   string // MVT type
	values []colValue
}

func NewEmptyByteColumn(name string, colType string, size int) *ByteColumn {
	return &ByteColumn{
		Name:   name,
		Type:   colType,
		values: make([]colValue, size),
	}
}

func NewByteColumn(name string, colType string, values []colValue) *ByteColumn {
	return &ByteColumn{
		Name:   name,
		Type:   colType,
		values: values,
	}
}

func (c *ByteColumn) GetValue(i int) []byte {
	return c.values[i]
}

func (c *ByteColumn) SetValue(i int, value []byte) {
	c.values[i] = value
}

func (c *ByteColumn) Size() int {
	return len(c.values)
}

// Take creates a new ByteColumn by taking values from the ByteColumn
// specified by integer indexes.  Out of bounds indexes will cause a panic.
func (c *ByteColumn) Take(indexes []int) *ByteColumn {
	values := make([]colValue, len(indexes), len(indexes)+1)
	for i, index := range indexes {
		values[i] = c.values[index]
	}

	return &ByteColumn{
		Name:   c.Name,
		values: values,
	}
}
