package mbtiles

import (
	"bytes"
	"compress/gzip"
	"context"
	"encoding/json"
	"fmt"
	"os"
	"path/filepath"

	"crawshaw.io/sqlite"
	"crawshaw.io/sqlite/sqlitex"
	"github.com/brendan-ward/arrowtiler/mvt"
	"github.com/brendan-ward/arrowtiler/tiles"
)

var emptyContext context.Context

type MBtilesWriter struct {
	pool *sqlitex.Pool
}

const init_sql = `
CREATE TABLE metadata (name text, value text);
CREATE TABLE tiles (zoom_level integer, tile_column integer, tile_row integer, tile_data blob);
CREATE UNIQUE INDEX name on metadata (name);
CREATE UNIQUE INDEX tile_index on tiles (zoom_level, tile_column, tile_row);
`

func NewMBtilesWriter(path string, poolsize int) (*MBtilesWriter, error) {
	ext := filepath.Ext(path)
	if ext != ".mbtiles" {
		return nil, fmt.Errorf("path must end in .mbtiles")
	}

	// always overwrite
	if _, err := os.Stat(path); !os.IsNotExist(err) {
		os.Remove(path)
	}

	// check flags: this may not be safe for multiple goroutines (only one write  per connection though)
	pool, err := sqlitex.Open(path, sqlite.SQLITE_OPEN_CREATE|sqlite.SQLITE_OPEN_READWRITE|sqlite.SQLITE_OPEN_NOMUTEX|sqlite.SQLITE_OPEN_WAL, poolsize)
	if err != nil {
		return nil, err
	}

	db := &MBtilesWriter{
		pool: pool,
	}

	con, err := db.GetConnection()
	if err != nil {
		return nil, err
	}
	defer db.CloseConnection(con)

	// create tables
	err = sqlitex.ExecScript(con, init_sql)
	if err != nil {
		return nil, fmt.Errorf("could not initialize database: %q", err)
	}

	return db, nil
}

func (db *MBtilesWriter) Close() {
	if db.pool != nil {

		// make sure that anything pending is written
		con, err := db.GetConnection()
		if err != nil {
			panic(err)
		}
		err = sqlitex.Exec(con, `PRAGMA wal_checkpoint;`, nil)
		if err != nil {
			panic(err)
		}
		db.CloseConnection(con)

		db.pool.Close()
	}
}

// GetConnection gets a sqlite.Conn from an open connection pool.
// CloseConnection(con) must be called to release the connection.
func (db *MBtilesWriter) GetConnection() (*sqlite.Conn, error) {
	con := db.pool.Get(emptyContext)
	if con == nil {
		return nil, fmt.Errorf("connection could not be opened")
	}
	return con, nil
}

// CloseConnection closes an open sqlite.Conn and returns it to the pool.
func (db *MBtilesWriter) CloseConnection(con *sqlite.Conn) {
	if con != nil {
		db.pool.Put(con)
	}
}

func writeMetadataItem(con *sqlite.Conn, key string, value interface{}) error {
	return sqlitex.Exec(con, "INSERT INTO metadata (name,value) VALUES (?, ?)", nil, key, value)
}

func (db *MBtilesWriter) WriteMetadata(name string, description string, minZoom uint32, maxZoom uint32, bounds [4]float64, layersInfo map[string][]*mvt.LayerInfo) (err error) {
	if db == nil || db.pool == nil {
		return fmt.Errorf("cannot write to closed mbtiles database")
	}

	con, e := db.GetConnection()
	if e != nil {
		return e
	}
	defer db.CloseConnection(con)

	// create savepoint
	defer sqlitex.Save(con)(&err)

	if err = writeMetadataItem(con, "name", name); err != nil {
		return err
	}
	if err = writeMetadataItem(con, "description", description); err != nil {
		return err
	}
	if err = writeMetadataItem(con, "minzoom", minZoom); err != nil {
		return err
	}
	if err = writeMetadataItem(con, "maxzoom", maxZoom); err != nil {
		return err
	}
	if err = writeMetadataItem(con, "center", fmt.Sprintf("%.5f,%.5f,%v", (bounds[2]-bounds[0])/2.0, (bounds[3]-bounds[1])/2.0, minZoom)); err != nil {
		return err
	}
	if err = writeMetadataItem(con, "bounds", fmt.Sprintf("%.5f,%.5f,%.5f,%.5f", bounds[0], bounds[1], bounds[2], bounds[3])); err != nil {
		return err
	}
	if err = writeMetadataItem(con, "type", "overlay"); err != nil {
		return err
	}
	if err = writeMetadataItem(con, "format", "pbf"); err != nil {
		return err
	}
	if err = writeMetadataItem(con, "version", 2); err != nil {
		return err
	}

	layerInfoJSON, err := json.Marshal(layersInfo)
	if err != nil {
		return err
	}
	if err = writeMetadataItem(con, "json", layerInfoJSON); err != nil {
		return err
	}

	return nil
}

func (db *MBtilesWriter) WriteTile(tile *tiles.TileID, data []byte) error {
	con, err := db.GetConnection()
	if err != nil {
		return err
	}
	defer db.CloseConnection(con)

	return WriteTile(con, tile, data)
}

// Write the tile to the open connection
func WriteTile(con *sqlite.Conn, tile *tiles.TileID, data []byte) error {
	// GZIP tile data
	var b bytes.Buffer
	gz := gzip.NewWriter(&b)

	_, err := gz.Write(data)
	if err != nil {
		return err
	}
	if err = gz.Flush(); err != nil {
		return err
	}
	if err = gz.Close(); err != nil {
		return err
	}

	// flip tile Y to match mbtiles spec
	y := (1 << tile.Zoom) - 1 - tile.Y

	err = sqlitex.Exec(con, "INSERT INTO tiles (zoom_level, tile_column, tile_row, tile_data) VALUES (?, ?, ?, ?)", nil, tile.Zoom, tile.X, y, b.Bytes())
	if err != nil {
		return fmt.Errorf("could not write tile %v to mbtiles: %q", tile, err)
	}

	return nil
}
