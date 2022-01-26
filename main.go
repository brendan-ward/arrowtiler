package main

import (
	"fmt"
	"sync"

	"github.com/brendan-ward/arrowtiler/mbtiles"
	"github.com/brendan-ward/arrowtiler/mvt"
	"github.com/brendan-ward/arrowtiler/tiles"
	"github.com/schollz/progressbar/v3"
)

func produce(minZoom uint32, maxZoom uint32, bounds [4]float64, queue chan<- *tiles.TileID) {
	defer close(queue)

	for zoom := minZoom; zoom <= maxZoom; zoom++ {
		minTile, maxTile := tiles.TileRange(zoom, bounds)
		count := (maxTile.X - minTile.X) * (maxTile.Y - minTile.Y)

		// simplistic bar based on rate that we are producing tiles rather than consuming them
		// TODO: put this on a channel with callbacks from consumers instead
		bar := progressbar.NewOptions(int(count), progressbar.OptionSetWidth(25), progressbar.OptionSetDescription(fmt.Sprintf("zoom: %v", zoom)))
		for x := minTile.X; x <= maxTile.X; x++ {
			for y := minTile.Y; y <= maxTile.Y; y++ {
				queue <- &tiles.TileID{Zoom: zoom, X: x, Y: y}
				bar.Add(1)
			}
		}
		bar.Clear()
	}
}

func main() {

	layerName := "test"
	description := "description"
	var minZoom uint32 = 0
	var maxZoom uint32 = 5
	numWorkers := 4
	filename := "/tmp/test.feather"

	// coordinates projected to Mercator on read
	features, err := mvt.ReadFeather(filename, "geometry", "")
	if err != nil {
		panic(err)
	}
	defer features.Geometry().Release()

	db, err := mbtiles.NewMBtilesWriter("/tmp/tiles/out.mbtiles", numWorkers)
	if err != nil {
		panic(err)
	}
	defer db.Close()

	fields, err := features.GetLayerInfo(layerName, description, minZoom, maxZoom)
	if err != nil {
		panic(err)
	}

	layersInfo := make(map[string][]*mvt.LayerInfo)
	layersInfo["vector_layers"] = []*mvt.LayerInfo{
		fields,
	}

	mercatorBounds, err := features.Geometry().TotalBounds()
	if err != nil {
		panic(err)
	}
	// project Mercator bound to geo bounds
	geoBounds := mvt.MercatorBoundsToGeoBounds(mercatorBounds)

	db.WriteMetadata("name", "description", minZoom, maxZoom, geoBounds, layersInfo)

	//	fmt.Printf("read and encoded %d records in %v\n", t.Size(), time.Since(start))

	queue := make(chan *tiles.TileID, 1000)
	var wg sync.WaitGroup

	go produce(minZoom, maxZoom, geoBounds, queue)

	for i := 0; i < numWorkers; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()

			con, err := db.GetConnection()
			if err != nil {
				panic(err)
			}
			defer db.CloseConnection(con)

			for tileID := range queue {
				tile, err := features.EncodeToLayer(layerName, tileID)
				if err != nil {
					panic(err)
				}

				if tile != nil {
					mbtiles.WriteTile(con, tileID, tile)
				}
			}
		}()
	}

	wg.Wait()

	fmt.Println("Done!")
}
