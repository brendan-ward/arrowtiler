package cmd

import (
	"errors"
	"fmt"
	"os"
	"path"
	"path/filepath"
	"strings"
	"sync"

	"github.com/brendan-ward/arrowtiler/mbtiles"
	"github.com/brendan-ward/arrowtiler/mvt"
	"github.com/brendan-ward/arrowtiler/tiles"
	"github.com/spf13/cobra"

	"github.com/gosuri/uiprogress"
)

var minzoom uint16
var maxzoom uint16
var tilesetName string
var layerName string
var description string
var numWorkers int
var idColumm string
var extent uint16
var buffer uint16
var precision uint8
var simplification uint8

var createCmd = &cobra.Command{
	Use:   "create [IN.feather] [OUT.mbtiles]",
	Short: "Create a MVT tileset from a GeoArrow file",
	Args: func(cmd *cobra.Command, args []string) error {
		if len(args) < 2 {
			return errors.New("feather and mbtiles filenames are required")
		}
		if _, err := os.Stat(args[0]); errors.Is(err, os.ErrNotExist) {
			return fmt.Errorf("input file '%s' does not exist", args[0])
		}
		outDir, _ := path.Split(args[1])
		if outDir != "" {
			if _, err := os.Stat(outDir); errors.Is(err, os.ErrNotExist) {
				return fmt.Errorf("output directory '%s' does not exist", outDir)
			}
		}
		if path.Ext(args[1]) != ".mbtiles" {
			return errors.New("mbtiles filename must end in '.mbtiles'")
		}
		return nil
	},
	RunE: func(cmd *cobra.Command, args []string) error {
		// validate flags
		if numWorkers < 1 {
			numWorkers = 1
		}
		if maxzoom < minzoom {
			return errors.New("maxzoom must be no smaller than minzoom")
		}
		if precision < 1 {
			precision = 1
		}
		if simplification < 1 {
			simplification = 1
		}

		return create(args[0], args[1])
	},
	SilenceUsage: true,
}

func init() {
	createCmd.Flags().Uint16VarP(&minzoom, "minzoom", "Z", 0, "minimum zoom level")
	createCmd.Flags().Uint16VarP(&maxzoom, "maxzoom", "z", 0, "maximum zoom level")
	createCmd.Flags().StringVarP(&layerName, "layer", "l", "", "layer name")
	createCmd.Flags().StringVarP(&tilesetName, "name", "n", "", "tileset name")
	createCmd.Flags().StringVar(&description, "description", "", "tileset description")
	createCmd.Flags().StringVar(&idColumm, "id", "", "integer column name to use as feature ID")
	createCmd.Flags().IntVarP(&numWorkers, "workers", "w", 4, "number of workers to create tiles")
	createCmd.Flags().Uint16Var(&extent, "extent", 4096, "extent of tile")
	createCmd.Flags().Uint16Var(&buffer, "buffer", 256, "number of pixels to buffer outside extent before clipping lines / polygons")
	createCmd.Flags().Uint8VarP(&precision, "precision", "p", 1, "precision used to snap coordinates, in pixels; points are deduplicated to this value")
	createCmd.Flags().Uint8VarP(&simplification, "simplify", "s", 1, "simplification factor used to simplify lines / polygons, in pixels")
}

func produce(minZoom uint16, maxZoom uint16, bounds [4]float64, queue chan<- *tiles.TileID) {
	defer close(queue)

	fmt.Println("Creating tiles")

	uiprogress.Start()

	for zoom := minZoom; zoom <= maxZoom; zoom++ {
		z := zoom
		minTile, maxTile := tiles.TileRange(zoom, bounds)
		count := ((maxTile.X - minTile.X) * (maxTile.Y - minTile.Y)) + 1
		bar := uiprogress.AddBar(int(count)).AppendCompleted().PrependElapsed()
		bar.PrependFunc(func(b *uiprogress.Bar) string {
			return fmt.Sprintf("zoom %2v (%8v/%8v)", z, b.Current(), count)
		})

		for x := minTile.X; x <= maxTile.X; x++ {
			for y := minTile.Y; y <= maxTile.Y; y++ {
				queue <- &tiles.TileID{Zoom: zoom, X: x, Y: y}
				bar.Incr()
			}
		}
	}
	uiprogress.Stop()
}

func create(infilename string, outfilename string) error {
	// set defaults
	if layerName == "" {
		layerName = strings.TrimSuffix(path.Base(infilename), filepath.Ext(infilename))
	}
	if tilesetName == "" {
		tilesetName = layerName
	}

	// coordinates projected to Mercator on read
	fmt.Printf("Reading features from %v\n", infilename)
	features, err := mvt.ReadFeather(infilename, idColumm)
	if err != nil {
		return err
	}
	defer features.Geometry().Release()

	db, err := mbtiles.NewMBtilesWriter(outfilename, numWorkers)
	if err != nil {
		panic(err)
	}
	defer db.Close()

	fields, err := features.GetLayerInfo(layerName, description, minzoom, maxzoom)
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

	db.WriteMetadata(tilesetName, description, minzoom, maxzoom, geoBounds, layersInfo)

	config := tiles.NewEncodingConfig(extent, buffer, precision, simplification)

	queue := make(chan *tiles.TileID)
	var wg sync.WaitGroup

	go produce(minzoom, maxzoom, mercatorBounds, queue)

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
				tile, err := features.EncodeToLayer(layerName, tileID, config)
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

	return nil
}
