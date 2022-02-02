package cmd

import (
	"github.com/spf13/cobra"
)

var rootCmd = &cobra.Command{
	Use:   "arrowtiler",
	Short: "A Mapbox Vector Tile generator for GeoArrow files",
}

// Execute executes the root command.
func Execute() error {
	return rootCmd.Execute()
}

func init() {

	rootCmd.AddCommand(createCmd)
}
