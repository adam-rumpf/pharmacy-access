# Graphs

This subdirectory includes pure graph representations of the road maps in the `maps/` subdirectory. These files are created by the `travel_times.py` script in the `scrips/` subdirectory.

The node indices used throughout all files in this subdirectory correspond to those used in the `.graphml` files in the `maps/` subdirectory.

## Arc File Format

Each location's graph is defined using a TSV file. Each row defines an arc with the following fields:

* `tail`: Tail node ID.
* `head`: Head node ID.
* `weight`: Time to traverse the arc (in seconds).

## Population and Facility File Formats

In order to associate the population centers and vaccine providers with nodes in the graph file, a population file and a facility file are also generated for each location. Both consist of two fields:

* `id`: ID number of the population center or vaccination facility in its respective population or facility file.
* `node`: ID number of the graph node corresponding to this location.

In particular, each location is snapped to its nearest graph node (using the osmnx `nearest_nodes` function).
