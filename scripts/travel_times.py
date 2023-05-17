"""COVID-19 Vaccine accessibility project travel time scripts.

The code below includes various scripts for generating required data from the
preprocessed data files, notably including origin/destination travel time
matrices.
"""

import os.path

import osmnx as ox
import tqdm

#==============================================================================
# Global Constants
#==============================================================================

# Distance center file header (including column labels)
DIST_HEADER = "pid\tfid\tpftime\tfptime\t\n"

# Location boundaries (large enough to include neighboring counties)
SANTA_CLARA_N = 38.428
SANTA_CLARA_S = 35.791
SANTA_CLARA_E = -120.045
SANTA_CLARA_W = -123.162

#==============================================================================
# Common Functions
#==============================================================================

def download_map_box(outfile, north, south, east, west):
    """Downloads and saves OpenStreetMap data for a rectangular region.
    
    Positional arguments:
        outfile (str) -- Output file path.
        north (float) -- Northern latitude bound.
        south (float) -- Southern latitude bound.
        east (float) -- Eastern longitude bound.
        west (float) -- Western longitude bound.
    
    This function downloads and saves a map of a given bounding box as a
    .graphml file. The map is downloaded using OpenStreetMap driving network
    data.
    """
    
    # Verify that the output file has the .graphml extension
    if os.path.splitext(outfile)[1] != ".graphml":
        raise ValueError("output file must have the '.graphml' extension")
    
    # Download the map
    print("Downloading map data.")
    G = ox.graph.graph_from_bbox(north, south, east, west, network_type="drive")
    
    # Impute missing edge speeds then calculate edge travel times
    print("Imputing edge speeds.")
    G = ox.add_edge_speeds(G) # kph
    G = ox.add_edge_travel_times(G)
    print("Edge speeds added.")
    
    # Save the map
    ox.save_graphml(G, outfile)
    print("Map successfully saved as " + str(outfile))

#------------------------------------------------------------------------------

def travel_time(origin, destination, mapfile=None):
    """Computes the travel time between a given pair of coordinates.
    
    Positional arguments:
        origin (tuple(float)) -- Origin latitude/longitude tuple.
        destination (tuple(float)) -- Destination latitude/longitude tuple.
    
    Keyword arguments:
        mapfile (str) -- File path to a pre-downloaded .graphml map file.
            Defaults to None, in which case new map data is downloaded and
            immediately discarded. The downloaded data consists of the bounding
            box defined by the origin and destination coordinates, plus a 15%
            buffer.
    
    Returns:
        (float) -- Travel time (minutes).
    """
    
    # Open specified map file
    if mapfile != None:
        
        # Verify file extension
        if os.path.splitext(mapfile)[1] != ".graphml":
            raise ValueError("output file must have the '.graphml' extension")
        
        # Open file
        G = ox.load_graphml(mapfile)
    
    else:
        
        # If no map file exists, find the origin/destination bounding box
        north = max(origin[1], destination[1])
        south = min(origin[1], destination[1])
        east = max(origin[0], destination[0])
        west = min(origin[0], destination[0])
        dim = max(north-south, east-west)
        north += 0.15*dim
        south -= 0.15*dim
        east += 0.15*dim
        west -= 0.15*dim
        
        # Generate graph from bounding box
        G = ox.graph.graph_from_bbox(north, south, east, west,
                                     network_type="drive")
        
        # Impute missing edge speeds then calculate edge travel times
        G = ox.add_edge_speeds(G) # kph
        G = ox.add_edge_travel_times(G)
    
    # Get nodes nearest to origin and destination coordinates
    onode = ox.distance.nearest_nodes(G, origin[1], origin[0])
    dnode = ox.distance.nearest_nodes(G, destination[1], destination[0])
    
    # Find the shortest path by travel time
    route = ox.distance.shortest_path(G, onode, dnode, weight="travel_time")
    
    # Get all edge times
    edge_times = ox.utils_graph.get_route_edge_attributes(G, route,
                                                          "travel_time")
    
    # Return the total edge times, converted to minutes
    return sum(edge_times)/60.0

#==============================================================================
# Location-Specific Processing Scripts
#==============================================================================

###
### Use tqdm.tqdm for OD pair progress bars

#==============================================================================
# Execution
#==============================================================================

# Comment or uncomment the function calls below to process each location.
#download_map_box(os.path.join("..", "maps", "santa_clara", "santa_clara_driving.graphml"), 37.4323, 37.1986, -121.7395, -122.0924)
print(travel_time((37.4, -122.0), (37.2, -121.7), os.path.join("..", "maps", "santa_clara", "santa_clara_driving.graphml")))
