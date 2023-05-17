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

def _temporary_map(orig, dest, buffer=0.15):
    """Download OpenStreetMap data for a rectangular region without saving.
    
    Positional arguments:
        orig (tuple(float)) -- Origin latitude/longitude tuple.
        dest (tuple(float)) -- Destination latitude/longitude tuple.
    
    Keyword arguments:
        buffer (float) -- Fractional buffer to apply the boundaries of the
            bounding box. Defaults to 0.15, which indicates that the bounding
            box will be 15% larger than the smallest box containing the given
            origin and destination.
    
    Returns:
        (ox.graph) -- The graph corresponding to the selected bounding box,
            plus a given buffer.
    """
    
    # Find the origin/destination bounding box
    north = max(orig[1], dest[1])
    south = min(orig[1], dest[1])
    east = max(orig[0], dest[0])
    west = min(orig[0], dest[0])
    dim = max(north - south, east - west)
    north += buffer*dim
    south -= buffer*dim
    east += buffer*dim
    west -= buffer*dim
    
    # Generate graph from bounding box
    G = ox.graph.graph_from_bbox(north, south, east, west, network_type="drive")
    
    # Impute missing edge speeds then calculate edge travel times
    G = ox.add_edge_speeds(G) # kph
    G = ox.add_edge_travel_times(G)
    
    # Return the graph
    return G

#==============================================================================
# Travel Time Computation Scripts
#==============================================================================

def travel_time(orig, dest, mapfile=None):
    """Computes the travel time between a given pair of coordinates.
    
    Positional arguments:
        orig (tuple(float)) -- Origin latitude/longitude tuple.
        dest (tuple(float)) -- Destination latitude/longitude tuple.
    
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
        
        # Generate graph from file
        G = ox.load_graphml(mapfile)
    
    else:
        
        # If no map file exists, generate the graph from scratch
        G = _temporary_map(orig, dest)
    
    # Get nodes nearest to origin and destination coordinates
    onode = ox.distance.nearest_nodes(G, orig[1], orig[0])
    dnode = ox.distance.nearest_nodes(G, dest[1], dest[0])
    
    # Find the shortest path by travel time
    route = ox.distance.shortest_path(G, onode, dnode, weight="travel_time")
    
    # Get all edge times
    edge_times = ox.utils_graph.get_route_edge_attributes(G, route,
                                                          "travel_time")
    
    # Return the total edge times, converted to minutes
    return sum(edge_times)/60.0

#------------------------------------------------------------------------------

def travel_time_destination_list(orig, dlist, mapfile=None):
    """Computes the travel times from one origin to a list of destinations.
    
    Positional arguments:
        orig (tuple(float)) -- Origin latitude/longitude tuple.
        dlist (list(tuple(float))) -- List of destination latitude/longitude
            tuples.
    
    Keyword arguments:
        mapfile (str) -- File path to a pre-downloaded .graphml map file.
            Defaults to None, in which case new map data is downloaded and
            immediately discarded. The downloaded data consists of the bounding
            box defined by the origin and destination coordinates, plus a 15%
            buffer.
    
    Returns:
        (list(float)) -- List of travel times from the single origin to all
            destinations in the list.
    """
    
    # Open specified map file
    if mapfile != None:
        
        # Verify file extension
        if os.path.splitext(mapfile)[1] != ".graphml":
            raise ValueError("output file must have the '.graphml' extension")
        
        # Generate graph from file
        G = ox.load_graphml(mapfile)
    
    else:
        
        # If no map file exists, generate the graph from scratch
        G = _temporary_map(orig, dest)
    
    # Get size of destination list
    n = len(dlist)
    
    # Generate lists of origin and destination nodes (origin list is constant)
    onodes = [ox.distance.nearest_nodes(G, orig[1], orig[0]) for i in range(n)]
    dnodes = [ox.distance.nearest_nodes(G, dlist[i][1], dlist[i][0])
              for i in range(n)]
    
    # Find all shortest path routes
    routes = ox.distance.shortest_path(G, onodes, dnodes, weight="travel_time",
                                       cpus=None)
    
    # Return a list of all pairwise travel times
    return [sum(ox.utils_graph.get_route_edge_attributes(G, routes[i],
            "travel_time"))/60.0 for i in range(n)]

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
#print(travel_time((37.4, -122.0), (37.2, -121.7), os.path.join("..", "maps", "santa_clara", "santa_clara_driving.graphml")))
print(travel_time_destination_list((37.4, -122.0), [(37.2, -121.7), (37.4, -121.7)], os.path.join("..", "maps", "santa_clara", "santa_clara_driving.graphml")))
