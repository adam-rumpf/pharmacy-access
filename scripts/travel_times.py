"""COVID-19 Vaccine accessibility project travel time scripts.

The code below includes various scripts for generating required data from the
preprocessed data files, notably including origin/destination travel time
matrices.
"""

import os.path
import time

import osmnx as ox
import tqdm

#==============================================================================
# Global Constants
#==============================================================================

# Distance center file header (including column labels)
DIST_HEADER = "pid\tfid\tpftime\tfptime\t\n"

# Preprocessed file column numbers
FAC_LAT = 2 # facility latitude
FAC_LON = 3 # facility longitude
FAC_CAP = 4 # facility capacity
POP_LAT = 2 # population latitude
POP_LON = 3 # population longitude
POP_POP = 4 # population total population value

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
    t = time.time()
    G = ox.graph.graph_from_bbox(north, south, east, west, network_type="drive")
    print(f"Map downloaded after {time.time()-t} seconds.")
    
    # Impute missing edge speeds then calculate edge travel times
    print("Imputing edge speeds.")
    t = time.time()
    G = ox.add_edge_speeds(G) # kph
    G = ox.add_edge_travel_times(G)
    print(f"Edge speeds added after {time.time()-t} seconds.")
    
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

#------------------------------------------------------------------------------

def _read_popfile(popfile):
    """Reads a population file and returns dictionaries of data.
    
    Positional arguments:
        popfile (str) -- Preprocessed population file path, which should
            include the coordinates and population of each population
            cetner.
    
    Returns:
        pop (dict(int)) -- Dictionary of population counts.
        coord (dict((float,float))) -- Dictionary of population center
            coordinates, as (latitude,longitude) tuples.
    
    All dictionaries are indexed by the population center IDs contained in the
    first column of the population file.
    """
    
    # Initialize dictionaries
    pop = dict() # dictionary of populations by popfile index
    coord = dict() # dictionary of latitutde/longitude pairs by popfile index
    
    # Read file
    with open(popfile, 'r') as f:
        
        for line in f:
            
            # Skip the comment line
            if line[0].isdigit() == False:
                continue
            
            # Get numbers from line
            s = line.strip().split('\t')
            pop[int(s[0])] = int(s[POP_POP])
            coord[int(s[0])] = (float(s[POP_LAT]), float(s[POP_LON]))
    
    return (pop, coord)

#------------------------------------------------------------------------------

def _read_facfile(facfile):
    """Reads a facility file and returns dictionaries of data.
    
    Positional arguments:
        facfile (str) -- Preprocessed facility file path, which should include
            the coordinates and capacity of each vaccination facility.
    
    Returns:
        cap (dict(int)) -- Dictionary of facility capacities.
        coord (dict((float,float))) -- Dictionary of facility coordinates,
            as (latitude,longitude) tuples.
    
    All dictionaries are indexed by the facility IDs contained in the first
    column of the facility file.
    """
    
    # Initialize dictionaries
    cap = dict() # dictionary of facility capacities by facfile index
    coord = dict() # dictionary of latitude/longitude pairs by facfile index
    
    # Read file
    with open(facfile, 'r') as f:
        
        for line in f:
            
            # Skip the comment line
            if line[0].isdigit() == False:
                continue
            
            # Get numbers from line
            s = line.strip().split('\t')
            cap[int(s[0])] = float(s[FAC_CAP])
            coord[int(s[0])] = (float(s[FAC_LAT]), float(s[FAC_LON]))
    
    return (cap, coord)

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
# Batch Processing Scripts
#==============================================================================

def generate_distance_file(popfile, facfile, mapfile, distfile):
    """Generates a travel time file for a batch of population/facility files.
    
    Positional arguments:
        popfile (str) -- Preprocessed population file path, which should
            include the coordinates and population of each population center.
        facfile (str) -- Preprocessed facility file path, which should include
            the coordinates and capacity of each vaccination facility.
        mapfile (str) -- File path to a pre-downloaded .graphml map file.
        distfile (str) -- File path for output distance file.
    
    The distance file includes a table of all population/facility pairs and
    the (directional) pairwise distances between each pair.
    """
    
    # Get population file coordinates
    pdic = _read_popfile(popfile)[1] # population coordinate dictionary
    
    # Get facility file coordinates
    fdic = _read_facfile(facfile)[1] # facility coordinate dictionary
    
    # Initialize dictionaries of distances
    pfdist = dict() # population to facility distance, indexed by (p,f) tuple
    fpdist = dict() # facility to population distance, indexed by (p,f) tuple
    
    # Find all population-to-facility distances
    print("Computing all population-to-facility distances.")
    t = time.time()
    for pid in tqdm.tqdm(pdic):
        
        # Find all distances from the current population center
        ###dlist = [fdic[fid] for fid in fdic] # list of all facility coords
        first_ids = list(fdic.keys())[0:2]###
        dlist = [fdic[i] for i in first_ids]###
        d = travel_time_destination_list(pdic[pid], dlist, mapfile)
        
        # Save distances in dictionary
        fids = list(fdic.keys()) # list of facility IDs
        ###for j in range(len(fdic)):
        for j in range(2):###
            pfdist[(pid, fids[j])] = d[j]
        
        break###
    
    print(f"All pairs processed after {time.time()-t} seconds.")
    
    # Find all facility-to-population distances
    print("Computing all facility-to-population distances.")
    t = time.time()
    for fid in tqdm.tqdm(fdic):
        
        # Find all distances from the current facility
        ###dlist = [pdic[pid] for pid in pdic] # list of all population coords
        first_ids = list(pdic.keys())[0:2]###
        dlist = [pdic[i] for i in first_ids]###
        d = travel_time_destination_list(fdic[fid], dlist, mapfile)
        
        # Save distances in dictionary
        pids = list(pdic.keys()) # list of population IDs
        ###for j in range(len(pdic)):
        for j in range(2):###
            fpdist[(pid, fids[j])] = d[j]
        
        break###
    
    # Write dictionary contents to distance file
    with open(distfile, 'w') as f:
        f.write(DIST_HEADER)
        for pair in pfdist:
            f.write(str(pair[0]) + "\t" + str(pair[1]) + "\t" +
                str(pfdist(pair)) + "\t" + str(fpdist(pair)) + "\t\n")

#==============================================================================
# Execution
#==============================================================================

# Comment or uncomment the function calls below to process each location.
#download_map_box(os.path.join("..", "maps", "santa_clara", "santa_clara_driving_small.graphml"), 37.4323, 37.1986, -121.7395, -122.0924)
#print(travel_time((37.4, -122.0), (37.2, -121.7), os.path.join("..", "maps", "santa_clara", "santa_clara_driving_small.graphml")))
#print(travel_time_destination_list((37.4, -122.0), [(37.2, -121.7), (37.4, -121.7)], os.path.join("..", "maps", "santa_clara", "santa_clara_driving_small.graphml")))
#download_map_box(os.path.join("..", "maps", "santa_clara", "santa_clara_driving.graphml"), SANTA_CLARA_N, SANTA_CLARA_S, SANTA_CLARA_E, SANTA_CLARA_W)
#print(travel_time((37.4, -122.0), (37.2, -121.7), os.path.join("..", "maps", "santa_clara", "santa_clara_driving.graphml")))
generate_distance_file(os.path.join("..", "processed", "santa_clara", "santa_clara_pop.tsv"), os.path.join("..", "processed", "santa_clara", "santa_clara_fac.tsv"), os.path.join("..", "maps", "santa_clara", "santa_clara_driving.graphml"), os.path.join("..", "processed", "santa_clara", "santa_clara_dist.tsv"))
