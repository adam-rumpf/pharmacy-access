"""COVID-19 Vaccine accessibility project travel time scripts.

The code below includes various scripts for generating required data from the
preprocessed data files, notably including origin/destination travel time
matrices.
"""

import heapq
import os.path
import time
import xml.etree.ElementTree as ET

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

# Counties neighboring Santa Clara
SANTA_CLARA_NEIGHBORS = ["Alameda", "Merced", "Monterey", "San Benito",
                         "San Francisco", "San Joaquin", "San Mateo",
                         "Santa Cruz", "Stanislaus"]

#==============================================================================
# Data Gathering and Conversion
#==============================================================================

def download_map_places(outfile, places):
    """Downloads and saves OpenStreetMap data for a list of counties.
    
    Positional arguments:
        outfile (str) -- Output file path.
        places (list(dict)) -- List of dictionaries containing county queries.
            Each such query should be a dictionary of the form:
            {"county": "COUNTYNAME", "state": "STATENAME"}
    
    This function downloads and saves a map as.graphml file. The map is
    downloaded using OpenStreetMap driving network data.
    """
    
    # Verify that the output file has the .graphml extension
    if os.path.splitext(outfile)[1] != ".graphml":
        raise ValueError("output file must have the '.graphml' extension")
    
    # Download the map
    print("Downloading map data.")
    t = time.time()
    G = ox.graph.graph_from_place(places, network_type="drive")
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

#------------------------------------------------------------------------------

def graphml_to_tsv(mapfile, arcfile, weight="travel_time"):
    """Converts a .graphml map into a simplified pure graph file.
    
    Positional arguments:
        mapfile (str) -- Path to input .graphml file.
        arcfile (str) -- Path to output arc data file.
    
    Keyword arguments:
        weight (str) -- Edge attribute to use for defining arc weights.
            Defaults to "travel_time".
    
    The output graph file simply defines arc-level data since there we have no
    need for node-level data. The arc file includes the following tab-separated
    columns, in order:
        tail -- Node index of arc's tail (origin).
        head -- Node index of arc's head (destination).
        weight -- Numerical weight to use in computing distances.
    
    The node indices in the output file correspond to those used in the
    original .graphml file.
    """
    
    # Define default .graphml namespace
    namespace = "{http://graphml.graphdrawing.org/xmlns}"
    
    # Read .graphml file as an XML element tree
    tree = ET.parse(mapfile)
    root = tree.getroot()
    
    # Find the ID number of the given weight attribute
    wid = None
    for child in root.findall(namespace + "key"):
        if (child.get("for") == "edge") and (child.get("attr.name") == weight):
            wid = child.get("id")
            break
    
    # Report an error if the weight field was not found
    if wid == None:
        ValueError("weight field '" + weight + "' not found in graphml file")
    
    # Find the "graph" tag
    for child in root:
        if child.tag != namespace + "graph":
            continue
        
        # Initialize the edge list
        edges = [[None, None, None] for i in
                  range(len(child.findall(namespace + "edge")))]
        
        # Go through each edge one-at-a-time
        i = 0
        for e in tqdm.tqdm(child.findall(namespace + "edge")):
            
            edges[i][0] = e.attrib["source"] # tail
            edges[i][1] = e.attrib["target"] # head
            
            for f in e:
                if f.attrib["key"] == wid:
                    edges[i][2] = f.text
                    break
            
            i += 1
    
    # Write output file
    with open(arcfile, 'w') as f:
        
        # Write header
        f.write("tail\thead\tweight\n")
        
        # Write array contents to file
        for i in range(len(edges)):
            f.write(str(edges[i][0]) + "\t" + str(edges[i][1]) + "\t" +
                    str(edges[i][2]) + "\n")

#------------------------------------------------------------------------------

def map_node_locations(mapfile, popfile, pnodefile, facfile, fnodefile):
    """Finds the nodes that correspond to population and facility locations.
    
    Positional arguments:
        mapfile (str) -- Path to input .graphml file.
        popfile (str) -- Path to preprocessed population file.
        pnodefile (str) -- Path to population node output file.
        facfile (str) -- Path to preprocessed facility file.
        fnodefile (str) -- Path to facility node output file.
    
    Each population center is snapped to its nearest node in the .graphml file
    based on the geographic coordinates in the popfile. The output pnodefile
    is a tab-separated file containing the following columns:
        id -- Population center index from the population file.
        node -- Corresponding node index from the .graphml file.
    
    The facfile and fnodefile arguments are analogous, but are for mapping
    facilities to their nearest nodes.
    """
    
    # Open graph file
    G = ox.load_graphml(mapfile)
    
    # Get population center coordinates
    (_, pcoord) = _read_popfile(popfile)
    pid = list(pcoord) # get population center IDs
    
    # Find nodes nearest each population center
    print("Mapping population nodes.")
    pnodes = ox.distance.nearest_nodes(G, [pcoord[i][1] for i in pid],
                                          [pcoord[i][0] for i in pid])
    
    # Write output file
    with open(pnodefile, 'w') as f:
        f.write("id\tnode\n")
        for i in range(len(pid)):
            f.write(str(pid[i]) + "\t" + str(pnodes[i]) + "\n")
    
    del pid
    del pcoord
    
    # Get facility coordinates
    (_, fcoord) = _read_facfile(facfile)
    fid = list(fcoord) # get facility IDs
    
    # Find nodes nearest each facility
    print("Mapping facilities.")
    fnodes = ox.distance.nearest_nodes(G, [fcoord[i][1] for i in fid],
                                          [fcoord[i][0] for i in fid])
    
    # Write output file
    with open(fnodefile, 'w') as f:
        f.write("id\tnode\n")
        for i in range(len(fid)):
            f.write(str(fid[i]) + "\t" + str(fnodes[i]) + "\n")
    
    del fid
    del fcoord

#------------------------------------------------------------------------------

def _generate_adjacency_list(arcfile, pnodefile, fnodefile, factor=10):
    """Generates an adjacency list representation of a graph and its O/D nodes.
    
    Positional arguments:
        arcfile (str) -- Path to the network's arc definition file, which
            should contain the tail node, head node, and travel time for each
            arc.
        pnodefile (str) -- Path to the network's population node definition
            file, which should contain a list of origin nodes.
        fnodefile (str) -- Path to the network's facility node definition file,
            which should contain a list of destination nodes.
    
    Keyword arguments:
        factor (float) -- Factor by which to multiply all arc weights before
            casting them as integers. Defaults to 10.
    
    Returns:
        adj (tuple(tuple(int))) -- Tuple of tuples of out-neighbors. adj[i][j]
            indicates the jth out-neighbor of the ith node.
        weight (tuple(tuple(int))) -- Tuple of tuples of arc weights.
            weight[i][j] indicates the weight of the arc from node i to the
            jth out-neighbor of i (i.e. the weight of the arc from i to
            adj[i][j]).
        size (int) -- Number of nodes.
        onodes (list(int)) -- List of origin node indices.
        dnodes (list(int)) -- List of destination node indices.
    
    This function generates an adjacency list representation of a graph defined
    in an arc data file. It also performs a few preprocessing tasks to
    simplify the shortest path generation process, including relabeling the
    nodes and altering the arc weights.
    
    More specifically, the node indices in the arc file are arbitrary positive
    integers from the original .graphml file. This script begins by relabeling
    the indices to become a list of zero-indexed consecutive integers (0, 1, 2,
    ...) so that the node indices can be used as list positions. The relabeling
    simply associates each node ID with its position in the sorted node ID
    list. The same relabeling is performed on the origin and destination node
    lists, which is why the two relabeled lists are returned.
    
    The arc weights are expected to be float values from the original .graphml
    file. This script multiplies all weights by a given factor and then casts
    them as integers in order to speed up the shortest path computations. Any
    resulting path lengths should be divided by this same multiplicative
    factor.
    """
    
    # Generate a map of node IDs to their positions in the sorted ID list
    
    # Find all unique node labels
    labelset = set() # set of all unique node labels used
    with open(arcfile, 'r') as f:
        # Read all unique labels from arc file
        for line in f:
            # Skip comment line
            if line[0].isdigit() == False:
                continue
            s = line.strip().split()
            labelset.add(int(s[0]))
            labelset.add(int(s[1]))
    
    # Sort the node label list and use to generate a position map
    size = len(labelset) # number of nodes
    labels = list(labelset)
    del labelset
    labels.sort()
    position = {labels[i]: i for i in range(len(labels))} # relabeling map
    del labels
    
    # Relabel origin nodes
    onodes = [] # origin node list
    with open(pnodefile, 'r') as f:
        for line in f:
            # Skip comment line
            if line[0].isdigit() == False:
                continue
            s = line.strip().split()
            onodes.append(position[int(s[1])])
    
    # Relabel destination nodes
    dnodes = [] # destination node list
    with open(fnodefile, 'r') as f:
        for line in f:
            # Skip comment line
            if line[0].isdigit() == False:
                continue
            s = line.strip().split()
            dnodes.append(position[int(s[1])])
    
    
    
    
    
    
    
    ###
    return (None, None, size, onodes, dnodes)

#==============================================================================
# Shortest Path Computation
#==============================================================================

def pairwise_distances(arcfile, pnodefile, fnodefile, factor=10):
    """Computes all pairwise distances from an origin set to a destination set.
    
    Positional arguments:
        arcfile (str) -- Path to the network's arc definition file, which
            should contain the tail node, head node, and travel time for each
            arc.
        pnodefile (str) -- Path to the network's population node definition
            file, which should contain a list of origin nodes.
        fnodefile (str) -- Path to the network's facility node definition file,
            which should contain a list of destination nodes.
    
    Keyword arguments:
        factor (float) -- Factor by which to multiply all arc weights before
            casting them as integers. Defaults to 10.
    
    Returns:
        ### Decide on output format. Probably a list of lists.
    """
    
    # Generate adjacency list representation
    (adj, weight, size, onodes, dnodes) = _generate_adjacency_list(arcfile,
                                           pnodefile, fnodefile, factor=factor)
    
    ### Eventually edit this to write distances to a log file as it goes.
    
    pass###

#==============================================================================
# Batch Processing Scripts
#==============================================================================

def generate_distance_file(popfile, facfile, mapfile, distfile,
    symmetric=False):
    """Generates a travel time file for a batch of population/facility files.
    
    Positional arguments:
        popfile (str) -- Preprocessed population file path, which should
            include the coordinates and population of each population center.
        facfile (str) -- Preprocessed facility file path, which should include
            the coordinates and capacity of each vaccination facility.
        mapfile (str) -- File path to a pre-downloaded .graphml map file.
        distfile (str) -- File path for output distance file.
    
    Keyword arguments:
        symmetric (bool) -- Whether to treat pairwise travel times as
            symmetric. Defaults to False, in which case travel times are
            computed independently in both directions. If True, the population-
            to-facility travel time is used for both directions. This saves
            computation time at the cost of some realism.
    
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
        dlist = [fdic[fid] for fid in fdic] # list of all facility coords
        d = travel_time_destination_list(pdic[pid], dlist, mapfile)
        
        # Save distances in dictionary
        fids = list(fdic.keys()) # list of facility IDs
        for j in range(len(fdic)):
            pfdist[(pid, fids[j])] = d[j]
    
    print(f"All pairs processed after {time.time()-t} seconds.")
    
    ###
    print(pfdist)
    print(pfdist.keys())
    
    # Find all facility-to-population distances
    print("Computing all facility-to-population distances.")
    t = time.time()
    for fid in tqdm.tqdm(fdic):
        
        # Copy distances if symmetric
        if symmetric == True:
            pids = list(pdic.keys()) # list of population IDs
            for j in range(len(pdic)):
                fpdist[(pids[j], fid)] = pfdist[(pids[j], fid)]
            continue
        
        # Find all distances from the current facility
        dlist = [pdic[pid] for pid in pdic] # list of all population coords
        d = travel_time_destination_list(fdic[fid], dlist, mapfile)
        
        # Save distances in dictionary
        pids = list(pdic.keys()) # list of population IDs
        for j in range(len(pdic)):
            fpdist[(pids[j], fid)] = d[j]
    
    ###
    print(fpdist)
    print(fpdist.keys())
    
    # Write dictionary contents to distance file
    with open(distfile, 'w') as f:
        f.write(DIST_HEADER)
        for pair in pfdist:
            f.write(str(pair[0]) + "\t" + str(pair[1]) + "\t" +
                str(pfdist[pair]) + "\t" + str(fpdist[pair]) + "\t\n")

#==============================================================================
# Execution
#==============================================================================

# Comment or uncomment the function calls below to process each location.

# Generate Santa Clara place list
santa_clara_mapfile = os.path.join("..", "maps", "santa_clara", "santa_clara_map.graphml")
santa_clara_arcfile = os.path.join("..", "graphs", "santa_clara", "santa_clara_arcs.tsv")
santa_clara_popfile = os.path.join("..", "processed", "santa_clara", "santa_clara_pop.tsv")
santa_clara_facfile = os.path.join("..", "processed", "santa_clara", "santa_clara_fac.tsv")
santa_clara_pnodefile = os.path.join("..", "graphs", "santa_clara", "santa_clara_popnodes.tsv")
santa_clara_fnodefile = os.path.join("..", "graphs", "santa_clara", "santa_clara_facnodes.tsv")
#santa_clara_places = [{"county": c, "state": "California"} for c in SANTA_CLARA_NEIGHBORS] + [{"county": "Santa Clara", "state": "California"}]

# Download: 501.8285789489746 seconds
# Edge speeds: 34.205276012420654 seconds
#download_map_places(, santa_clara_places)

#graphml_to_tsv(santa_clara_mapfile, santa_clara_arcfile, weight="travel_time")
#map_node_locations(santa_clara_mapfile, santa_clara_popfile, santa_clara_pnodefile, santa_clara_facfile, santa_clara_fnodefile)
pairwise_distances(santa_clara_arcfile, santa_clara_pnodefile, santa_clara_fnodefile, factor=10)
