"""Pharmacy accessibility project travel time scripts.

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

# Counties neighboring Santa Clara, CA
SANTA_CLARA_NEIGHBORS = ["Alameda", "Merced", "Monterey", "San Benito",
                         "San Francisco", "San Joaquin", "San Mateo",
                         "Santa Cruz", "Stanislaus"]

# Counties neighboring Polk, FL
POLK_NEIGHBORS = ["Hardee", "Highlands", "Hillsborough", "Lake", "Manatee",
                  "Okeechobee", "Orange", "Osceola", "Pasco", "Sumter"]

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

def graphml_to_tsv(mapfile, arcfile, weight="travel_time", directed=True):
    """Converts a .graphml map into a simplified pure graph file.
    
    Positional arguments:
        mapfile (str) -- Path to input .graphml file.
        arcfile (str) -- Path to output arc data file.
    
    Keyword arguments:
        weight (str) -- Edge attribute to use for defining arc weights.
            Defaults to "travel_time".
        directed (bool) -- Whether to treat the network as directed. Defaults
            to True, in which case the one-way designation of roads is taken
            into account. If False, all roads are taken as two-way.
    
    The output graph file simply defines arc-level data since there we have no
    need for node-level data. The arc file includes the following tab-separated
    columns, in order:
        tail -- Node index of arc's tail (origin).
        head -- Node index of arc's head (destination).
        weight -- Numerical weight to use in computing distances.
    
    The node indices in the output file correspond to those used in the
    original .graphml file.
    
    Two-way roads are handled by duplicating every arc in the output file with
    its head and tail reversed.
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
    
    # Find the ID number of the "oneway" attribute
    oneid = None
    for child in root.findall(namespace + "key"):
        if (child.get("for") == "edge") and (child.get("attr.name") == "oneway"):
            oneid = child.get("id")
            break
    
    # Initialize edge weight dictionary, indexed by (tail, head) pair
    edges = dict()
    
    # Find the "graph" tag
    for child in root:
        if child.tag != namespace + "graph":
            continue
        
        # Go through each edge one-at-a-time
        for e in tqdm.tqdm(child.findall(namespace + "edge")):
            
            # Get endpoints
            tail = e.attrib["source"]
            head = e.attrib["target"]
            
            # Initialize edge weight and one-way status
            wt = 0.0
            oneway = False
            
            # Read fields
            for f in e:
                if f.attrib["key"] == wid:
                    wt = f.text
                if (directed == True and f.attrib["key"] == oneid and
                    f.text == "True"):
                    oneway = True
            
            # If the edge is new, create a new dictionary entry
            if (tail, head) not in edges:
                edges[(tail, head)] = wt
            else:
                # Otherwise update the weight if smaller
                if wt < edges[(tail, head)]:
                    edges[(tail, head)] = wt
            
            # If the edge is not one-way, repeat the process for its reverse
            if oneway == False:
                if (head, tail) not in edges:
                    edges[(head, tail)] = wt
                else:
                    if wt < edges[(head, tail)]:
                        edges[(head, tail)] = wt
    
    # Write output file
    with open(arcfile, 'w') as f:
        
        # Write header
        f.write("tail\thead\tweight\n")
        
        # Write array contents to file
        for pair in edges:
            f.write(str(pair[0]) + "\t" + str(pair[1]) + "\t" +
                    str(edges[pair]) + "\n")

#------------------------------------------------------------------------------

def map_node_locations(mapfile, popfile, pnodefile, facfile, fnodefile,
                       popnbrfile=None, facnbrfile=None):
    """Finds the nodes that correspond to population and facility locations.
    
    Positional arguments:
        mapfile (str) -- Path to input .graphml file.
        popfile (str) -- Path to preprocessed population file.
        pnodefile (str) -- Path to population node output file.
        facfile (str) -- Path to preprocessed facility file.
        fnodefile (str) -- Path to facility node output file.
    
    Keyword arguments:
        popnbrfile (str) -- Preprocessed neighboring county population file
            path. Defaults to None.
        facnbrfile (str) -- Preprocessed neighboring county facility file path.
            Defaults to None.
    
    Each population center is snapped to its nearest node in the .graphml file
    based on the geographic coordinates in the popfile. The output pnodefile
    is a tab-separated file containing the following columns:
        id -- Population center index from the population file.
        node -- Corresponding node index from the .graphml file.
    
    If a popnbrfile is given, its included locations are included in the
    population node list (it is assumed here that no neighbor node ID is the
    same as a main county node ID, which should matches the format produced by
    the preprocessing scripts).
    
    The facfile, fnodefile, and facnbrfile arguments are analogous, but are for
    mapping facilities to their nearest nodes.
    """
    
    # Open graph file
    G = ox.load_graphml(mapfile)
    
    t = time.time()
    
    # Get population center coordinates
    (_, pcoord) = _read_popfile(popfile)
    
    # If a neighboring population file is provided, merge its contents
    if popnbrfile != None:
        (_, pcoordnbr) = _read_popfile(popnbrfile)
        pcoord = {**pcoord, **pcoordnbr}
        del pcoordnbr
    
    # Get population center IDs
    pid = list(pcoord.keys())
    pid.sort()
    
    # Find nodes nearest each population center
    print(f"Mapping {len(pcoord)} population nodes.")
    t = time.time()
    pnodes = ox.distance.nearest_nodes(G, [pcoord[i][1] for i in pid],
                                          [pcoord[i][0] for i in pid])
    print(f"Nodes mapped after {time.time()-t:.1} seconds.")
    
    # Write output file
    with open(pnodefile, 'w') as f:
        f.write("id\tnode\n")
        for i in range(len(pid)):
            f.write(str(pid[i]) + "\t" + str(pnodes[i]) + "\n")
    
    del pid
    del pcoord
    
    # Get facility coordinates
    (_, fcoord) = _read_facfile(facfile)
    
    # If a neighboring facility file is provided, merge its contents
    if facnbrfile != None:
        (_, fcoordnbr) = _read_popfile(facnbrfile)
        fcoord = {**fcoord, **fcoordnbr}
        del fcoordnbr
    
    # Get facility IDs
    fid = list(fcoord.keys())
    fid.sort()
    
    # Find nodes nearest each facility
    print(f"Mapping {len(fcoord)} facility nodes.")
    t = time.time()
    fnodes = ox.distance.nearest_nodes(G, [fcoord[i][1] for i in fid],
                                          [fcoord[i][0] for i in fid])
    print(f"Nodes mapped after {time.time()-t:.1} seconds.")
    
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
    
    # Build adjacency list and gather weights from arc file
    adj = [[] for i in range(size)] # adjacency table
    weight = [[] for i in range(size)] # weight table
    with open(arcfile, 'r') as f:
        for line in f:
            # Skip comment line
            if line[0].isdigit() == False:
                continue
            s = line.strip().split()
            # Remap endpoints
            tail = position[int(s[0])]
            head = position[int(s[1])]
            # Multiply weight and cast as integer
            wt = int(factor*float(s[2]))
            # Add data to adjacency and weight lists
            adj[tail].append(head)
            weight[tail].append(wt)
    
    # Recursively convert lists to tuples
    adj = tuple(tuple(adj[i]) for i in range(len(adj)))
    weight = tuple(tuple(weight[i]) for i in range(len(weight)))
    
    return (adj, weight, onodes, dnodes)

#==============================================================================
# Shortest Path Computation
#==============================================================================

def _dijkstra_single(adj, weight, s, dnodes):
    """Carries out Dijkstra's algorithm for a single source.
    
    Positional arguments:
        adj (tuple(tuple(int))) -- Tuple of tuples of out-neighbors. adj[i][j]
            indicates the jth out-neighbor of the ith node.
        weight (tuple(tuple(float))) -- Tuple of tuples of arc weights.
            weight[i][j] indicates the weight of the arc from node i to the
            jth out-neighbor of i (i.e. the weight of the arc from i to
            adj[i][j]).
        s (int) -- Single origin node.
        dnodes (list(int)) -- List of destination nodes.
    
    Returns:
        (list(float)) -- List of distances from the origin to all destinations,
            respectively.
    
    Note that this script assumes zero-indexed consecutive node IDs, in which
    case adjacency list positions correspond to tail node IDs.
    """
    
    # Initialize distance list
    dist = [None for i in range(len(adj))]
    dist[s] = 0 # origin's distance to self is zero
    
    # Initialize explored node set
    explored = set([s])
    
    # Initialize unexplored destination set
    udest = set(dnodes)
    if s in udest:
        udest.remove(s)
    
    # Initialize tentative distance min-priority queue (dist, node)
    tentative = []
    for i in range(len(adj[s])):
        heapq.heappush(tentative, (weight[s][i], adj[s][i]))
        dist[adj[s][i]] = weight[s][i]
    
    # Enter the main Dijkstra loop
    while len(tentative) > 0 and len(udest) > 0:
        
        # Pick the minimum-distance tentative node
        (udist, u) = heapq.heappop(tentative)
        
        # Skip if already explored
        if u in explored:
            continue
        
        # Mark as explored (and remove from destination list if needed)
        explored.add(u)
        if u in udest:
            udest.remove(u)
        
        # Update all neighboring tentative distances
        for i in range(len(adj[u])):
            v = adj[u][i]
            newdist = udist + weight[u][i]
            if dist[v] == None or newdist < dist[v]:
                dist[v] = newdist
                heapq.heappush(tentative, (newdist, v))
    
    # Return the distances
    return [dist[dnodes[i]] for i in range(len(dnodes))]

#------------------------------------------------------------------------------

def distance_table(arcfile, pnodefile, fnodefile, factor=10):
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
        (list(list(float))) -- Table of origin-to-destination pairwise
            distances. Element [i][j] of this table is the distance from the
            ith origin defined in pnodefile to the jth destination defined in
            fnodefile.
    """
    
    # Generate adjacency list representation
    (adj, weight, onodes, dnodes) = _generate_adjacency_list(arcfile,
                                           pnodefile, fnodefile, factor=factor)
    
    # Initialize distance table
    dist = [[None for j in range(len(dnodes))] for i in range(len(onodes))]
    
    # Fill rows of table one source at a time
    for i in tqdm.tqdm(range(len(onodes))):
        dist[i] = _dijkstra_single(adj, weight, onodes[i], dnodes)
    
    ### Eventually edit this to write distances to a log file as it goes.
    
    # Undo multiplicative factor
    for i in range(len(onodes)):
        for j in range(len(dnodes)):
            dist[i][j] /= factor
    
    return dist

#==============================================================================
# Main Frontend
#==============================================================================

def distance_file(arcfile, pnodefile, fnodefile, distfile, factor=10,
                  multiplier=float(1.0/60)):
    """Generates a complete origin-to-destination distance file.
    
    Positional arguments:
        arcfile (str) -- Path to the network's arc definition file, which
            should contain the tail node, head node, and travel time for each
            arc.
        pnodefile (str) -- Path to the network's population node definition
            file, which should contain a list of origin nodes.
        fnodefile (str) -- Path to the network's facility node definition file,
            which should contain a list of destination nodes.
        distfile (str) -- Path to output distance file.
    
    Keyword arguments:
        factor (float) -- Factor by which to multiply all arc weights before
            casting them as integers during the distance computations. Defaults
            to 10, which is appropriate for an arc file with weights that have
            only one digit past the decimal place.
        multiplier (float) -- Factor by which to multiply the final distance
            values. Defaults to 1/60, which is appropriate for generating
            travel times in minutes given arc weights in seconds.
    
    The format of the output distance file is a TSV containing the following
    columns:
        pid -- Index of a population center from the population file.
        fid -- Index of a facility from the facility file.
        time -- Travel time from the population center to the facility.
    """
    
    # Gather population center IDs
    pid = []
    with open(pnodefile, 'r') as f:
        for line in f:
            if line[0].isdigit() == False:
                continue
            pid.append(line.strip().split()[0])
    
    # Gather facility IDs
    fid = []
    with open(fnodefile, 'r') as f:
        for line in f:
            if line[0].isdigit() == False:
                continue
            fid.append(line.strip().split()[0])
    
    # Generate travel times
    dist = distance_table(arcfile, pnodefile, fnodefile, factor=factor)
    
    # Write distance file
    with open(distfile, 'w') as f:
        f.write("pid\tfid\ttime\n")
        for i in range(len(pid)):
            for j in range(len(fid)):
                # Apply multiplicative factor to distance
                d = multiplier*dist[i][j]
                f.write(str(pid[i]) + "\t" + str(fid[j]) + "\t" + str(d) + "\n")

#==============================================================================
# Execution
#==============================================================================

# Comment or uncomment the function calls below to process each location.

# Generate Santa Clara place list
#santa_clara_mapfile = os.path.join("..", "maps", "santa_clara", "santa_clara_map.graphml")
#santa_clara_arcfile = os.path.join("..", "graphs", "santa_clara", "santa_clara_arcs.tsv")
#santa_clara_popfile = os.path.join("..", "processed", "santa_clara", "santa_clara_pop.tsv")
#santa_clara_popfile_nbr = os.path.join("..", "processed", "santa_clara", "santa_clara_pop_nbr.tsv")
#santa_clara_facfile = os.path.join("..", "processed", "santa_clara", "santa_clara_fac.tsv")
#santa_clara_facfile_nbr = os.path.join("..", "processed", "santa_clara", "santa_clara_fac_nbr.tsv")
#santa_clara_pnodefile = os.path.join("..", "graphs", "santa_clara", "santa_clara_popnodes.tsv")
#santa_clara_fnodefile = os.path.join("..", "graphs", "santa_clara", "santa_clara_facnodes.tsv")
#santa_clara_distfile = os.path.join("..", "processed", "santa_clara", "santa_clara_dist.tsv")
#santa_clara_places = [{"county": c, "state": "California"} for c in SANTA_CLARA_NEIGHBORS] + [{"county": "Santa Clara", "state": "California"}]

# Download: 501.8285789489746 seconds
# Edge speeds: 34.205276012420654 seconds
#download_map_places(santa_clara_mapfile, santa_clara_places)

#graphml_to_tsv(santa_clara_mapfile, santa_clara_arcfile, weight="travel_time", directed=False)
#map_node_locations(santa_clara_mapfile, santa_clara_popfile, santa_clara_pnodefile, santa_clara_facfile, santa_clara_fnodefile, popnbrfile=santa_clara_popfile_nbr, facnbrfile=santa_clara_facfile_nbr)
#distance_file(santa_clara_arcfile, santa_clara_pnodefile, santa_clara_fnodefile, santa_clara_distfile, factor=10, multiplier=1.0/60)

# Generate Polk place list
polk_mapfile = os.path.join("..", "maps", "polk", "polk_map.graphml")
#polk_arcfile = os.path.join("..", "graphs", "polk", "polk_arcs.tsv")
#polk_popfile = os.path.join("..", "processed", "polk", "polk_pop.tsv")
#polk_popfile_nbr = os.path.join("..", "processed", "polk", "polk_pop_nbr.tsv")
#polk_facfile = os.path.join("..", "processed", "polk", "polk_fac.tsv")
#polk_facfile_nbr = os.path.join("..", "processed", "polk", "polk_fac_nbr.tsv")
#polk_pnodefile = os.path.join("..", "graphs", "polk", "polk_popnodes.tsv")
#polk_fnodefile = os.path.join("..", "graphs", "polk", "polk_facnodes.tsv")
#polk_distfile = os.path.join("..", "processed", "polk", "polk_dist.tsv")
#polk_places = [{"county": c, "state": "Florida"} for c in POLK_NEIGHBORS] + [{"county": "Polk", "state": "Florida"}]

# Download: 267.9547348022461 seconds
# Edge speeds: 24.942673921585083 seconds
#download_map_places(polk_mapfile, polk_places)
