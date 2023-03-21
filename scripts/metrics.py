"""COVID-19 Vaccine accessibility project accessibility metrification scripts.

The code below defines various functions that work with the preprocessed data
files, all of which should have a standardized format. Various accessibility
metrics can be computed, and various distance measurements can be used.
"""

import os.path

import geopy.distance
import tqdm

#==============================================================================
# Global Constants
#==============================================================================

# Metric file header (including column labels)
METRIC_HEADER = "id\tmetric\t\n"

# Preprocessed file column numbers
FAC_LAT = 2 # facility latitude
FAC_LON = 3 # facility longitude
FAC_CAP = 4 # facility capacity
POP_LAT = 2 # population latitude
POP_LON = 3 # population longitude
POP_POP = 4 # population total population value

#==============================================================================
# Common Functions
#==============================================================================

def geodesic_distance(p1, p2):
    """Measures the geodesic distance between two coordinate pairs.
    
    Positional arguments:
        p1 (tuple(float)) -- First latitude/longitude pair.
        p2 (tuple(float)) -- Second latitude/longitude pair.
    
    Returns:
        (float) -- Geodesic distance (mi) between the coordinates.
    """
    
    return geopy.distance.geodesic(p1, p2).miles

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

def _augment_file(outfile, infile, column, label, default="-1"):
    """Augments a population or facility file with a new column.
    
    Positional arguments:
        outfile (str) -- Output population/facility file path.
        infile (str) -- Preprocessed population/facility file path.
        column (dict(float)) -- Dictionary of entries for the new table
            column. The keys of the dictionary should correspond to the
            indices in the first column of the input file.
        label (str) -- Label for the new column. Since the output file is
            a tab-separated value table, avoid the use of tab characters.
    
    Optional keyword arguments:
        default (str) -- Default value for rows with no corresponding
            dictionary value. Defaults to "-1".
    
    This is meant for use in adding a column of population accessibility
    metrics to a population file, or facility crowding metrics to a facility
    file. It works by copying each row of the original file and adding a new
    element to the end of each row corresponding to that row's index.
    """
    
    with open(poutfile, 'w') as f:
        with open(popfile, 'r') as g:
            for line in g:
                
                # Copy header with a new column label
                if line[0].isdigit() == False:
                    f.write(line.strip() + "\t" + label + "\t\n")
                    continue
                
                # Append the new column element to the following rows
                i = int(line.split('\t')[0]) # current row index
                if i in column:
                    c = column[i] # dictionary element of current line
                else:
                    c = default # default element for missing indices
                f.write(line.strip() + '\t' + str(c) + "\t\n")

#==============================================================================
# Gravity Metric Scripts
#==============================================================================

def gravity_metric(poutfile, foutfile, popfile, facfile, distfile=None,
                   beta=1.0):
    """Computes a table of gravitational metrics for a given community.
    
    Positional arguments:
        poutfile (str) -- Output population file path.
        foutfile (str) -- Output facility file path.
        popfile (str) -- Preprocessed population file path, which should
            include the coordinates and population of each population center.
        facfile (str) -- Preprocessed facility file path, which should include
            the coordinates and capacity of each vaccination facility.
    
    Optional keyword arguments:
        distfile (str) -- Preprocessed distance file path, which should include
            the travel times between each population center/facility pair.
            Defaults to None, in which case geodesic distances are computed and
            used as needed.
        beta (float) -- Gravitational decay parameter. Defaults to 1.0.
    
    This function implements a gravitational accessibility metric as described
    in Luo and Wang 2003 (doi:10.1068/b29120). We begin by computing a
    distance-weighted "crowdedness" metric for each facility consisting of a
    sum of all populations, each scaled by a weight that decays over distance
        
        V_j = sum_(k=1)^m P_k d_kj^(-beta)
    
    where V_j is the crowdedness metric for facility j, P_k is the population
    of center k, d_kj is the distance (or travel time) from k to j, and beta
    is a parameter that controls how quickly the weight should decay over
    distance (larger beta causes more distant locations to receive less
    weight).
    
    We then compute the accessibility metric for each community as a
    distance-weighted sum of "capacity-to-crowdedness" ratios for each
    facility consisting of a sum of all capacity to crowdedness ratios over
    all facilities, each scaled by a weight that decays over distance
        
        A_i = sum_(j=1)^n (S_j d_ij^(-beta))/V_j
    
    where A_i is the metric of population center i, S_j is the capacity of
    facility j, and the remaining notation is the same as above.
    
    The output files are augmented versions of the input files, with a new
    column of crowdedness or accessibility metrics appended to the original
    tables.
    """
    
    # Call an appropriate subroutine depending on distance specification
    if distfile == None:
        (fmet, pmet) = _gravity_metric_geodesic(popfile, facfile, beta=beta)
    else:
        (fmet, pmet) = _gravity_metric_file(popfile, facfile, distfile,
                                            beta=beta)
    
    # Write facility output file with a new crowdedness metric column
    print("Writing facility metric file.")
    _augment_file(foutfile, facfile, fmet, "crowding")
    
    # Write population output file with a new accessibility metric column
    print("Writing population metric file.")
    _augment_file(poutfile, popfile, pmet, "access")

#------------------------------------------------------------------------------

def _gravity_metric_geodesic(popfile, facfile, beta=1.0, speed=45.0):
    """The geodesic distance version of gravity_metric.
    
    This script computes gravity metrics using geodesic distances computed as-
    needed. It is automatically called by gravity_metric when a distance file
    is not provided. It returns dictionaries of facility and population metrics
    rather than writing the results directly to a file.
    
    Geodesic distances are used to compute travel times in minutes, assuming a
    constant speed (in mph).
    """
    
    # Convert speed to miles per minute
    speed /= 60.0
    
    # Read population file
    (pop, pcoord) = _read_popfile(popfile)
    
    # Read facility file
    (cap, fcoord) = _read_facfile(facfile)
    
    # Compute the facility crowdedness metrics
    fmet = dict() # facility crowdedness metrics by facfile index
    print("Computing facility crowdedness metrics using geodesic distance.")
    for j in tqdm.tqdm(cap):
        fmet[j] = 0.0
        for k in pop:
            d = geodesic_distance(pcoord[k], fcoord[j])/speed # time (minutes)
            fmet[j] += pop[k]*(d**(-beta))
    
    # Compute the population accessibility metrics
    pmet = dict() # population accessibility metrics by popfile index
    print("Computing population accessibility metrics using geodesic distance.")
    for i in tqdm.tqdm(pop):
        pmet[i] = 0.0
        for j in cap:
            d = geodesic_distance(pcoord[i], fcoord[j])/speed # time (minutes)
            pmet[i] += (cap[j]*(d**(-beta)))/fmet[j]
    
    # Return facility and population metric dictionaries
    return (fmet, pmet)

#------------------------------------------------------------------------------

def _gravity_metric_file(popfile, facfile, distfile, beta=1.0):
    """The distance file version of gravity_metric.
    
    This script computes gravity metrics using a predefined distance file. It
    is automatically called by gravity_metric when a distance file is provided.
    It returns dictionaries of facility and population metrics rather than
    writing the results directly to a file.
    """
    
    ###
    
    # Return facility and population metric dictionaries
    return (None, None)

#==============================================================================
# 2SFCA Metric Scripts
#==============================================================================

def fca_metric(poutfile, foutfile, popfile, facfile, distfile=None,
               cutoff=30.0):
    """Computes a table of 2SFCA metrics for a given community.
    
    Positional arguments:
        poutfile (str) -- Output population file path.
        foutfile (str) -- Output facility file path.
        popfile (str) -- Preprocessed population file path, which should
            include the coordinates and population of each population center.
        facfile (str) -- Preprocessed facility file path, which should include
            the coordinates and capacity of each vaccination facility.
    
    Optional keyword arguments:
        distfile (str) -- Preprocessed distance file path, which should include
            the travel times between each population center/facility pair.
            Defaults to None, in which case geodesic distances are computed and
            used as needed.
        cutoff (float) -- Travel time cutoff for defining catchment areas (in
            minutes). Defaults to 30.0.
    
    This function implements the two-step floating catchment area (2SFCA)
    metric as described in Luo and Wang 2003 (doi:10.1068/b29120). We begin by
    computing a "crowdedness" metric for each facility by dividing its
    capacity by the population within that facility's catchment area
        
        R_j = S_j/(sum_{k:d_kj<=d_0} P_k)
    
    where R_j is the crowdedness metric for facility j, P_k is the population
    of center k, d_kj is the distance (or travel time) from k to j, and d_0 is
    a distance (or travel time) cutoff used to define catchment areas. A
    larger value of d_0 results in larger catchment area definitions.
    
    We then compute the accessibiility metric for each community as a sum of
    crowdedness metrics for all facilities within the community's catchment
    area
        
        A_i = sum_{k:d_ij<=d_0} R_j
    
    where A_i is the metric of population center i and the remaining notation
    is the same as above. The metric is called "two-step" because of how R_j
    and A_i are computed in a two-step process.
    
    The output files are augmented versions of the input files, with a new
    column of crowdedness or accessibility metrics appended to the original
    tables.
    """
    
    # Call an appropriate subroutine depending on distance specification
    if distfile == None:
        (fmet, pmet) = _fca_metric_geodesic(popfile, facfile, cutoff=cutoff)
    else:
        (fmet, pmet) = _fca_metric_file(popfile, facfile, distfile,
                                            cutoff=cutoff)
    
    # Write facility output file with a new crowdedness metric column
    print("Writing facility metric file.")
    _augment_file(foutfile, facfile, fmet, "crowding")
    
    # Write population output file with a new accessibility metric column
    print("Writing population metric file.")
    _augment_file(poutfile, popfile, pmet, "access")

#------------------------------------------------------------------------------

def _fca_metric_geodesic(popfile, facfile, cutoff=30.0, speed=45.0):
    """The geodesic distance version of fca_metric.
    
    This script computes 2SFCA metrics using geodesic distances computed as-
    needed. It is automatically called by fca_metric when a distance file is
    not provided. It returns dictionaries of facility and population metrics
    rather than writing the results directly to a file.
    
    Geodesic distances are used to compute travel times in minutes, assuming a
    constant speed (in mph).
    """
    
    # Convert speed to miles per minute
    speed /= 60.0
    
    # Read population file
    (pop, pcoord) = _read_popfile(popfile)
    
    # Read facility file
    (cap, fcoord) = _read_facfile(facfile)
    
    # Compute the facility crowdedness metrics
    fmet = dict() # facility crowdedness metrics by facfile index
    print("Computing facility crowdedness metrics using geodesic distance.")
    for j in tqdm.tqdm(cap):
        fmet[j] = 0.0
        for k in pop:
            d = geodesic_distance(pcoord[k], fcoord[j])/speed # time (minutes)
            if d <= cutoff:
                fmet[j] += pop[k]
    
    # Compute the population accessibility metrics
    pmet = dict() # population accessibility metrics by popfile index
    print("Computing population accessibility metrics using geodesic distance.")
    for i in tqdm.tqdm(pop):
        pmet[i] = 0.0
        for j in cap:
            d = geodesic_distance(pcoord[i], fcoord[j])/speed # time (minutes)
            if d <= cutoff:
                pmet[i] += fmet[j]
    
    # Return facility and population metric dictionaries
    return (fmet, pmet)

#------------------------------------------------------------------------------

def _fca_metric_file(popfile, facfile, distfile, cutoff=30.0):
    """The distance file version of gravity_metric.
    
    This script computes gravity metrics using a predefined distance file. It
    is automatically called by gravity_metric when a distance file is provided.
    It returns dictionaries of facility and population metrics rather than
    writing the results directly to a file.
    """
    
    ###
    
    # Return facility and population metric dictionaries
    return (None, None)

#==============================================================================
# Execution
#==============================================================================

# Comment or uncomment the function calls below to process each location.
#gravity_metric(os.path.join("..", "results", "santa_clara", "santa_clara_pop_gravity_1-00.tsv"), os.path.join("..", "results", "santa_clara", "santa_clara_fac_gravity_1-00.tsv"), os.path.join("..", "processed", "santa_clara", "santa_clara_pop.tsv"), os.path.join("..", "processed", "santa_clara", "santa_clara_fac.tsv"), distfile=None, beta=1.0)
fca_metric(os.path.join("..", "results", "santa_clara", "santa_clara_pop_cutoff_30.tsv"), os.path.join("..", "results", "santa_clara", "santa_clara_fac_cutoff_30.tsv"), os.path.join("..", "processed", "santa_clara", "santa_clara_pop.tsv"), os.path.join("..", "processed", "santa_clara", "santa_clara_fac.tsv"), distfile=None, beta=1.0)
