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
    
    We then compute the accessibility metric for each facility as a distance-
    weighted sum of "capacity-to-crowdedness" ratios for each facility
    consisting of a sum of all capacity to crowdedness ratios over all
    facilities, each scaled by a weight that decays over distance
        
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
    with open(foutfile, 'w') as f:
        with open(facfile, 'r') as g:
            for line in g:
                
                # Copy header with a new metric column
                if line[0].isdigit() == False:
                    f.write(line.strip() + "\tcrowding\t\n")
                    continue
                
                # Append the metric line to the following rows
                i = int(line.split('\t')[0]) # current facility index
                m = fmet[i] # current facility metric
                f.write(line.strip() + '\t' + str(m) + "\t\n")
    
    # Write population output file with a new accessibility metric column
    print("Writing population metric file.")
    with open(poutfile, 'w') as f:
        with open(popfile, 'r') as g:
            for line in g:
                
                # Copy header with a new metric column
                if line[0].isdigit() == False:
                    f.write(line.strip() + "\taccess\t\n")
                    continue
                
                # Append the metric line to the following rows
                i = int(line.split('\t')[0]) # current population index
                m = pmet[i] # current population metric
                f.write(line.strip() + '\t' + str(m) + "\t\n")

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
    pop = dict() # dictionary of populations by popfile index
    pcoord = dict() # dictionary of latitutde/longitude pairs by popfile index
    with open(popfile, 'r') as f:
        
        for line in f:
            
            # Skip the comment line
            if line[0].isdigit() == False:
                continue
            
            # Get numbers from line
            s = line.strip().split('\t')
            pop[int(s[0])] = int(s[POP_POP])
            pcoord[int(s[0])] = (float(s[POP_LAT]), float(s[POP_LON]))
    
    # Read facility file
    cap = dict() # dictionary of facility capacities by facfile index
    fcoord = dict() # dictionary of latitude/longitude pairs by facfile index
    with open(facfile, 'r') as f:
        
        for line in f:
            
            # Skip the comment line
            if line[0].isdigit() == False:
                continue
            
            # Get numbers from line
            s = line.strip().split('\t')
            cap[int(s[0])] = float(s[FAC_CAP])
            fcoord[int(s[0])] = (float(s[FAC_LAT]), float(s[FAC_LON]))
    
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
# Execution
#==============================================================================

# Comment or uncomment the function calls below to process each location.
gravity_metric(os.path.join("..", "results", "santa_clara", "santa_clara_pop_gravity_1-00.tsv"), os.path.join("..", "results", "santa_clara", "santa_clara_fac_gravity_1-00.tsv"), os.path.join("..", "processed", "santa_clara", "santa_clara_pop.tsv"), os.path.join("..", "processed", "santa_clara", "santa_clara_fac.tsv"), distfile=None, beta=1.0)
