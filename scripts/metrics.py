"""COVID-19 Vaccine accessibility project accessibility metrification scripts.

The code below defines various functions that work with the preprocessed data
files, all of which should have a standardized format. Various accessibility
metrics can be computed, and various distance measurements can be used.
"""

import os.path

import geopy
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
# Metric Definitions
#==============================================================================

def gravity_metric(outfile, popfile, facfile, distfile=None, beta=1.0):
    """Computes a table of gravitational metrics for a given community.
    
    Positional arguments:
        outfile (str) -- Output table file path.
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
    """
    
    pass

#==============================================================================
# Execution
#==============================================================================

# Comment or uncomment the function calls below to process each location.
gravity_metric(os.path.join("..", "results", "santa_clara", "santa_clara_accessibility_gravity_1-00.tsv"), os.path.join("..", "processed", "santa_clara", "santa_clara_pop.tsv"), os.path.join("..", "processed", "santa_clara", "santa_clara_fac.tsv"), distfile=None, beta=1.0)
