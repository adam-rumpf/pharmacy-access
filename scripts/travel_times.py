"""COVID-19 Vaccine accessibility project data processing scripts.

The code below includes various scripts for generating required data from the
preprocessed data files, notably including origin/destination travel time
matrices.
"""

import traveltimepy
from tqdm import tqdm

#==============================================================================
# Global Constants
#==============================================================================

# Distance center file header (including column labels)
DIST_HEADER = "pid\tfid\tpftime\tfptime\t\n"

#==============================================================================
# Common Functions
#==============================================================================

def travel_time(origin, destination):
    """Computes the travel time between a given pair of coordinates.
    
    Positional arguments:
        origin (tuple(float)) -- Origin latitude/longitude tuple.
        destination (tuple(float)) -- Destination latitude/longitude tuple.
    
    Returns:
        (float) -- Travel time (minutes).
    """
    
    ### Add options for things like travel mode and time of day
    
    pass

#==============================================================================
# Location-Specific Processing Scripts
#==============================================================================

###
### Use tqdm for OD pair progress bars

#==============================================================================
# Execution
#==============================================================================

# Comment or uncomment the function calls below to process each location.
###
