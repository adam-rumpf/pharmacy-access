"""COVID-19 Vaccine accessibility project data processing scripts.

The code below includes various scripts for generating required data from the
preprocessed data files, notably including origin/destination travel time
matrices.

Travel times are computed using the Google Maps Distance Matrix API.
"""

import requests

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
    
    ### Add options for thigns like travel mode and time of day
    
    pass

#==============================================================================
# Location-Specific Processing Scripts
#==============================================================================

###

#==============================================================================
# Execution
#==============================================================================

# Comment or uncomment the function calls below to process each location.
###
