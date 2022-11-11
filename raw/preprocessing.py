"""COVID-19 Vaccine accessibility project preprocessing scripts.

The code below includes various scripts for collecting data from a variety of
assorted raw data files and collecting them into summary files with a standard
format. The Google Maps API is used to convert addresses to coordinates where
needed.

Since each location included in the study organizes its data slightly
differently, a different function has been defined to perform the preprocessing
for each location.
"""

import requests

#==============================================================================
# Common Functions
#==============================================================================

def address_to_coords(address):
    """Computes the latitude/longitude of a given address.
    
    Positional arguments:
        address (str) -- Address to search for.
    
    Returns:
        (tuple(float)) -- Latitude/longitude of the given address.
    """
    
    pass

#==============================================================================
# Location-Specific Preprocessing Scripts
#==============================================================================

def process_chicago(popfile="chicago_pop.tsv", facfile="chicago_fac.tsv"):
    """Preprocessing scripts for the Chicago data.
    
    Optional keyword arguments:
        popfile (str) -- Population output file path. Defaults to a local file
            named "chicago_pop.tsv".
        facfile (str) -- Facility output file path. Defaults to a local file
            named "chicago_fac.tsv".
    """
    
    ###
    
    # Write population output file
    with open(popfile, 'w') as f:
        pass
    
    # Write facility output file
    with open(facfile, 'w') as f:
        pass

#------------------------------------------------------------------------------

def process_santa_clara(popfile="santa_clara_pop.tsv",
                        facfile="santa_clara_fac.tsv"):
    """Preprocessing scripts for the Santa Clara data.
    
    Optional keyword arguments:
        popfile (str) -- Population output file path. Defaults to a local file
            named "santa_clara_pop.tsv".
        facfile (str) -- Facility output file path. Defaults to a local file
            named "santa_clara_fac.tsv".
    """
    
    ###
    
    # Write population output file
    with open(popfile, 'w') as f:
        pass
    
    # Write facility output file
    with open(facfile, 'w') as f:
        pass

#==============================================================================
# Execution
#==============================================================================

# Comment or uncomment the function calls below to process each location.
process_chicago()
process_santa_clara()
