"""COVID-19 Vaccine accessibility project preprocessing scripts.

The code below includes various scripts for collecting data from a variety of
assorted raw data files and collecting them into summary files with a standard
format. The Google Maps API is used to convert addresses to coordinates where
needed.

Since each location included in the study organizes its data slightly
differently, a different function has been defined to perform the preprocessing
for each location.
"""

import os.path
import requests

#==============================================================================
# Global Constants
#==============================================================================

# Population center file header (including column labels)
POP_HEADER = "id\tname\tlat\tlon\tpop\tvacc\tadi\t\n"

# Facility file header (including column labels)
FAC_HEADER = "id\tname\tlat\tlon\tcap\t\n"

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
    
    # Define location-specific file names
    case_file = os.path.join("chicago",
                           "COVID-19_Cases__Tests__and_Deaths_by_ZIP_Code.csv")
    fac_file = os.path.join("chicago", "COVID-19_Vaccination_Locations.csv")
    vacc_file = os.path.join("chicago","COVID-19_Vaccinations_by_ZIP_Code.csv")
    adi_file = os.path.join("chicago", "IL_2020_ADI_9 Digit Zip Code_v3.2.csv")
    
    ###
    
    # Write population output file
    with open(popfile, 'w') as f:
        f.write(POP_HEADER)
        ###
    
    # Write facility output file
    with open(facfile, 'w') as f:
        f.write(FAC_HEADER)
        ###

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
    
    # Define location-specific file names
    adi_file = os.path.join("santa_clara",
                            "CA_2020_ADI_Census Block Group_v3.2.csv")
    census_file = os.path.join("santa_clara", "CensusTract2020.csv")
    vacc_file = os.path.join("santa_clara",
             "COVID-19_Vaccination_among_County_Residents_by_Census_Tract.csv")
    
    ###
    
    # Write population output file
    with open(popfile, 'w') as f:
        f.write(POP_HEADER)
        ###
    
    # Write facility output file
    with open(facfile, 'w') as f:
        f.write(FAC_HEADER)
        ###

#==============================================================================
# Execution
#==============================================================================

# Comment or uncomment the function calls below to process each location.
process_chicago()
process_santa_clara()
