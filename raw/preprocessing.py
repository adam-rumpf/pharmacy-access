"""COVID-19 Vaccine accessibility project preprocessing scripts.

The code below includes various scripts for collecting data from a variety of
assorted raw data files and collecting them into summary files with a standard
format.

Since each location included in the study organizes its data slightly
differently, a different function has been defined to perform the preprocessing
for each location.
"""

import os.path
import re

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

def point_to_coords(point):
    """Extracts the latitude/longitude from a GIS POINT object.
    
    Positional arguments:
        point (str) -- String of the format "POINT (LON LAT)".
    
    Returns:
        (tuple(float)) -- Latitude/longitude.
    """
    
    s = re.findall("[-.\d]+", point)
    return (float(s[1]), float(s[0]))

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
    
    # Initialize population center dictionary
    pdic = dict()
    
    # Gather ZIP code locations and populations
    with open(case_file, 'r') as f:
        
        for line in f:
            
            # Skip comment line
            if line[0].isdigit() == False:
                continue
            
            s = line.strip().split(',')
            zc = int(s[0]) # current row's ZIP code
            
            # Initialize empty entry for a new ZIP code
            if zc not in pdic:
                pdic[zc] = [0 for i in range(5)]
            
            # Gather coordinates and population
            pdic[zc][0], pdic[zc][1] = point_to_coords(s[-1])
            pdic[zc][2] = max(pdic[zc][2], int(s[18]))
    
    # Gather vaccination rates
    with open(vacc_file, 'r') as f:
        
        for line in f:
            
            # Skip comment line
            if line[0].isdigit() == False:
                continue
            
            s = line.strip().split(',')
            zc = int(s[0]) # current row's ZIP code
            
            # Gather cumulative vaccinations
            ### Fix this, since this data does not include vaccination rates
            ### per ZIP code
            try:
                pdic[zc][3] += int(s[4])
            except ValueError:
                pass
    
    # Gather ADI rankings
    with open(adi_file, 'r') as f:
        
        # Initialize dictionary to group 9-digit ZIP codes by 5-digit code
        adi = dict([(zc, [0, 0]) for zc in pdic])
        
        for line in f:
            
            # Skip comment line
            if line[1].isdigit() == False:
                continue
            
            s = line.replace('"', '').strip().split(',')
            
            # Take 5-digit header
            zc = int(s[0][:5])
            if zc not in pdic:
                continue
            
            # Add ADI ranking to tally
            try:
                adi[zc][0] += int(s[4])
                adi[zc][1] += 1
            except ValueError:
                pass
        
        # Average 9-digit values across 5-digit codes
        for zc in pdic:
            if adi[zc][1] > 0:
                pdic[zc][4] = adi[zc][0]/adi[zc][1]
    
    # Write population output file
    with open(popfile, 'w') as f:
        f.write(POP_HEADER)
        sk = sorted(pdic.keys())
        for i in range(len(sk)):
            line = str(i) + '\t' + str(sk[i]) + '\t'
            for item in pdic[sk[i]]:
                line += str(item) + '\t'
            f.write(line + '\n')
    
    # Initialize facility dictionary
    fdic = dict()
    
    # Gather vaccination facility locations
    with open(fac_file, 'r') as f:
        
        for line in f:
            
            # Skip comment line
            if line[0].isdigit() == False:
                continue
            
            s = line.strip().split(',')
            fi = int(s[0]) # current row's facility ID
            
            # Initialize empty entry for a new facility
            if fi not in fdic:
                fdic[fi] = [0 for i in range(3)]
            
            # Gather coordinates
            try:
                fdic[fi][0], fdic[fi][1] = point_to_coords(s[-1])
            except IndexError:
                # If no POINT is given, search for the coordinates
                ### Fill in auto coordinate searching if no coordinates.
                pass
            
            # Gather capacity
            ### Find a way to measure capacity.
            fdic[fi][2] = 1
    
    # Write facility output file
    with open(facfile, 'w') as f:
        f.write(FAC_HEADER)
        sk = sorted(fdic.keys())
        for i in range(len(sk)):
            line = str(i) + '\t' + str(sk[i]) + '\t'
            for item in fdic[sk[i]]:
                line += str(item) + '\t'
            f.write(line + '\n')

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
    
    # Initialize population center dictionary
    pdic = dict()
    
    # Gather tract locations
    with open(census_file, 'r') as f:
        
        for line in f:
            
            # Skip comment line
            if line.strip()[-1].isdigit() == False:
                continue
            
            s = line.strip().split(',')
            tid = int(s[-13]) # current row's tract ID
            
            # Initialize empty entry for a new tract ID
            if tid not in pdic:
                pdic[tid] = [0 for i in range(5)]
            
            # Gather coordinates
            pdic[tid][0] = float(s[-6])
            pdic[tid][1] = float(s[-5])
    
    # Write population output file
    with open(popfile, 'w') as f:
        f.write(POP_HEADER)
        sk = sorted(pdic.keys())
        for i in range(len(sk)):
            line = str(i) + '\t' + str(sk[i]) + '\t'
            for item in pdic[sk[i]]:
                line += str(item) + '\t'
            f.write(line + '\n')
    
    # Write facility output file
    with open(facfile, 'w') as f:
        f.write(FAC_HEADER)
        ###

#==============================================================================
# Execution
#==============================================================================

# Comment or uncomment the function calls below to process each location.
#process_chicago()
process_santa_clara()
