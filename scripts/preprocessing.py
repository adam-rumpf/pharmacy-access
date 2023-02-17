"""COVID-19 Vaccine accessibility project preprocessing scripts.

The code below includes various scripts for collecting data from a variety of
assorted raw data files and collecting them into summary files with a standard
format.

Since each location included in the study organizes its data slightly
differently, a different function has been defined to perform the preprocessing
for each location.
"""

import csv
import os.path
import re
import time

import geopy
import tqdm

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

#------------------------------------------------------------------------------

def address_to_coords(address, geocoder=None):
    """Finds the latitude/longitude coordinates from a given address string.
    
    Positional arguments:
        address (str) -- The address string to search for.
    
    Optional keyword arguments:
        geocoder (Nominatim) -- A geopy geocoder object. Defaults to None, in
            which case a new geocoder is temporarily created.
    
    Returns:
        (tuple(float)) -- Latitude/longitude.
    """
    
    if geocoder == None:
        # Try to get email address
        try:
            with open("email.txt", 'r') as f:
                user_agent = f.read().strip().split()[0]
        except FileNotFoundError:
            user_agent = "user_agent"
        
        # Initialize a geocoder
        geocoder = geopy.geocoders.Nominatim(user_agent=user_agent)
    
    location = geocoder.geocode(address)
    
    return (location.latitude, location.longitude)

#------------------------------------------------------------------------------

def pharmacy_table_coords(infile, outfile=None, user_agent=None):
    """Augments a pharmacy CSV file by adding latitude/longitude fields.
    
    Positional arguments:
        infile (str) -- Pharmacy file path. This is a custom input file based
            on data gathered by hand.
    
    Optional keyword arguments:
        outfile (str) -- Output file path for the augmented pharmacy file.
            Default None, in which case the original file is overwritten.
        user_agent (str) -- User agent string to give to the geocoder. Default
            None, in which case this script will attempt to read an email
            address from a local "email.txt" file. If none is found, defaults
            to the string "user_agent".
    """
    
    if outfile == None:
        outfile = infile
    
    # Try to get email address
    if user_agent == None:
        try:
            with open("email.txt", 'r') as f:
                user_agent = f.read().strip().split()[0]
        except FileNotFoundError:
            user_agent = "user_agent"
    
    # Initialize a geocoder
    gc = geopy.geocoders.Nominatim(user_agent=user_agent)
    
    # Read pharmacy fields from CSV file header
    with open(infile, 'r') as f:
        fields = f.readline().strip().split(',')
    
    # Add new lat/lon fields
    fields += ["latitude", "longitude"]
    
    # Read contents of pharmacy CSV file
    with open(infile, 'r') as f:
        
        # Create a list of dictionaries for each row
        pdic = list(csv.DictReader(f, delimiter=',', quotechar='"'))
        
        # Look up the address on each row
        i = 0
        retries = [] # list of rows that caused errors
        for row in tqdm.tqdm(pdic):
            i += 1
            
            # Get the address string
            address = (row["address line 1"] + ", " + row["city"] + ", "
                       + row["state"] + " " + row["zipcode"])
            
            # Geocode the address
            try:
                (lat, lon) = address_to_coords(address, gc)
            except geopy.exc.GeocoderUnavailable:
                retries.append(i)
                (lat, lon) = (None, None)
            except AttributeError:
                retries.append(i)
                (lat, lon) = (None, None)
            time.sleep(1) # wait 1 second between Nominatim requests
            
            # Add address fields to dictionary line
            row["latitude"] = lat
            row["longitude"] = lon
    
    if len(retries) > 0:
        print("Errors in the following rows:\n" +
              ", ".join([str(i) for i in retries]))
    
    del gc
    
    # Write the new pharmacy CSV file
    with open(outfile, 'w', newline="") as f:
        writer = csv.DictWriter(f, delimiter=',', quotechar='"',
                                fieldnames=fields)
        writer.writeheader()
        for row in pdic:
            writer.writerow(row)

#==============================================================================
# Location-Specific Preprocessing Scripts
#==============================================================================

def process_chicago(popfile=os.path.join("..", "processed", "chicago",
                                         "chicago_pop.tsv"),
                    facfile=os.path.join("..", "processed", "chicago",
                                         "chicago_fac.tsv")):
    """Preprocessing scripts for the Chicago data.
    
    Optional keyword arguments:
        popfile (str) -- Population output file path. Defaults to a file in the
            processed/ directory named "chicago_pop.tsv".
        facfile (str) -- Facility output file path. Defaults to a file in the
            processed/ directory named "chicago_fac.tsv".
    """
    
    # Define location-specific file names
    case_file = os.path.join("..", "data", "chicago",
                           "COVID-19_Cases__Tests__and_Deaths_by_ZIP_Code.csv")
    fac_file = os.path.join("..", "data", "chicago",
                            "COVID-19_Vaccination_Locations.csv")
    vacc_file = os.path.join("..", "data", "chicago",
                             "COVID-19_Vaccinations_by_ZIP_Code.csv")
    adi_file = os.path.join("..", "data", "chicago",
                            "IL_2020_ADI_9 Digit Zip Code_v3.2.csv")
    
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
            ### Fix this, since this data do not include vaccination rates
            ### per ZIP code of residence, which leads to many values
            ### greater than 100%.
            try:
                pdic[zc][3] += int(s[4])
            except ValueError:
                pass
            
        # Compute vaccination rates
        ### Currently using a 95% cap as a placeholder.
        for zc in pdic:
            
            if pdic[zc][2] <= 0:
                pdic[zc][3] = 0.0
            else:
                pdic[zc][3] = min(0.95, pdic[zc][3]/pdic[zc][2])
    
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
        index = 0
        for i in range(len(sk)):
            # Skip lines with no coordinates
            if pdic[sk[i]][0] == 0 or pdic[sk[i]][1] == 0:
                continue
            line = str(index) + '\t' + str(sk[i]) + '\t'
            for item in pdic[sk[i]]:
                line += str(item) + '\t'
            f.write(line + '\n')
            index += 1
    
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

def process_santa_clara(popfile=os.path.join("..", "processed", "santa_clara",
                                             "santa_clara_pop.tsv"),
                        facfile=os.path.join("..", "processed", "santa_clara",
                                             "santa_clara_fac.tsv")):
    """Preprocessing scripts for the Santa Clara data.
    
    Optional keyword arguments:
        popfile (str) -- Population output file path. Defaults to a file in the
            processed/ directory named "santa_clara_pop.tsv".
        facfile (str) -- Facility output file path. Defaults to a file in the
            processed/ directory named "santa_clara_fac.tsv".
    """
    
    # Define location-specific file names
    adi_file = os.path.join("..", "data", "santa_clara",
                            "CA_2020_ADI_Census Block Group_v3.2.csv")
    fac_file = os.path.join("..", "data", "santa_clara",
                            "Santa_Clara_County_Pharmacies.csv")
    census_file = os.path.join("..", "data", "santa_clara",
                               "2022_gaz_tracts_06.txt")
    vacc_file = os.path.join("..", "data", "santa_clara",
             "COVID-19_Vaccination_among_County_Residents_by_Census_Tract.csv")
    
    # Initialize population center dictionary
    pdic = dict()
    
    # Gather vaccination rates
    with open(vacc_file, 'r') as f:
        
        for line in f:
            
            # Skip comment line
            if line[0].isdigit() == False:
                continue
            
            s = line.strip().split(',')
            fips = int(s[0]) # current row's FIPS
            
            # Initialize empty entry for a new FIPS
            if fips not in pdic:
                pdic[fips] = [0 for i in range(5)]
            
            # Gather population and vaccinations
            try:
                pdic[fips][2] = int(s[2])
                v = int(s[3]) # Vaccination count
                # Compute vaccination rate, capped at 0.95
                if pdic[fips][2] <= 0:
                    pdic[fips][3] = 0.0
                else:
                    pdic[fips][3] = min(0.95, v/pdic[fips][2])
                    ###
                    # Decide how to handle cases >0.95
            except ValueError:
                pdic[fips][2] = 0
                pdic[fips][3] = 0.0
    
    # Gather tract locations
    with open(census_file, 'r') as f:
        
        for line in f:
            
            # Skip comment line
            if line.strip()[-1].isdigit() == False:
                continue
            
            s = line.strip().split('\t')
            fips = int(s[1]) # current row's FIPS
            
            # Skip tracts not logged in the vaccination file
            if fips not in pdic:
                continue
            
            # Gather coordinates
            pdic[fips][0] = float(s[6])
            pdic[fips][1] = float(s[7])
    
    # Gather ADI rankings
    with open(adi_file, 'r') as f:
        
        # Initialize dictionary to group 12-digit FIPS by 11-digit FIPS
        adi = dict([(fips, [0, 0]) for fips in pdic])
        
        for line in f:
            
            # Skip comment line
            if line[2].isdigit() == False:
                continue
            
            s = line.replace('"', '').strip().split(',')
            
            # Take 11-digit header
            fips = int(s[3][:11])
            if fips not in pdic:
                continue
            
            # Add ADI ranking to tally
            try:
                adi[fips][0] += int(s[2])
                adi[fips][1] += 1
            except ValueError:
                pass
        
        # Average 12-digit values across 11-digit codes
        for fips in pdic:
            if adi[fips][1] > 0:
                pdic[fips][4] = adi[fips][0]/adi[fips][1]
    
    # Write population output file
    with open(popfile, 'w') as f:
        f.write(POP_HEADER)
        sk = sorted(pdic.keys())
        index = 0
        for i in range(len(sk)):
            # Skip lines with no coordinates
            if pdic[sk[i]][0] == 0 or pdic[sk[i]][1] == 0:
                continue
            line = str(index) + '\t' + str(sk[i]) + '\t'
            for item in pdic[sk[i]]:
                line += str(item) + '\t'
            f.write(line + '\n')
            index += 1
    
    # Initialize facility dictionary
    fdic = dict()
    
    # Gather vaccination facility locations
    with open(fac_file, 'r') as f:
        
        # Create a list of dictionaries for each row
        pdic = list(csv.DictReader(f, delimiter=',', quotechar='"'))
        namekey = list(pdic[0].keys())[0] # first key should be pharmacy name
        
        # Look up the address on each row
        for row in pdic:
            
            # Get facility name
            fi = row[namekey]
            
            # Initialize empty entry for a new facility
            fdic[fi] = [0 for i in range(3)]
            
            # Try to get coordinates
            try:
                fdic[fi][0] = row["latitude"]
                fdic[fi][1] = row["longitude"]
            except KeyError:
                fdic[fi][0] = None
                fdic[fi][1] = None
            
            # Get capacity
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

#==============================================================================
# Execution
#==============================================================================

# Comment or uncomment the function calls below to process each location.
#process_chicago()
#process_santa_clara()
#pharmacy_table_coords(os.path.join("..", "data", "santa_clara", "Santa_Clara_County_Pharmacies.csv"))
