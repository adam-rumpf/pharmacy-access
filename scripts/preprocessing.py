"""COVID-19 Vaccine accessibility project preprocessing scripts.

The code below includes various scripts for collecting data from a variety of
assorted raw data files and collecting them into summary files with a standard
format.

Since each location included in the study organizes its data slightly
differently, a different function has been defined to perform the preprocessing
for each location.
"""

import configparser
import csv
import os.path
import re
import time

import geopy
import geopy.distance
import tqdm

#==============================================================================
# Global Constants
#==============================================================================

# Population center file header (including column labels)
POP_HEADER = "id\tname\tlat\tlon\tpop\tvacc\tadi\tsvi\turban\tpov150\tnoveh\n"

# Shorter population center file header (for neighbor files)
POP_HEADER_SHORT = "id\tname\tlat\tlon\tpop\n"

# Facility file header (including column labels)
FAC_HEADER = "id\tname\tlat\tlon\tcap\t\n"

# Counties neighboring Santa Clara
SANTA_CLARA_NEIGHBORS = ["Alameda", "Merced", "Monterey", "San Benito",
                         "San Francisco", "San Joaquin", "San Mateo",
                         "Santa Cruz", "Stanislaus"]

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

#------------------------------------------------------------------------------

def address_test(fname, tol=0.09469697, outfile=None):
    """Checks for missing coordinates and likely address duplicates.
    
    Positional arguments:
        fname (str) -- Path to an input file to test. This should be a raw
            pharmacy location file including latitude and longitude fields.
    
    Optional keyword arguments:
        tol (float) -- Distance tolerance (in miles). Defaults to 0.09469697,
            which is approximately 500 ft. Pairs of coordinates within this
            distance tolerance will be flagged as likely duplicates.
        outfile (str) -- Path to an output report file. Defaults to None, in
            which case the results of the search are printed to the terminal.
    
    This script is meant for use in preprocessing a raw pharmacy location file
    to find entries with missing coordinates or pairs of likely duplicate
    locations.
    """
    
    # Define standard location file field names
    latfield = "latitude"
    lonfield = "longitude"
    namefield = "pharmacy name"
    a1field = "address line 1"
    a2field = "address line 2"
    cityfield = "city"
    statefield = "state"
    zipfield = "zipcode"
    
    # Initialize report string
    s = ("Testing file: '" + 
         str(os.path.basename(fname)) + "'" + "\nTolerance: " + str(tol) +
         " mi (" + f"{5280*tol:.2f}" + " ft)")
    
    # Read file into a list of dictionaries
    with open(fname, 'r') as f:
        dic = list(csv.DictReader(f, delimiter=',', quotechar='"'))
    
    # Verify that all required fields are present
    if (latfield not in dic[0] or lonfield not in dic[0] or
        namefield not in dic[0] or a1field not in dic[0] or
        a2field not in dic[0] or cityfield not in dic[0] or
        statefield not in dic[0] or zipfield not in dic[0]):
        print("input file missing required fields")
        return
    
    # Process every unique pair of rows
    print("Searching address pairs.")
    missing = 0 # running total of missing coordinates
    dupes = 0 # running total of likely duplicates
    skip = False # whether to skip the current outer loop
    for i in tqdm.tqdm(range(len(dic))):
        if skip == True:
            skip = False
            continue
        for j in range(i+1, len(dic)):
            
            # Get first coordinates (if they exist)
            try:
                c1 = (float(dic[i][latfield]), float(dic[i][lonfield]))
            except ValueError:
                # Report missing coordinates
                missing += 1
                s += "\n\nRow " + str(i+2) + "\nMissing coordinates"
                
                # Log address and continue
                a1 = (dic[i][namefield] + ", " + dic[i][a1field] + ", " +
                      dic[i][a2field] + ", " + dic[i][cityfield] + ", " +
                      dic[i][statefield] + ", " + dic[i][zipfield])
                s += "\n" + a1
                
                # Skip the rest of the current outer loop
                skip = True
                break
            
            # Get second coordinates (if they exist)
            try:
                c2 = (float(dic[j][latfield]), float(dic[j][lonfield]))
            except ValueError:
                # Skip missing j coordinates
                continue
            
            # Log pairs with distance below tolerance
            dist = geodesic_distance(c1, c2)
            if dist <= tol:
                
                dupes += 1
                s += ("\n\nRows " + str(i+2) + " and " + str(j+2) +
                    "\nCoordinates closer than " +
                    f"{dist:f} mi ({5280*dist:.2f} ft)")
                
                # Gather addresses
                a1 = (dic[i][namefield] + ", " + dic[i][a1field] + ", " +
                      dic[i][a2field] + ", " + dic[i][cityfield] + ", " +
                      dic[i][statefield] + ", " + dic[i][zipfield])
                a2 = (dic[j][namefield] + ", " + dic[j][a1field] + ", " +
                      dic[j][a2field] + ", " + dic[j][cityfield] + ", " +
                      dic[j][statefield] + ", " + dic[j][zipfield])
                s += "\n" + a1 + "\n" + a2
    
    print(str(missing) + " missing coordinates found.")
    print(str(dupes) + " likely duplicates found.")
    
    # Write report file or print report to screen
    if outfile == None:
        print("\n\n" + s + "\n")
    else:
        with open(outfile, 'w') as f:
            f.write(s + "\n")

#------------------------------------------------------------------------------

def county_tract_info(county, cenfile, append=[], geocode=9, countycode=14,
                      basename=86, population=90, lat=92, lon=93, lsadc=94,
                      sumlev=(2,"140")):
    """Extracts tract-level population info from a Census Redistricting file.
    
    Positional arguments:
        counties (str|list(str)) -- County name or list of county names.
        cenfile (str) -- Path to a Census Redistricting (PL) file for the state
            containing the specified counties.
    
    Keyword arguments:
        append (list) -- List to append to the end of each dictionary entry.
            Defaults to the empty list. Meant to leave room for adding extra
            entries later.
        geocode (int) -- Column number (starting from 0) of the row's FIPS
            code. Defaults to 9.
        countycode (int) -- Column number of county-level FIPS code. Defaults
            to 14.
        basename (int) -- Column number of the record's name. Defaults to 86.
        population (int) -- Column number of the record's population number.
            Defaults to 90.
        lat (int) -- Column number of the record's latitude. Defaults to 92.
        lon (int) -- Column number of the record's longitude. Defaults to 93.
        lsadc (int) -- Column number of the record's Legal/Statistical Area
            Description Code. Defaults to 94.
        sumlev (tuple) -- Tuple indicating the column number of the
            summary level code (default 2) and the summary level corresponding
            to a census tract (default 140).
    
    Returns:
        (dict) -- Dictionary of tract-level information for all of the
            specified counties. The dictionary is indexed by the tract's FIPS
            code, and each entry is a list containing, respectively, the
            tract's latitude, longitude, and population (followed by the
            contents of the "append" argument)
    
    This function extracts tract-level population information from a US Census
    National Redistricting Data Summary File, which is a pipe-delimited table
    that contains a large amount of data for various levels of statistical
    region.
    
    The keyword arguments indicate the column numbers of the fields to extract,
    and are based on the 2020 summary file format. The lsadc field is used to
    determine which rows correspond to counties (code 06) and which correspond
    to census tracts (code CT).
    
    Counties are processed one-by-one. For each county, we begin by finding
    that county's row by finding a match in the basename field, and then we
    extract its county FIPS code from the countycode field. We then go through
    each census tract with a matching county code and extract the information
    from the lat, lon, and population fields.
    
    County names that are not found in the census file are skipped.
    """
    
    # Ensure that we have a list of county names
    if type(county) is not list and type(county) is not tuple:
        county = [county]
    
    # Initialize output dictionary
    pdic = dict()
    
    # Read census file into a table
    with open(cenfile, 'r') as f:
        tab = [line.strip().split('|') for line in f]
    
    # Process each county one-at-a-time
    for c in county:
        
        # Find the county's FIPS code
        cfips = ""
        for row in tab:
            if row[lsadc] == "06" and row[basename] == c:
                cfips = row[countycode]
                break
        
        # Skip the county if no code was found
        if cfips == "":
            continue
        
        # Process each tract that matches the county's FIPS code
        for row in tab:
            if (row[lsadc] != "CT" or row[countycode] != cfips or
                row[sumlev[0]] != sumlev[1]):
                continue
            
            # Store data in output dictionary
            pdic[row[geocode]] = ([float(row[lat]), float(row[lon]),
                                   int(row[population])] + append)
    
    # Return the data dictionary
    return pdic

#------------------------------------------------------------------------------

def ini_section_keys(inifile, section):
    """Reads all keys from a given INI section.
    
    Positional arguments:
        inifile (str) -- INI file path.
        section (str|list(str)) -- Section or list of sections.
    
    Returns:
        (list) -- List of all keys in the given INI section(s).
    """
    
    # Convert singleton section into a list
    if type(section) is not list and type(section) is not tuple:
        section = [section]
    
    # Initialize config file parser
    parser = configparser.ConfigParser(allow_no_value=True)
    parser.read(inifile)
    
    # Gather keys from all sections
    keylist = []
    for s in section:
        keylist.extend([k for k in parser[s]])
    
    return keylist

#------------------------------------------------------------------------------

def _pharmacy_type(name):
    """Attempts to classify the type of pharmacy from its name.
    
    Positional arguments:
        name (str) -- Name of the pharmacy.
    
    Returns:
        (str) -- Type of the pharmacy, or empty string if unclassified.
    """
    
    if "Costco".lower() in name.lower():
        return "Costco"
    elif (("CVS".lower() in name.lower()) or
          ("MinuteClinic".lower() in name.lower())):
        return "CVS"
    elif "Kroger".lower() in name.lower():
        return "Kroger"
    elif "Publix".lower() in name.lower():
        return "Publix"
    elif "Rite Aid".lower() in name.lower():
        return "Rite Aid"
    elif "Sam's Club".lower() in name.lower():
        return "Sam's Club"
    elif "Safeway".lower() in name.lower():
        return "Safeway"
    elif "Walgreens".lower() in name.lower():
        return "Walgreens"
    else:
        return ""

#------------------------------------------------------------------------------

def filter_providers(ziplist, filterfile):
    """Creates a filtered list of vaccine providers for a set of ZIP codes.
    
    Positional arguments:
        ziplist (list(int)) -- List of ZIP codes. Only facilities matching one
            of these ZIP codes will be included in the output file.
        filterfile (str) -- File path for the filtered copy of the vaccine
            provider file.
    
    This script is made for processing the vaccines.gov master provider list to
    produce filtered lists of providers within a given location. Keep in mind
    that the output file may still require hand-processing due to the
    possibility of duplicate entries, which is what the above address_test()
    script is for.
    """
    
    # Define master provider file name and field indices
    master_file = os.path.join("..", "data", "_general",
                   "Vaccines.gov__COVID-19_vaccinating_provider_locations.csv")
    
    # Read contents of master CSV file
    with open(master_file, 'r') as f:
        reader = list(csv.DictReader(f, delimiter=',', quotechar='"'))
    
    # Write output CSV file
    with open(filterfile, 'w', newline='') as f:
        
        # Define field names
        fields = ["pharmacy name", "pharmacy type", "address line 1",
                  "address line 2", "city", "state", "zipcode", "min age",
                  "latitude", "longitude"]
        
        # Initialize CSV parser
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        
        # Initialize name list
        names = []
        
        # Process each row of master file
        for row in tqdm.tqdm(reader):
            
            # Skip rows that don't match the ZIP code list
            if row["loc_admin_zip"][:5] not in ziplist:
                continue
            
            # Skip rows whose names match a previously-logged name
            if row["loc_name"] in names:
                continue
            
            # Record the new pharmacy name
            names.append(row["loc_name"])
            
            # Attempt to classify the location type
            type = _pharmacy_type(row["loc_name"])
            
            # Write the new row to the output file
            writer.writerow({"pharmacy name": row["loc_name"],
                             "pharmacy type": type,
                             "address line 1": row["loc_admin_street1"],
                             "address line 2": row["loc_admin_street2"],
                             "city": row["loc_admin_city"],
                             "state": row["loc_admin_state"],
                             "zipcode": row["loc_admin_zip"][:5],
                             "min age": row["min_age_years"],
                             "latitude": row["latitude"],
                             "longitude": row["longitude"]})

#------------------------------------------------------------------------------

def pharmacy_to_facility(pharmfile, facfile, offset=0):
    """Converts a pharmacy data file into a standardized facility file.
    
    Positional arguments:
        pharmfile (str) -- Path to pharmacy data file (the type output by the
            provider filtration script above).
        facfile (str) -- Path to the output facility file.
    
    Keyword arguments:
        offset (int) -- Offset for the first row's index. Defaults to 0. It may
            be desirable to set it to a different number if generating several
            facility files in sequence.
    
    Returns:
        (int) -- The index of the last row.
    
    The input and output files are both expected to follow a specific format.
    In particular the input pharmacy file is meant to be a .tsv following the
    format produced by the facility filtration script above, with the following
    fields (in order, with these names):
        pharmacy name
        pharmacy type
        address line 1
        address line 2
        city
        state
        zipcode
        min age
        latitude
        longitude
    
    The output facility file (for use in the metric generation scripts) is a
    simplified version including the following fields (in order, with these
    names):
        id (unique integer index)
        name (full name from the "pharmacy name" field)
        lat (latitude)
        lon (longitude)
        cap (capacity)
    """
    
    # Initialize facility dictionary
    fdic = dict()
    
    # Gather vaccination facility locations
    with open(pharmfile, 'r') as f:
        
        # Create a list of dictionaries for each row
        dic = list(csv.DictReader(f, delimiter=',', quotechar='"'))
        
        # Look up the address on each row
        for row in dic:
            
            # Get facility name (removing tabs if needed)
            fi = row["pharmacy name"].replace('\t', '')
            
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
    index = offset
    with open(facfile, 'w') as f:
        f.write(FAC_HEADER)
        sk = sorted(fdic.keys())
        for i in range(len(sk)):
            line = str(index) + '\t' + str(sk[i]) + '\t'
            for item in fdic[sk[i]]:
                line += str(item) + '\t'
            f.write(line + '\n')
            index += 1
    
    del fdic
    
    return index

#------------------------------------------------------------------------------

def county_svi(fips_prefix, svifile=os.path.join("..", "data", "_general",
               "SVI2020_US.csv"), fields=["RPL_THEMES"]):
    """Reads fields from the master SVI data file.
    
    Positional arguments:
        fips_prefix (str) -- FIPS code prefix. Every tract with this prefix
            will be logged in the output dictionary. For reference, the first
            2 digits correspond to the state while the next 3 correspond to the
            county.
    
    Keyword arguments:
        svifile (str) -- Path to SVI file. Defaults to this repo's
            data/_general/SVI2020_US.csv file.
        fields (list(str)) -- Field names to read. Defaults to ["RPL_THEMES"],
            which is the overall SVI percentile.
    
    Returns:
        (dict(tuple(float))) -- Dictionary of tract-level entries under the
            "field" columns for every census tract that matches the given FIPS
            prefix.
    """
    
    # Convert singleton fields argument to a list
    if type(fields) is not list and type(fields) is not tuple:
        fields = [fields]
    
    # Initialize output dictionary
    entries = dict()
    
    # Read contents of input file as a CSV
    with open(svifile, 'r') as f:
        reader = csv.DictReader(f, delimiter=',', quotechar='"')
        for row in reader:
            # Filter rows with the correct FIPS prefix
            if row["FIPS"][:len(fips_prefix)] != fips_prefix:
                continue
            # Gather all requested fields in the output dictionary
            entries[row["FIPS"]] = tuple(float(row[s]) for s in fields)
    
    return entries

#==============================================================================
# Location-Specific Preprocessing Scripts
#==============================================================================

def process_santa_clara(popfile=os.path.join("..", "processed", "santa_clara",
                                             "santa_clara_pop.tsv"),
                        nbrpopfile=os.path.join("..", "processed", "santa_clara",
                                                "santa_clara_pop_nbr.tsv"),
                        facfile=os.path.join("..", "processed", "santa_clara",
                                             "santa_clara_fac.tsv"),
                        nbrfacfile=os.path.join("..", "processed", "santa_clara",
                                                "santa_clara_fac_nbr.tsv")):
    """Preprocessing scripts for the Santa Clara data.
    
    Keyword arguments:
        popfile (str) -- Population output file path. Defaults to a file in the
            processed/ directory named "santa_clara_pop.tsv".
        nbrpopfile (str) -- Population neighbor output file path. Defaults to a
            file in the processed/ directory named "santa_clara_pop_nbr.tsv".
        facfile (str) -- Facility output file path. Defaults to a file in the
            processed/ directory named "santa_clara_fac.tsv".
        nbrfacfile (str) -- Facility neighbor output file path. Defaults to a
            file in the processed/ directory named "santa_clara_fac_nbr.tsv".
    """
    
    # Define location-specific file names
    adi_file = os.path.join("..", "data", "santa_clara",
                            "CA_2020_ADI_Census Block Group_v3.2.csv")
    svi_file = os.path.join("..", "data", "_general", "SVI2020_US.csv")
    fac_file = os.path.join("..", "data", "santa_clara",
                            "Santa_Clara_County_Pharmacies.csv")
    fac_nbr_file = os.path.join("..", "data", "santa_clara",
                                "Santa_Clara_County_Neighbor_Pharmacies.csv")
    census_file = os.path.join("..", "data", "ca", "cageo2020.pl")
    vacc_file = os.path.join("..", "data", "santa_clara",
             "COVID-19_Vaccination_among_County_Residents_by_Census_Tract.csv")
    urban_file = os.path.join("..", "data", "santa_clara",
        "santa_clara_urban.csv")
    
    # Define location-specific parameters
    neighbors = SANTA_CLARA_NEIGHBORS
    
    # Gather population data from census file (leaving room for extra fields)
    pdic = county_tract_info("Santa Clara", census_file, append=[-1]*6)
    
    # Gather vaccination rates
    with open(vacc_file, 'r') as f:
        
        for line in f:
            
            # Skip comment line
            if line[0].isdigit() == False:
                continue
            
            s = line.strip().split(',')
            fips = s[0] # current row's FIPS
            if fips not in pdic:
                continue
            
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
            fips = s[3][:11]
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
    
    # Gather SVI fields
    santa_clara_fips = list(pdic.keys())[0][:5] # county's FIPS prefix
    fields = ["RPL_THEMES", "EP_POV150", "EP_NOVEH"] # field names to gather
    fieldids = [5, 7, 8] # corresponding columns in final table
    fieldpercent = [7, 8] # fields to be converted from [0,100] to [0.0,1.0]
    svi = county_svi(santa_clara_fips, fields=fields) # SVI dictionary
    for fips in pdic:
        if fips in svi:
            for i in range(len(svi[fips])):
                pdic[fips][fieldids[i]] = svi[fips][i]
                if fieldids[i] in fieldpercent:
                    pdic[fips][fieldids[i]] /= 100.0
    
    del svi
    
    # Gather urban/rural fractions
    with open(urban_file, 'r') as f:
        reader = csv.DictReader(f, delimiter=',', quotechar='"')
        for row in reader:
            # Find any and all FIPS codes that match the current row
            prefix = row["GEOID"]
            flist = []
            for fips in pdic:
                if fips[:len(prefix)] == prefix:
                    flist.append(fips)
            # Skip FIPS codes with no matches
            if len(flist) < 1:
                continue
            # Assign urban fraction to all collected FIPS code
            for fips in flist:
                pdic[fips][6] = row["urbanPercent"]
    
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
            f.write(line[:-1] + '\n')
            index += 1
    
    # Gather neighboring county data
    pndic = county_tract_info(neighbors, census_file)
    
    # Delete duplicates
    for fips in pdic:
        if fips in pndic:
            del pndic[fips]
    
    # Write population neighbor output file
    with open(nbrpopfile, 'w') as f:
        f.write(POP_HEADER_SHORT)
        sk = sorted(pndic.keys())
        # "index" carries over from the main population file
        for i in range(len(sk)):
            # Skip lines with no coordinates
            if pndic[sk[i]][0] == 0 or pndic[sk[i]][1] == 0:
                continue
            line = str(index) + '\t' + str(sk[i]) + '\t'
            for item in pndic[sk[i]]:
                line += str(item) + '\t'
            f.write(line[:-1] + '\n')
            index += 1
    
    del pdic
    del pndic
    
    # Generate facility file
    index = pharmacy_to_facility(fac_file, facfile)
    
    # Generate facility neighbor file
    pharmacy_to_facility(fac_nbr_file, nbrfacfile, offset=index)

#==============================================================================
# Execution
#==============================================================================

# Comment or uncomment the function calls below to process each location.

#zips = ini_section_keys(os.path.join("..", "data", "santa_clara", "santa_clara_zips.ini"), "santa_clara")
#filter_providers(zips, os.path.join("..", "data", "santa_clara", "Santa_Clara_County_Pharmacies.csv"))
#address_test(os.path.join("..", "data", "santa_clara", "Santa_Clara_County_Pharmacies.csv"), outfile=os.path.join("..", "data", "santa_clara", "Report.txt"))

#zips = ini_section_keys(os.path.join("..", "data", "santa_clara", "santa_clara_zips.ini"), [s.lower().replace(' ','_') for s in SANTA_CLARA_NEIGHBORS])
#filter_providers(zips, os.path.join("..", "data", "santa_clara", "Santa_Clara_County_Neighbor_Pharmacies.csv"))
#address_test(os.path.join("..", "data", "santa_clara", "Santa_Clara_County_Neighbor_Pharmacies.csv"), outfile=os.path.join("..", "data", "santa_clara", "Report.txt"))

#county_tract_info("Santa Clara", os.path.join("..", "data", "ca", "cageo2020.pl"))
process_santa_clara(popfile=os.path.join("..", "processed", "santa_clara", "santa_clara_pop.tsv"), facfile=os.path.join("..", "processed", "santa_clara", "santa_clara_fac.tsv"))
