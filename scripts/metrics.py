"""Pharmacy accessibility project accessibility metrification scripts.

The code below defines various functions that work with the preprocessed data
files, all of which should have a standardized format. Various accessibility
metrics can be computed, and various distance measurements can be used.
"""

import os.path

import tqdm

#==============================================================================
# Global Constants
#==============================================================================

# Preprocessed file column numbers
FAC_LAT = 2 # facility latitude
FAC_LON = 3 # facility longitude
FAC_CAP = 4 # facility capacity
POP_LAT = 2 # population latitude
POP_LON = 3 # population longitude
POP_POP = 4 # population total population value
POP_SVI = 7 # population SVI ranking
POP_URBAN = 8 # population urban fraction
SCHED_OFFSET = 2 # first column of the schedule file to contain time info

# Processed file roots
POLK_PROCESSED = os.path.join("..", "processed", "polk")

# Result file roots
POLK_RESULTS = os.path.join("..", "results", "polk")

#==============================================================================
# Common Functions
#==============================================================================

def _read_popfile(popfile):
    """Reads a population file and returns dictionaries of data.
    
    Positional arguments:
        popfile (str) -- Preprocessed population file path, which should
            include the coordinates and population of each population
            center.
    
    Returns:
        pop (dict(int)) -- Dictionary of population counts.
        coord (dict((float,float))) -- Dictionary of population center
            coordinates, as (latitude,longitude) tuples.
        svi (dict(float)) -- Dictionary of SVI rankings.
        urban (dict(float)) -- Dictionary of urban fractions.
    
    All dictionaries are indexed by the population center IDs contained in the
    first column of the population file.
    """
    
    # Initialize dictionaries
    pop = dict() # dictionary of populations by popfile index
    coord = dict() # dictionary of latitutde/longitude pairs by popfile index
    svi = dict() # dictionary of SVI rankings
    urban = dict() # dictionary of urban/rural splits (fraction urban out of 1)
    
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
            svi[int(s[0])] = float(s[POP_SVI])
            urban[int(s[0])] = float(s[POP_URBAN])
    
    return (pop, coord, svi, urban)

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

def _read_schedfile(schedfile):
    """Reads a schedule file and returns dictionaries of data.
    
    Positional arguments:
        schedfile (str) -- Preprocessed schedule file path, which should
            include columns indicating binary indicators of whether a facility
            is open during each time slot.
    
    Returns:
        sched (dict(dict(bool))) -- Dictionary of dictionaries of time slot
            availability indicators (True or False).
    
    The returned dictionary contains one entry for every time slot specified
    in the schedule file, indexed by the name of the time slot. Each entry is,
    itself, a dictionary, indexed by facility ID and containing a True if the
    facility is open during that time slot and a False otherwise.
    """
    
    # Initialize dictionaries
    sched = dict() # dictionary of dictionaries of schedule items
    
    # Read file
    with open(schedfile, 'r') as f:
        
        first = True
        slots = []
        
        for line in f:
            
            # Gather the time slot names from the first line
            if first == True:
                first = False
                # Initialize a dictionary entry for each time slot
                slots = line.strip().split('\t')[SCHED_OFFSET:]
                for slot in slots:
                    sched[slot] = dict()
                continue
            
            # Read the availability numbers for the current facility
            s = line.strip().split('\t')
            id = int(s[0])
            for i in range(len(slots)):
                sched[slots[i]][id] = bool(s[i+SCHED_OFFSET])
    
    return sched

#------------------------------------------------------------------------------

def _read_distfile(distfile):
    """Reads a distance file and returns dictionaries of data.
    
    Positional arguments:
        distfile (str) -- Preprocessed distance file path, which should include
            a column of population center IDs followed by a column of facility
            IDs followed by a column of travel times.
    
    Returns:
        dist (dict(dict(float))) -- Dictionary of dictionaries of travel times.
    
    The returned dictionary contains one entry for every population center,
    indexed by center ID. Each of these entries is, itself, a dictionary,
    indexed by facility ID, containing the travel time from that center to that
    facility.
    """
    
    # Initialize dictionary
    dist = dict() # dictionary of dictionaries of travel times
    
    # Read file
    with open(distfile, 'r') as f:
        
        first = True
        
        for line in f:
            
            # Skip the first line
            if first == True:
                first = False
                continue
            
            # Get entries
            s = line.strip().split('\t')
            pid = int(s[0]) # population ID
            fid = int(s[1]) # facility ID
            t = float(s[2]) # travel time
            
            # Add new population entry if needed
            if pid not in dist:
                dist[pid] = dict()
            
            # Log travel time
            dist[pid][fid] = t
    
    return dist

#------------------------------------------------------------------------------

def _augment_file(infile, outfile, column, label, default="-1"):
    """Augments a population or facility file with a new column.
    
    Positional arguments:
        infile (str) -- Preprocessed population/facility file path.
        outfile (str) -- Output population/facility file path.
        column (dict(float)) -- Dictionary of entries for the new table
            column. The keys of the dictionary should correspond to the
            indices in the first column of the input file.
        label (str) -- Label for the new column. Since the output file is
            a tab-separated value table, avoid the use of tab characters.
    
    Keyword arguments:
        default (str) -- Default value for rows with no corresponding
            dictionary value. Defaults to "-1".
    
    This is meant for use in adding a column of population accessibility
    metrics to a population file, or facility crowding metrics to a facility
    file. It works by copying each row of the original file and adding a new
    element to the end of each row corresponding to that row's index.
    """
    
    # Append new columns to input file lines
    with open(infile, 'r') as f:
        lines = ""
        
        for line in f:
            
            # Copy header with new column label
            if line[0].isdigit() == False:
                lines += line.strip() + '\t' + label + '\n'
                continue
            
            # Append the new column element to the following rows
            i = int(line.split('\t')[0]) # current row index
            if i in column:
                c = column[i] # dictionary element of current line
            else:
                c = default # default element for missing indices
            lines += line.strip() + '\t' + str(c) + '\n'
    
    # Write output file
    with open(outfile, 'w') as f:
        for line in lines:
            f.write(line)

#==============================================================================
# Specific Metrics
#==============================================================================

### Common dictinaries:
# pop[i] - population of tract i
# svi[i] - SVI of tract i
# sched[k][j] - boolean if facility j is open during time slot k
# dist[i][j] - travel time from tract i to facility j



#==============================================================================
# Metric Compilation Scripts
#==============================================================================

def all_metrics(pinfile, poutfile, distfile, schedfile, tstart, tfinish):
    """Driver to generate all population metrics.
    
    Positional arguments:
        pinfile (str) -- Path to input population file.
        poutfile (str) -- Path to output population file.
        distfile (str) -- Path to input distance file.
        schedfile (str) -- Path to input schedule file.
        tstart (ind) -- Starting time index for computing metrics within a
            limited window.
        tfinish (ind) -- Ending time index for computing metrics within a
            limited window. This final time index *is* included.
    
    This function runs through a set of the above metric generation functions.
    All results are appended to the given population file as new fields.
    """
    
    # Read dictionaries from input files
    popdic, _, _, _ = _read_popfile(pinfile)
    distdic = _read_distfile(distfile)
    scheddic = _read_schedfile(schedfile)
    
    # Get time increment from schedule file
    dt = 60 # number of minutes between schedule slots
    with open(schedfile, 'r') as f:
        # Get first two time strings
        line = f.readline().split('\t')
        t1 = int(line[SCHED_OFFSET][-2:])
        t2 = int(line[SCHED_OFFSET+1][-2:])
        if t2 - t1 > 0:
            dt = t2 - t1
    
    #

#==============================================================================

#print(_read_distfile(os.path.join(POLK_PROCESSED, "polk_dist_pharmacy.tsv")))

# Polk file paths
polk_pop = os.path.join(POLK_PROCESSED, "polk_pop.tsv")
polk_sched_pharm = os.path.join(POLK_PROCESSED, "polk_schedule_pharmacy_15.tsv")
polk_sched_uc = os.path.join(POLK_PROCESSED, "polk_schedule_uc_15.tsv")
polk_dist_pharm = os.path.join(POLK_PROCESSED, "polk_dist_pharmacy.tsv")
polk_dist_uc = os.path.join(POLK_PROCESSED, "polk_dist_uc.tsv")
polk_pop_pharm_results = os.path.join(POLK_RESULTS, "polk_pop_pharm_results.tsv")
polk_pop_uc_results = os.path.join(POLK_RESULTS, "polk_pop_uc_results.tsv")

# Starting and ending time indices
tstart = 262 # WED 5:00pm-5:15pm
tfinish = 281 # WED 9:45pm-10:00pm

all_metrics(polk_pop, polk_pop_pharm_results, polk_dist_pharm, polk_sched_pharm, tstart, tfinish)
