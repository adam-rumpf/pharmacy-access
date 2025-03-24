"""Pharmacy accessibility project accessibility metrification scripts.

The code below defines various functions that work with the preprocessed data
files, all of which should have a standardized format. Various accessibility
metrics can be computed, and various distance measurements can be used.
"""

import csv
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
            include columns indicating the fraction of each time slot during
            which a facility is open.
    
    Returns:
        sched (dict(dict(float))) -- Dictionary of dictionaries of time slot
            availability numbers.
    
    The returned dictionary contains one entry for every time slot specified
    in the schedule file, indexed by the name of the time slot. Each entry is,
    itself, a dictionary, indexed by facility ID and containing the fraction of
    that time slot during which that facility is open.
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
                sched[slots[i]][id] = float(s[i+SCHED_OFFSET])
    
    return sched

#------------------------------------------------------------------------------

def _augment_file(outfile, infile, column, label, default="-1"):
    """Augments a population or facility file with a new column.
    
    Positional arguments:
        outfile (str) -- Output population/facility file path.
        infile (str) -- Preprocessed population/facility file path.
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

##

#==============================================================================
# Metric Compilation Scripts
#==============================================================================

##
