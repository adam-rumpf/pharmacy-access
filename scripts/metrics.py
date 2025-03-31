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

def _read_schedfile(schedfile, num=True):
    """Reads a schedule file and returns dictionaries of data.
    
    Positional arguments:
        schedfile (str) -- Preprocessed schedule file path, which should
            include columns indicating binary indicators of whether a facility
            is open during each time slot.
    
    Optional keyword arguments:
        num (bool) -- Whether to convert time slot titles to numerical indices.
            Defaults to True, in which case dictionary keys are integer
            indices that correspond to time slot names. If False, the time slot
            names, themselves, are the dictionary keys.
    
    Returns:
        sched (dict(dict(bool))) -- Dictionary of dictionaries of time slot
            availability indicators (True or False).
         names (list(str)) -- List of time slot names.
    
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
                for i in range(len(slots)):
                    # Determine key type
                    k = i
                    if num == False:
                        k = slots[i]
                    sched[k] = dict()
                continue
            
            # Read the availability numbers for the current facility
            s = line.strip().split('\t')
            id = int(s[0])
            for i in range(len(slots)):
                # Determine key type
                k = i
                if num == False:
                    k = slots[i]
                sched[k][id] = bool(int(s[i+SCHED_OFFSET]))
    
    return sched, slots

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

def _is_open(sdic, fid, tstart, tfinish=None):
    """Determines whether a facility is open during a given time slot or range.
    
    Positional arguments:
        sdic (dict(dict)) -- Schedule dictionary of dictionaries, indexed first
            by facility ID and second by time ID. Defaults to None, in which
            case hours are ignored.
        fid (int) -- Facility index.
        tstart (int) -- Index of time slot, or index of first time slot in a
            range.
    
    Optional keyword argument:
        tfinish (int) -- Index of final time slot in a range (inclusive).
            Defaults to None, in which case only the first time slot is
            considered.
    
    Returns:
        (bool) -- True if the facility is open during the given time slot, or
            at any point during the given range of time slots.
    """
    
    # Set default finish time slot
    if tfinish == None:
        tfinish = tstart
    
    # Loop through time slots to find one where the facility is open
    for t in range(tstart, tfinish+1):
        if sdic[t][fid] == True:
            return True
    return False

#------------------------------------------------------------------------------

def _count_open(sdic, fid, tstart, tfinish=None):
    """Counts number of time slots in a range during which a facility is open.
    
    Positional arguments:
        sdic (dict(dict)) -- Schedule dictionary of dictionaries, indexed first
            by facility ID and second by time ID. Defaults to None, in which
            case hours are ignored.
        fid (int) -- Facility index.
        tstart (int) -- Index of time slot, or index of first time slot in a
            range.
    
    Optional keyword argument:
        tfinish (int) -- Index of final time slot in a range (inclusive).
            Defaults to None, in which case only the first time slot is
            considered.
    
    Returns:
        (int) -- The number of time slots during the specified range during
            which the facility is open.
    """
    
    # Set default finish time slot
    if tfinish == None:
        tfinish = tstart
    
    # Loop through time slots and count number during which facility is open
    total = 0
    for t in range(tstart, tfinish+1):
        if sdic[t][fid] == True:
            total += 1
    return total

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

def cutoff_count(pdic, fdic, ddic, cutoff, sdic=None, hours=None, avg=False):
    """Counts facilities within a given travel time cutoff.
    
    Positional arguments:
        pdic (dict) -- Population dictionary indexed by population center ID.
        fdic (dict) -- Facility dictionary indexed by facility ID.
        ddic (dict(dict)) -- Distance dictionary of dictionaries, indexed first
            by population center ID and second by facility ID.
        cutoff (float) -- Travel time cutoff.
    
    Optional keyword arguments:
        sdic (dict(dict)) -- Schedule dictionary of dictionaries, indexed first
            by facility ID and second by time ID. Defaults to None, in which
            case hours are ignored.
        hours (None|tuple(int,int)) -- Whether to filter to a specific time
            window. Defaults to None, in which case hours are ignored.
            Otherwise should be a tuple containing the first and last schedule
            indices (inclusive).
        avg (bool) -- Whether to compute an average facility count over the
            given time period. Defaults to False, in which case the result is
            a simple count of facilities available at any point within the
            time window.
    
    Returns:
        (dict(int)) -- Dictionary of facility counts, indexed by population
            center ID.
    """
    
    # Initialize dictionary
    counts = dict([(pid, 0) for pid in pdic])
    
    # Shchedule-free case
    if sdic == None or hours == None:
        # Go through each population/facility distance pair
        for pid in pdic:
            for fid in fdic:
                # Increment facility count if below distance cutoff
                if ddic[pid][fid] <= cutoff:
                    counts[pid] += 1
    
    # Scheduled count
    elif avg == False:
        # Go through each population/facility distance pair
        for pid in pdic:
            for fid in fdic:
                # Skip if facility is closed during all time slots
                if _is_open(sdic, fid, hours[0], hours[1]) == False:
                    continue
                
                # Otherwise, increment facility count if below distance cutoff
                if ddic[pid][fid] <= cutoff:
                    counts[pid] += 1
    
    # Scheduled average
    else:
        # Go through each population/facility distance pair
        for pid in pdic:
            total = 0
            for fid in fdic:
                # Skip facilities out of range
                if ddic[pid][fid] > cutoff:
                    continue
                
                # Compute number of time slots during which facility is open
                total += _count_open(sdic, fid, hours[0], hours[1])
            
            # Compute average
            counts[pid] = total/(hours[1] - hours[0] + 1)
    
    return counts

#------------------------------------------------------------------------------

def fraction_above_threshold(pdic, fdic, ddic, cutoff, sdic, hours, count):
    """Computes fraction of time period above a given facility count.
    
    Positional arguments:
        pdic (dict) -- Population dictionary indexed by population center ID.
        fdic (dict) -- Facility dictionary indexed by facility ID.
        ddic (dict(dict)) -- Distance dictionary of dictionaries, indexed first
            by population center ID and second by facility ID.
        cutoff (float) -- Travel time cutoff.
        sdic (dict(dict)) -- Schedule dictionary of dictionaries, indexed first
            by facility ID and second by time ID.
        hours (tuple(int,int)) -- Tuple containing the first and last schedule
            indices (inclusive).
        count (int) -- Facility count to use as the lower threshold.
    
    Returns:
        (dict(float)) -- Dictionary of fractions, indexed by population center
            ID. The number represents the fraction of the time period
            indicated by "hours" during which the population center has
            access to at least "count" facilities within travel time "cutoff".
    """
    
    # Generate a list of count dictionaries for each time slot in the range
    clist = [cutoff_count(pdic, fdic, ddic, cutoff, sdic, (h, h))
             for h in range(hours[0], hours[1]+1)]
    
    # For each population center, count times slots above the facility threshold
    frac = dict() # dictionary of fractions
    for pid in pdic:
        total = 0 # number of time slots meeting threshold
        for c in clist:
            if c[pid] >= count:
                total += 1
        frac[pid] = total/len(clist)
    
    return frac

#------------------------------------------------------------------------------

def avg_time_to_nearest(pdic, fdic, ddic, counts, sdic, hours):
    """Computes average travel time to k nearest facilities.
    
    Positional arguments:
        pdic (dict) -- Population dictionary indexed by population center ID.
        fdic (dict) -- Facility dictionary indexed by facility ID.
        ddic (dict(dict)) -- Distance dictionary of dictionaries, indexed first
            by population center ID and second by facility ID.
        counts (tuple(int)) -- Tuple of facility counts to include in the
            averaging.
        sdic (dict(dict)) -- Schedule dictionary of dictionaries, indexed first
            by facility ID and second by time ID.
        hours (tuple(int,int)) -- Tuple containing the first and last schedule
            indices (inclusive).
    
    Returns:
        TBD###
    
    Due to the computational overhead required to compute distances to k-nearest
    facilities, this function is made to compute a collection of metrics all at
    once from the same lists of sorted travel times.
    """
    
    pass###

#==============================================================================
# Metric Compilation Scripts
#==============================================================================

def all_metrics(pinfile, poutfile, facfile, distfile, schedfile, cutoffs, hours,
                facnums):
    """Driver to generate all population metrics.
    
    Positional arguments:
        pinfile (str) -- Path to input population file.
        poutfile (str) -- Path to output population file.
        facfile (str) -- Path to facility file.
        distfile (str) -- Path to input distance file.
        schedfile (str) -- Path to input schedule file.
        cutoffs (tuple(float)) -- Tuple of travel time cutoffs to include.
        hours (tuple(int)) -- Tuple of starting and ending schedule file indices
            to define a limited time window (inclusive).
        facnums (tuple(int)) -- Tuple of numbers of facilities to average over.
            These are the k's for which some of the metrics measure things like
            average travel times to the k nearest facilities.
    
    This function runs through a set of the above metric generation functions.
    All results are appended to the given population file as new fields.
    """
    
    # Read dictionaries from input files
    pdic, _, _, _ = _read_popfile(pinfile)
    fdic, _ = _read_facfile(facfile)
    ddic = _read_distfile(distfile)
    sdic, hournames = _read_schedfile(schedfile)
    
    # Compute total population
    totpop = 0
    for pid in pdic:
        totpop += pdic[pid]
    print(f"Total population: {totpop}")
    
    ## Get time increment from schedule file
    #dt = 60 # number of minutes between schedule slots
    #with open(schedfile, 'r') as f:
    #    # Get first two time strings
    #    line = f.readline().split('\t')
    #    t1 = int(line[SCHED_OFFSET][-2:])
    #    t2 = int(line[SCHED_OFFSET+1][-2:])
    #    if t2 - t1 > 0:
    #        dt = t2 - t1
    
    # Get start and end time names
    tstart = hournames[hours[0]]
    tfinish = hournames[hours[1]+1]
    
    # Initialize list of result dictionaries
    labels = [] # labels for new columns
    results = [] # dictionaries containing results for each column
    
    # Generate counts, both with and without store hours
    for t0 in cutoffs:
        
        # Counts and averages
        labels.append(f"fac-count_all-times_cutoff-{t0}")
        results.append(cutoff_count(pdic, fdic, ddic, t0))
        countcol = len(results) - 1 # all-times count column for thresholds
        labels.append(f"fac-count_{tstart}-{tfinish}_cutoff-{t0}")
        results.append(cutoff_count(pdic, fdic, ddic, t0, sdic, hours))
        labels.append(f"fac-count-avg_{tstart}-{tfinish}_cutoff-{t0}")
        results.append(cutoff_count(pdic, fdic, ddic, t0, sdic, hours, True))
        
        # Indicators for facility count thresholds
        for fn in facnums:
            # Initialize total population
            total = 0
            # Generate result column
            labels.append(f"fac-count-below-{fn}_all-times_cutoff-{t0}")
            results.append(dict())
            for pid in pdic:
                results[-1][pid] = 0 # indicator for below threshold
                if results[countcol][pid] < fn:
                    results[-1][pid] = 1
                    total += pdic[pid]
            # Print global statistic
            print(f"Total population lacking access to {fn} facilities " +
                  f"within {t0} minutes: {total}; Fraction: {total/totpop:f}")
            # Fraction of time slot during which threshold is met
            labels.append(f"frac-time-above-{fn}_{tstart}-{tfinish}_cutoff-{t0}")
            results.append(fraction_above_threshold(pdic, fdic, ddic, t0, sdic,
                                                    hours, fn))
    
    ###
    
    # Write results
    _augment_file(pinfile, poutfile, results[0], labels[0])
    for i in range(1, len(labels)):
        _augment_file(poutfile, poutfile, results[i], labels[i])

#==============================================================================

#print(_read_distfile(os.path.join(POLK_PROCESSED, "polk_dist_pharmacy.tsv")))

# Polk file paths
polk_pop = os.path.join(POLK_PROCESSED, "polk_pop.tsv")
polk_pharm = os.path.join(POLK_PROCESSED, "polk_pharmacy_all.tsv")
polk_uc = os.path.join(POLK_PROCESSED, "polk_uc_all.tsv")
polk_sched_pharm = os.path.join(POLK_PROCESSED, "polk_schedule_pharmacy_15.tsv")
polk_sched_uc = os.path.join(POLK_PROCESSED, "polk_schedule_uc_15.tsv")
polk_dist_pharm = os.path.join(POLK_PROCESSED, "polk_dist_pharmacy.tsv")
polk_dist_uc = os.path.join(POLK_PROCESSED, "polk_dist_uc.tsv")
polk_pop_pharm_results = os.path.join(POLK_RESULTS, "polk_pop_pharm_results.tsv")
polk_pop_uc_results = os.path.join(POLK_RESULTS, "polk_pop_uc_results.tsv")

# Additional parameters
cutoffs = (15, 30) # travel time cutoffs
hours = (260, 279) # WED 5:00pm-5:15pm through WED 9:45pm-10:00pm
facnums = (1, 3, 5) # numbers of nearest facilities to average over

print("="*20 + "\nPharmacy tests\n" + "="*20)
all_metrics(polk_pop, polk_pop_pharm_results, polk_pharm, polk_dist_pharm, polk_sched_pharm, cutoffs, hours, facnums)
print("\n" + "="*20 + "\nUrgent care tests\n" + "="*20)
all_metrics(polk_pop, polk_pop_uc_results, polk_uc, polk_dist_uc, polk_sched_uc, cutoffs, hours, facnums)
