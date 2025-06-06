"""Pharmacy accessibility project accessibility metrification scripts.

This is the original (deprecated) set of scripts.

The code below defines various functions that work with the preprocessed data
files, all of which should have a standardized format. Various accessibility
metrics can be computed, and various distance measurements can be used.
"""

import os.path

import geopy.distance
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

#------------------------------------------------------------------------------

def _threshold(x, cutoff, a, b):
    """Returns a if x >= cutoff and b otherwise."""
    
    if x >= cutoff:
        return a
    else:
        return b

#------------------------------------------------------------------------------

def _convex(x, a, b):
    """Returns x*a + (1-x)*b."""
    
    return x*a + (1-x)*b

#==============================================================================
# Gravity Metric Scripts
#==============================================================================

def gravity_metric(poutfile, foutfile, popfile, facfile, distfile=None,
                   popnbrfile=None, facnbrfile=None, beta=1.0, crowding=True,
                   floor=0.0):
    """Computes a table of gravitational metrics for a given community.
    
    Deprecated for project revision (switching to 2SFCA metrics).
    
    Positional arguments:
        poutfile (str) -- Output population file path.
        foutfile (str) -- Output facility file path.
        popfile (str) -- Preprocessed population file path, which should
            include the coordinates and population of each population center.
        facfile (str) -- Preprocessed facility file path, which should include
            the coordinates and capacity of each vaccination facility.
    
    Keyword arguments:
        distfile (str) -- Preprocessed distance file path, which should include
            the travel times between each population center/facility pair.
            Defaults to None, in which case geodesic distances are computed and
            used as needed.
        popnbrfile (str) -- Preprocessed neighboring county population file
            path. Defaults to None.
        facnbrfile (str) -- Preprocessed neighboring county facility file path.
            Defaults to None.
        beta (float) -- Gravitational decay parameter. Defaults to 1.0.
        crowding (bool) -- Whether or not to take crowding into consideration.
            Defaults to True.
        floor (float) -- Travel time floor (minutes) to use in computing
            gravitational metrics. Defaults to 0.0. It may be desirable to
            increase the floor to a nonzero value to prevent excessively large
            metrics for extremely small distances.
    
    This function implements a gravitational accessibility metric as described
    in Luo and Wang 2003 (doi:10.1068/b29120). We begin by computing a
    distance-weighted "crowdedness" metric for each facility consisting of a
    sum of all populations, each scaled by a weight that decays over distance
        
        V_j = sum_(k=1)^m P_k d_kj^(-beta)
    
    where V_j is the crowdedness metric for facility j, P_k is the population
    of center k, d_kj is the distance (or travel time) from k to j, and beta
    is a parameter that controls how quickly the weight should decay over
    distance (larger beta causes more distant locations to receive less
    weight).
    
    We then compute the accessibility metric for each community as a
    distance-weighted sum of "capacity-to-crowdedness" ratios for each
    facility consisting of a sum of all capacity to crowdedness ratios over
    all facilities, each scaled by a weight that decays over distance
        
        A_i = sum_(j=1)^n (S_j d_ij^(-beta))/V_j
    
    where A_i is the metric of population center i, S_j is the capacity of
    facility j, and the remaining notation is the same as above.
    
    The formulas above describe the default behavior of this script. If the
    "crowding" flag is set to False, an alternate crowding-free metric will be
    computed for each population center by setting V_j = 1 for all facilities,
    resulting in accessibility metrics that represent distance-weighted counts
    of facility capacities.
    
    Optionally, population centers and vaccination facilities in neighboring
    counties can be provided, in which case both will be used in computing
    crowding and accessibility metrics. This can be used to combat edge
    effects.
    
    The output files are augmented versions of the input files, with a new
    column of crowdedness or accessibility metrics appended to the original
    tables.
    """
    
    # Call an appropriate subroutine depending on distance specification
    if distfile == None:
        (fmet, pmet) = _gravity_metric_geodesic(popfile, facfile, beta=beta,
                                                popnbrfile=popnbrfile,
                                                facnbrfile=facnbrfile,
                                                crowding=crowding, floor=floor)
    else:
        (fmet, pmet) = _gravity_metric_file(popfile, facfile, distfile,
                                            beta=beta, popnbrfile=popnbrfile,
                                            facnbrfile=facnbrfile,
                                            crowding=crowding, floor=floor)
    
    # Write facility output file with a new crowdedness metric column
    if crowding == True:
        print("Writing facility metric file.")
        _augment_file(foutfile, facfile, fmet, "crowding")
    
    # Write population output file with a new accessibility metric column
    print("Writing population metric file.")
    _augment_file(poutfile, popfile, pmet, "access")

#------------------------------------------------------------------------------

def _gravity_metric_geodesic(popfile, facfile, beta=1.0, popnbrfile=None,
                             facnbrfile=None, crowding=True, floor=0.0,
                             speed=45.0):
    """The geodesic distance version of gravity_metric.
    
    Deprecated for project revision (switching to 2SFCA metrics).
    
    This script computes gravity metrics using geodesic distances computed as-
    needed. It is automatically called by gravity_metric when a distance file
    is not provided. It returns dictionaries of facility and population metrics
    rather than writing the results directly to a file.
    
    Geodesic distances are used to compute travel times in minutes, assuming a
    constant speed (in mph).
    """
    
    # Convert speed to miles per minute
    speed /= 60.0
    
    # Read population file
    (pop, pcoord) = _read_popfile(popfile)
    pkeys = list(pop.keys()) # store main population keys
    
    # If a neighboring population file is provided, merge its contents
    if popnbrfile != None:
        (popnbr, pcoordnbr) = _read_popfile(popnbrfile)
        pop = {**pop, **popnbr}
        pcoord = {**pcoord, **pcoordnbr}
        del popnbr
        del pcoordnbr
    
    # Read facility file
    (cap, fcoord) = _read_facfile(facfile)
    fkeys = list(cap.keys()) # store main facility keys
    
    # If a neighboring facility file is provided, merge its contents
    if facnbrfile != None:
        (capnbr, fcoordnbr) = _read_popfile(facnbrfile)
        cap = {**cap, **capnbr}
        fcoord = {**fcoord, **fcoordnbr}
        del capnbr
        del fcoordnbr
    
    # Compute the facility crowdedness metrics
    fmet = dict() # facility crowdedness metrics by facfile index
    if crowding == True:
        print("Computing facility crowdedness metrics using geodesic distance.")
        for j in tqdm.tqdm(cap):
            fmet[j] = 0.0
            for k in pop:
                # Compute travel time (minutes) and bound below by time floor
                d = max(geodesic_distance(pcoord[k], fcoord[j])/speed, floor)
                fmet[j] += pop[k]*(d**(-beta))
    else:
        for j in cap:
            fmet[j] = 1.0
    
    # Compute the population accessibility metrics
    pmet = dict() # population accessibility metrics by popfile index
    print("Computing population accessibility metrics using geodesic distance.")
    for i in tqdm.tqdm(pop):
        pmet[i] = 0.0
        for j in cap:
                # Compute travel time (minutes) and bound below by time floor
            d = max(geodesic_distance(pcoord[i], fcoord[j])/speed, floor)
            pmet[i] += (cap[j]*(d**(-beta)))/fmet[j]
    
    # Return facility and population metric dictionaries (original indices only)
    return ({k: fmet[k] for k in fkeys}, {k: pmet[k] for k in pkeys})

#------------------------------------------------------------------------------

def _gravity_metric_file(popfile, facfile, distfile, beta=1.0, popnbrfile=None,
                         facnbrfile=None, crowding=True, floor=0.0):
    """The distance file version of gravity_metric.
    
    Deprecated for project revision (switching to 2SFCA metrics).
    
    This script computes gravity metrics using a predefined distance file. It
    is automatically called by gravity_metric when a distance file is provided.
    It returns dictionaries of facility and population metrics rather than
    writing the results directly to a file.
    """
    
    # Read population file
    (pop, pcoord) = _read_popfile(popfile)
    pkeys = list(pop.keys()) # store main population keys
    
    # If a neighboring population file is provided, merge its contents
    if popnbrfile != None:
        (popnbr, pcoordnbr) = _read_popfile(popnbrfile)
        pop = {**pop, **popnbr}
        pcoord = {**pcoord, **pcoordnbr}
        del popnbr
        del pcoordnbr
    
    # Read facility file
    (cap, fcoord) = _read_facfile(facfile)
    fkeys = list(cap.keys()) # store main facility keys
    
    # If a neighboring facility file is provided, merge its contents
    if facnbrfile != None:
        (capnbr, fcoordnbr) = _read_popfile(facnbrfile)
        cap = {**cap, **capnbr}
        fcoord = {**fcoord, **fcoordnbr}
        del capnbr
        del fcoordnbr
    
    # Read distance file
    dist = dict() # distance dict indexed by (origin, destination) ID pairs
    with open(distfile, 'r') as f:
        for line in f:
            if line[0].isdigit() == False:
                continue
            s = line.strip().split()
            dist[(int(s[0]), int(s[1]))] = float(s[2])
    
    # Compute the facility crowdedness metrics
    fmet = dict() # facility crowdedness metrics by facfile index
    if crowding == True:
        print("Computing facility crowdedness metrics using distance file.")
        for j in tqdm.tqdm(cap):
            fmet[j] = 0.0
            for k in pop:
                d = max(dist[(k,j)], floor) # travel time (min) bounded by floor
                fmet[j] += pop[k]*(d**(-beta))
    else:
        for j in cap:
            fmet[j] = 1.0
    
    # Compute the population accessibility metrics
    pmet = dict() # population accessibility metrics by popfile index
    print("Computing population accessibility metrics using distance file.")
    for i in tqdm.tqdm(pop):
        pmet[i] = 0.0
        for j in cap:
            d = max(dist[(i,j)], floor) # travel time (min) bounded by floor
            pmet[i] += (cap[j]*(d**(-beta)))/fmet[j]
    
    # Return facility and population metric dictionaries (original indices only)
    return ({k: fmet[k] for k in fkeys}, {k: pmet[k] for k in pkeys})
    
    # Return facility and population metric dictionaries
    return (None, None)

#==============================================================================
# 2SFCA Metric Scripts
#==============================================================================

def fca_metric(poutfile, foutfile, popfile, facfile, distfile=None,
               cutoff=30.0, popnbrfile=None, facnbrfile=None, crowding=True,
               piecewise=0.5, speed=45.0, schedfile=None):
    """Computes a table of 2SFCA metrics for a given community.
    
    Positional arguments:
        poutfile (str) -- Output population file path.
        foutfile (str) -- Output facility file path.
        popfile (str) -- Preprocessed population file path, which should
            include the coordinates and population of each population center.
        facfile (str) -- Preprocessed facility file path, which should include
            the coordinates and capacity of each vaccination facility.
    
    Keyword arguments:
        distfile (str) -- Preprocessed distance file path, which should include
            the travel times between each population center/facility pair.
            Defaults to None, in which case geodesic distances are computed and
            used as needed.
        cutoff (float|tuple(float)) -- Travel time cutoff for defining
            catchment areas (in minutes). Defaults to 30.0. If given a 2-tuple
            of cutoff times, the first is treated as a cutoff time for urban
            tracts while the second is for rural tracts.
        popnbrfile (str) -- Preprocessed neighboring county population file
            path. Defaults to None.
        facnbrfile (str) -- Preprocessed neighboring county facility file path.
            Defaults to None.
        crowding (bool) -- Whether or not to take crowding into consideration.
            Defaults to True.
        piecewise (None|float) -- Controls how to treat the urban/rural travel
            time cutoffs (and thus only relevant if "cutoff" is given as a
            tuple). Defaults to 0.5. If given a float between 0.0 and 1.0, this
            is treated as an urban percentage threshold above which the urban
            travel time should be used, and below which the rural one should be
            used. If "None", the cutoff is instead taken as a convex
            combination of the urban and rural travel time cutoffs based on
            the tracts urban/rural split.
        speed (float) -- Assumed travel speed (mph). Defaults to 45.0. Only
            needed if using geodesic distance, since otherwise the distance
            file directly provides the travel times.
        schedfile (str) -- Schedule file, containing a row for each facility
            and columns to indicate availability during various time slots.
            Defaults to None. If a schedule file is included, the output file
            formats are changed to include multiple crowding/access columns to
            correspond to each time slot.
    
    This function implements the two-step floating catchment area (2SFCA)
    metric as described in Luo and Wang 2003 (doi:10.1068/b29120). We begin by
    computing a "crowdedness" metric for each facility by dividing its
    capacity by the population within that facility's catchment area
        
        R_j = S_j/(sum_{k:d_kj<=d_0} P_k)
    
    where R_j is the crowdedness metric for facility j, P_k is the population
    of center k, d_kj is the distance (or travel time) from k to j, and d_0 is
    a distance (or travel time) cutoff used to define catchment areas. A
    larger value of d_0 results in larger catchment area definitions.
    
    We then compute the accessibiility metric for each community as a sum of
    crowdedness metrics for all facilities within the community's catchment
    area
        
        A_i = sum_{j:d_ij<=d_0} R_j
    
    where A_i is the metric of population center i and the remaining notation
    is the same as above. The metric is called "two-step" because of how R_j
    and A_i are computed in a two-step process.
    
    The formulas above describe the default behavior of this script. If the
    "crowding" flag is set to False, an alternate crowding-free metric will be
    computed for each population center by setting R_j = S_j for all
    facilities, resulting in accessibility metrics that represent counts of
    facility capacities.
    
    Optionally, population centers and vaccination facilities in neighboring
    counties can be provided, in which case both will be used in computing
    crowding and accessibility metrics. This can be used to combat edge
    effects.
    
    Additional options also exist to use different travel time cutoffs
    depending on how much of the tract's population is urban versus rural. If
    the "cutoff" argument is given as a 2-tuple, its elements will be
    interpreted as an urban travel time cutoff d_U and a rural travel time
    cutoff d_R, respectively. The "piecewise" argument determines how to treat
    these cutoffs.
    
    To be more precise, let u_i be the fraction of tract i that is urban. If
    the "piecewise" argument is given a float between 0.0 and 1.0, then the
    number is treated as a threshold for u_i above which the urban travel time
    cutoff is used and below which the rural one is used. For example, the
    argument "cutoff=0.5" results in using a cutoff time for tract i of
        
              { d_U     if u_i >= 0.5
        d_0 = {
              { d_R     if u_i < 0.5
    
    If instead the "piecewise" argument is set to "None", the travel time
    cutoff for tract i is taken as a convex combination of d_U and d_R
        
        d_0 = u_i*d_U + (1-u_i)*d_R
    
    Optionally, schedule files can also be included. A schedule file should
    contain a row for each facility and a sequence of columns corresponding to
    different time slots, with each entry being a float between 0.0 and 1.0
    indicating the fraction of that time slot during which the facility is
    open. Each facility's capacity will be weighted by that amount (e.g. 0.0
    indicates complete closure and zero capacity, 1.0 indicates full capacity,
    and 0.5 indicates being open during half the time slot and thus half
    capacity). Including such a schedule file changes the facility and
    population file output formats to include multiple crowding or access
    columns: one for each time slot.
    
    The output files are augmented versions of the input files, with a new
    column of crowdedness or accessibility metrics appended to the original
    tables.
    """
    
    # Call an appropriate subroutine depending on distance specification
    if distfile == None:
        (fmet, pmet) = _fca_metric_geodesic(popfile, facfile, cutoff=cutoff,
                                            popnbrfile=popnbrfile,
                                            facnbrfile=facnbrfile,
                                            crowding=crowding,
                                            piecewise=piecewise,
                                            speed=speed,
                                            schedfile=schedfile)
    else:
        (fmet, pmet) = _fca_metric_file(popfile, facfile, distfile,
                                        cutoff=cutoff, popnbrfile=popnbrfile,
                                        facnbrfile=facnbrfile,
                                        crowding=crowding,
                                        piecewise=piecewise,
                                        schedfile=schedfile)
    
    # Get list of schedule time slots
    if schedfile == None:
        slots = ["crowding"]
        for i in pmet:
            pmet[i] = [pmet[i]]
        for j in fmet:
            fmet[j] = [fmet[j]]
    else:
        sdic = _read_schedfile(schedfile)
        slots = list(sdic.keys())
        for i in range(len(slots)):
            slots[i] = "crowding_" + slots[i]
        del sdic
    
    # Write facility output file with a new crowdedness metric column
    if crowding == True or schedfile != None:
        print("Writing facility metric file.")
        for k in range(len(slots)):
            col = dict() # new column of entries
            for j in fmet:
                col[j] = fmet[j][k]
            if k == 0:
                _augment_file(foutfile, facfile, col, slots[k])
            else:
                _augment_file(foutfile, foutfile, col, slots[k])
  
    # Update column labels
    for i in range(len(slots)):
        slots[i] = slots[i].replace("crowding", "access")
    
    # Write population output file with a new accessibility metric column
    print("Writing population metric file.")
    for k in range(len(slots)):
        col = dict() # new column of entries
        for i in pmet:
            col[i] = pmet[i][k]
        if k == 0:
            _augment_file(poutfile, popfile, col, slots[k])
        else:
            _augment_file(poutfile, poutfile, col, slots[k])

#------------------------------------------------------------------------------

def _fca_metric_geodesic(popfile, facfile, cutoff=30.0, popnbrfile=None,
                         facnbrfile=None, crowding=True, piecewise=0.5,
                         speed=45.0, schedfile=None):
    """The geodesic distance version of fca_metric.
    
    This script computes 2SFCA metrics using geodesic distances computed as-
    needed. It is automatically called by fca_metric when a distance file is
    not provided. It returns dictionaries of facility and population metrics
    rather than writing the results directly to a file.
    
    Geodesic distances are used to compute travel times in minutes, assuming a
    constant speed (in mph).
    """
    
    ### Schedule file support not yet included
    if schedfile == None:
        raise NotImplementedError("schedule files not implemented for geodesic"
                                  + " metrics")
    
    # Convert speed to miles per minute
    speed /= 60.0
    
    # Read population file
    (pop, pcoord, _, urban) = _read_popfile(popfile)
    pkeys = list(pop.keys()) # store main population keys
    
    # If a neighboring population file is provided, merge its contents
    if popnbrfile != None:
        (popnbr, pcoordnbr) = _read_popfile(popnbrfile)
        pop = {**pop, **popnbr}
        pcoord = {**pcoord, **pcoordnbr}
        del popnbr
        del pcoordnbr
    
    # Define a travel time cutoff for each tract
    d0 = dict()
    if isinstance(cutoff, tuple) and len(cutoff) > 1:
        # If the cutoff is at least a 2-tuple, use variable cutoffs
        if piecewise == None:
            # Convex combination
            for k in pop:
                d0[k] = _convex(urban[k], cutoff[0], cutoff[1])
        else:
            # Binary urban/rural threshold
            for k in pop:
                d0[k] = _threshold(urban[k], piecewise, cutoff[0], cutoff[1])
    
    # Read facility file
    (cap, fcoord) = _read_facfile(facfile)
    fkeys = list(cap.keys()) # store main facility keys
    
    # If a neighboring facility file is provided, merge its contents
    if facnbrfile != None:
        (capnbr, fcoordnbr) = _read_popfile(facnbrfile)
        cap = {**cap, **capnbr}
        fcoord = {**fcoord, **fcoordnbr}
        del capnbr
        del fcoordnbr
    
    # Compute the facility crowdedness metrics
    fmet = dict() # facility crowdedness metrics by facfile index
    if crowding == True:
        print("Computing facility crowdedness metrics using geodesic distance.")
        for j in tqdm.tqdm(cap):
            p = 0 # total population in range
            for k in pop:
                d = geodesic_distance(pcoord[k], fcoord[j])/speed # time (min)
                if d <= d0[k]:
                    p += pop[k]
            if p > 0:
                fmet[j] = cap[j]/p
            else:
                fmet[j] = -1
    else:
        for j in cap:
            fmet[j] = cap[j]
    
    # Compute the population accessibility metrics
    pmet = dict() # population accessibility metrics by popfile index
    print("Computing population accessibility metrics using geodesic distance.")
    for i in tqdm.tqdm(pop):
        pmet[i] = 0.0
        for j in cap:
            if fmet[j] < 0:
                continue
            d = geodesic_distance(pcoord[i], fcoord[j])/speed # time (minutes)
            if d <= d0[i]:
                pmet[i] += fmet[j]
    
    # Return facility and population metric dictionaries (original indices only)
    return ({k: fmet[k] for k in fkeys}, {k: pmet[k] for k in pkeys})

#------------------------------------------------------------------------------

def _fca_metric_file(popfile, facfile, distfile, cutoff=30.0, popnbrfile=None,
                     facnbrfile=None, crowding=True, piecewise=0.5,
                     schedfile=None):
    """The distance file version of fca_metric.
    
    This script computes FCA metrics using a predefined distance file. It
    is automatically called by fca_metric when a distance file is provided.
    It returns dictionaries of facility and population metrics rather than
    writing the results directly to a file.
    
    If "crowding" is set to False, the metrics produced are equivalent to
    simple counts of facilities (weighted by capacity and possibly also
    schedule availability).
    
    The output facility and population dictionaries both include an entry for
    each facility or population ID defined in the input files. Each dictionary
    entry is a list of crowding or accessibility metrics, respectively, for each
    time listed in the schedule file (or just a list with a single element if
    no schedule file was specified).
    """
    
    # Read population file
    (pop, pcoord, _, urban) = _read_popfile(popfile)
    pkeys = list(pop.keys()) # store main population keys
    
    # If a neighboring population file is provided, merge its contents
    if popnbrfile != None:
        (popnbr, pcoordnbr, _, urbnbr) = _read_popfile(popnbrfile)
        pop = {**pop, **popnbr}
        pcoord = {**pcoord, **pcoordnbr}
        urban = {**urban, **urbnbr}
        del popnbr
        del pcoordnbr
        del urbnbr
    
    # Define a travel time cutoff for each tract
    d0 = dict()
    if isinstance(cutoff, tuple) and len(cutoff) > 1:
        # If the cutoff is at least a 2-tuple, use variable cutoffs
        if piecewise == None:
            # Convex combination
            for k in pop:
                d0[k] = _convex(urban[k], cutoff[0], cutoff[1])
        else:
            # Binary urban/rural threshold
            for k in pop:
                d0[k] = _threshold(urban[k], piecewise, cutoff[0], cutoff[1])
    else:
        # Otherwise use a constant cutoff
        for k in pop:
            d0[k] = float(cutoff)
    
    # Read facility file
    (cap, fcoord) = _read_facfile(facfile)
    fkeys = list(cap.keys()) # store main facility keys
    
    # If a neighboring facility file is provided, merge its contents
    if facnbrfile != None:
        (capnbr, fcoordnbr) = _read_facfile(facnbrfile)
        cap = {**cap, **capnbr}
        fcoord = {**fcoord, **fcoordnbr}
        del capnbr
        del fcoordnbr
    
    # Read distance file
    dist = dict() # distance dict indexed by (origin, destination) ID pairs
    with open(distfile, 'r') as f:
        for line in f:
            if line[0].isdigit() == False:
                continue
            s = line.strip().split()
            dist[(int(s[0]), int(s[1]))] = float(s[2])
    
    # Attempt to read schedule file (if any)
    if schedfile != None:
        # Read schedule file
        sdic = _read_schedfile(schedfile) # dict of dicts of availabilities
        slots = list(sdic.keys()) # list of time slot names
    else:
        # Construct dummy schedule dictionaries if no file is present
        slots = ["0"]
        sdic = {"0" : dict()}
        for j in cap:
            sdic["0"][j] = 1.0
    
    # Define crowding metrics for each time slot
    fmet = dict() # facility crowdedness metrics by facfile index
    for j in cap:
        fmet[j] = []
    print("Computing facility crowdedness metrics using distance file.")
    for i in range(len(slots)):
        
        if crowding == True:
            for j in cap:
                p = 0 # total population in range
                for k in pop:
                    d = dist[(k,j)] # travel time (min)
                    if d <= d0[k]:
                        p += pop[k]
                if p > 0:
                    fmet[j] += [cap[j]/p] # capacity over population
                else:
                    fmet[j] += [-1]
                fmet[j][-1] *= sdic[slots[i]][j] # schedule adjustment
        else:
            for j in cap:
                fmet[j] += [cap[j]] # facility capacity (if no crowding)
                fmet[j][-1] *= sdic[slots[i]][j] # schedule adjustment
    
    # Define accessibility metrics for each facility crowding metric
    pmet = dict() # population accessibility metric list by popfile index
    for i in pop:
        pmet[i] = []
    print("Computing population accessibility metrics using distance file.")
    for k in tqdm.tqdm(range(len(slots))):
        for i in pop:
            pmet[i] += [0.0] # append a new metric
            for j in cap:
                if fmet[j][k] < 0:
                    continue
                d = dist[(i,j)] # travel time (min)
                if d <= d0[i]:
                    pmet[i][-1] += fmet[j][k]
    
    # If no schedule file was included, return floats instead of lists
    if schedfile == None:
        for j in fkeys:
            fmet[j] = fmet[j][0]
        for i in pop:
            pmet[i] = pmet[i][0]
    
    # Return facility and population metric dictionaries (original indices only)
    return ({k: fmet[k] for k in fkeys}, {k: pmet[k] for k in pkeys})

#==============================================================================
# Execution
#==============================================================================

# Comment or uncomment the function calls below to process each location.
#bl = {0.50: "0-50", 0.75: "0-75", 1.00: "1-00", 1.25: "1-25", 1.50: "1-50", 1.75: "1-75", 2.00: "2-00"}
#for beta in bl:
#    gravity_metric(os.path.join("..", "results", "santa_clara", "santa_clara_pop_gravity_" + bl[beta] +".tsv"), os.path.join("..", "results", "santa_clara", "santa_clara_fac_gravity_" + bl[beta] + ".tsv"), os.path.join("..", "processed", "santa_clara", "santa_clara_pop.tsv"), os.path.join("..", "processed", "santa_clara", "santa_clara_fac.tsv"), distfile=None, beta=beta)
#for beta in bl:
#    gravity_metric(os.path.join("..", "results", "santa_clara", "santa_clara_pop_gravity_nocrowding_" + bl[beta] +".tsv"), os.path.join("..", "results", "santa_clara", "santa_clara_fac_gravity_nocrowding_" + bl[beta] + ".tsv"), os.path.join("..", "processed", "santa_clara", "santa_clara_pop.tsv"), os.path.join("..", "processed", "santa_clara", "santa_clara_fac.tsv"), distfile=None, beta=beta, crowding=False)
#co = {5: "005", 10: "010", 15: "015", 20: "020", 25: "025", 30: "030", 35: "035", 40: "040", 45: "045", 50: "050", 55: "055", 60: "060"}
#for cutoff in co:
#    fca_metric(os.path.join("..", "results", "santa_clara", "santa_clara_pop_cutoff_" + co[cutoff] + ".tsv"), os.path.join("..", "results", "santa_clara", "santa_clara_fac_cutoff_" + co[cutoff] + ".tsv"), os.path.join("..", "processed", "santa_clara", "santa_clara_pop.tsv"), os.path.join("..", "processed", "santa_clara", "santa_clara_fac.tsv"), distfile=None, cutoff=cutoff)
#for cutoff in co:
#    fca_metric(os.path.join("..", "results", "santa_clara", "santa_clara_pop_cutoff_nocrowding_" + co[cutoff] + ".tsv"), os.path.join("..", "results", "santa_clara", "santa_clara_fac_cutoff_nocrowding_" + co[cutoff] + ".tsv"), os.path.join("..", "processed", "santa_clara", "santa_clara_pop.tsv"), os.path.join("..", "processed", "santa_clara", "santa_clara_fac.tsv"), distfile=None, cutoff=cutoff, crowding=False)

#popfile = os.path.join("..", "processed", "santa_clara", "santa_clara_pop.tsv")
#popnbrfile = os.path.join("..", "processed", "santa_clara", "santa_clara_pop_nbr.tsv")
#facfile = os.path.join("..", "processed", "santa_clara", "santa_clara_fac.tsv")
#facnbrfile = os.path.join("..", "processed", "santa_clara", "santa_clara_fac_nbr.tsv")
#distfile = os.path.join("..", "processed", "santa_clara", "santa_clara_dist.tsv")

#poutfile = os.path.join("..", "results", "santa_clara", "santa_clara_pop_gravity_1-00_geo.tsv")
#foutfile = os.path.join("..", "results", "santa_clara", "santa_clara_fac_gravity_1-00_geo.tsv")
#gravity_metric(poutfile, foutfile, popfile, facfile, beta=1.0, popnbrfile=popnbrfile, facnbrfile=facnbrfile, crowding=True)
#poutfile = os.path.join("..", "results", "santa_clara", "santa_clara_pop_gravity_1-00_geo_nocrowding.tsv")
#gravity_metric(poutfile, foutfile, popfile, facfile, beta=1.0, popnbrfile=popnbrfile, facnbrfile=facnbrfile, crowding=False)
#poutfile = os.path.join("..", "results", "santa_clara", "santa_clara_pop_gravity_1-00_tt.tsv")
#foutfile = os.path.join("..", "results", "santa_clara", "santa_clara_fac_gravity_1-00_tt.tsv")
#gravity_metric(poutfile, foutfile, popfile, facfile, beta=1.0, popnbrfile=popnbrfile, facnbrfile=facnbrfile, crowding=True, distfile=distfile, floor=1.0)
#poutfile = os.path.join("..", "results", "santa_clara", "santa_clara_pop_gravity_1-00_tt_nocrowding.tsv")
#gravity_metric(poutfile, foutfile, popfile, facfile, beta=1.0, popnbrfile=popnbrfile, facnbrfile=facnbrfile, crowding=False, distfile=distfile, floor=1.0)

#poutfile = os.path.join("..", "results", "santa_clara", "santa_clara_pop_fca_030_geo.tsv")
#foutfile = os.path.join("..", "results", "santa_clara", "santa_clara_fac_fca_030_geo.tsv")
#fca_metric(poutfile, foutfile, popfile, facfile, cutoff=30.0, popnbrfile=popnbrfile, facnbrfile=facnbrfile, crowding=True)
#poutfile = os.path.join("..", "results", "santa_clara", "santa_clara_pop_fca_030_geo_nocrowding.tsv")
#fca_metric(poutfile, foutfile, popfile, facfile, cutoff=30.0, popnbrfile=popnbrfile, facnbrfile=facnbrfile, crowding=False)
#poutfile = os.path.join("..", "results", "santa_clara", "santa_clara_pop_fca_015-030_convex.tsv")
#foutfile = os.path.join("..", "results", "santa_clara", "santa_clara_fac_fca_015-030_convex.tsv")
#fca_metric(poutfile, foutfile, popfile, facfile, cutoff=(15.0, 30.0), piecewise=None, popnbrfile=popnbrfile, facnbrfile=facnbrfile, crowding=True, distfile=distfile)
#poutfile = os.path.join("..", "results", "santa_clara", "santa_clara_pop_fca_030_tt_nocrowding.tsv")
#fca_metric(poutfile, foutfile, popfile, facfile, cutoff=30.0, popnbrfile=popnbrfile, facnbrfile=facnbrfile, crowding=False, distfile=distfile)
#poutfile = os.path.join("..", "results", "santa_clara", "santa_clara_pop_fca_015-030_cutoff_50.tsv")
#foutfile = os.path.join("..", "results", "santa_clara", "santa_clara_fac_fca_015-030_cutoff_50.tsv")
#fca_metric(poutfile, foutfile, popfile, facfile, cutoff=(15.0, 30.0), piecewise=0.5, popnbrfile=popnbrfile, facnbrfile=facnbrfile, crowding=True, distfile=distfile)
#poutfile = os.path.join("..", "results", "santa_clara", "santa_clara_pop_fca_015-030_cutoff_25.tsv")
#foutfile = os.path.join("..", "results", "santa_clara", "santa_clara_fac_fca_015-030_cutoff_25.tsv")
#fca_metric(poutfile, foutfile, popfile, facfile, cutoff=(15.0, 30.0), piecewise=0.25, popnbrfile=popnbrfile, facnbrfile=facnbrfile, crowding=True, distfile=distfile)
#poutfile = os.path.join("..", "results", "santa_clara", "santa_clara_pop_fca_015-030_cutoff_75.tsv")
#foutfile = os.path.join("..", "results", "santa_clara", "santa_clara_fac_fca_015-030_cutoff_75.tsv")
#fca_metric(poutfile, foutfile, popfile, facfile, cutoff=(15.0, 30.0), piecewise=0.75, popnbrfile=popnbrfile, facnbrfile=facnbrfile, crowding=True, distfile=distfile)

#------------------------------------------------------------------------------

# Polk County urgent care access

popfile = os.path.join("..", "processed", "polk", "polk_pop.tsv")
popnbrfile = os.path.join("..", "processed", "polk", "polk_pop_nbr.tsv")

facfile = os.path.join("..", "processed", "polk", "polk_uc.tsv")
facnbrfile = os.path.join("..", "processed", "polk", "polk_uc_nbr.tsv")
distfile = os.path.join("..", "processed", "polk", "polk_dist_uc.tsv")
schedfilefull = os.path.join("..", "processed", "polk", "polk_schedule_uc.tsv")
schedfileabbrv = os.path.join("..", "processed", "polk", "polk_schedule_abbrv_uc.tsv")

#poutfile = os.path.join("..", "results", "polk", "polk_pop_uc_count_15-cutoff_noschedule.tsv")
#foutfile = os.path.join("..", "results", "polk", "polk_fac_uc_count_15-cutoff_noschedule.tsv")
#fca_metric(poutfile, foutfile, popfile, facfile, cutoff=15.0, piecewise=None, popnbrfile=popnbrfile, facnbrfile=facnbrfile, crowding=False, distfile=distfile, schedfile=None)

#poutfile = os.path.join("..", "results", "polk", "polk_pop_uc_count_15-cutoff_fullschedule.tsv")
#foutfile = os.path.join("..", "results", "polk", "polk_fac_uc_count_15-cutoff_fullschedule.tsv")
#fca_metric(poutfile, foutfile, popfile, facfile, cutoff=15.0, piecewise=None, popnbrfile=popnbrfile, facnbrfile=facnbrfile, crowding=False, distfile=distfile, schedfile=schedfilefull)

#poutfile = os.path.join("..", "results", "polk", "polk_pop_uc_count_15-cutoff_abbrvschedule.tsv")
#foutfile = os.path.join("..", "results", "polk", "polk_fac_uc_count_15-cutoff_abbrvschedule.tsv")
#fca_metric(poutfile, foutfile, popfile, facfile, cutoff=15.0, piecewise=None, popnbrfile=popnbrfile, facnbrfile=facnbrfile, crowding=False, distfile=distfile, schedfile=schedfileabbrv)

#poutfile = os.path.join("..", "results", "polk", "polk_pop_uc_count_30-cutoff_noschedule.tsv")
#foutfile = os.path.join("..", "results", "polk", "polk_fac_uc_count_30-cutoff_noschedule.tsv")
#fca_metric(poutfile, foutfile, popfile, facfile, cutoff=30.0, piecewise=None, popnbrfile=popnbrfile, facnbrfile=facnbrfile, crowding=False, distfile=distfile, schedfile=None)

#poutfile = os.path.join("..", "results", "polk", "polk_pop_uc_count_30-cutoff_fullschedule.tsv")
#foutfile = os.path.join("..", "results", "polk", "polk_fac_uc_count_30-cutoff_fullschedule.tsv")
#fca_metric(poutfile, foutfile, popfile, facfile, cutoff=30.0, piecewise=None, popnbrfile=popnbrfile, facnbrfile=facnbrfile, crowding=False, distfile=distfile, schedfile=schedfilefull)

#poutfile = os.path.join("..", "results", "polk", "polk_pop_uc_count_30-cutoff_abbrvschedule.tsv")
#foutfile = os.path.join("..", "results", "polk", "polk_fac_uc_count_30-cutoff_abbrvschedule.tsv")
#fca_metric(poutfile, foutfile, popfile, facfile, cutoff=30.0, piecewise=None, popnbrfile=popnbrfile, facnbrfile=facnbrfile, crowding=False, distfile=distfile, schedfile=schedfileabbrv)

# Travel time cutoffs to consider, respectively:
# 15 minute cutoff
# 30 minute cutoff
# 15-30 minute convex combination by urban/rural split
# 15-30 minute cutoff if more than 50% urban

cutoffnames = ["15-cutoff", "30-cutoff", "15-30-urban-convex", "15-30-urban-piecewise"]
cutoffs = [15.0, 30.0, (15.0, 30.0), (15.0, 30.0)]
piecewise = [None, None, None, 0.5]

# Schedules to consideer, respectively:
# None (just a facility count regardless of shcedule)
# Abbreviated

schedules = [None, schedfileabbrv]
schednames = ["all-times", "day-type-times"]

for i in range(len(cutoffnames)):
    for j in range(len(schedules)):
        poutfile = os.path.join("..", "results", "polk", "polk_pop_uc_count_" + cutoffnames[i] + '_' + schednames[j] + ".tsv")
        foutfile = os.path.join("..", "results", "polk", "polk_fac_uc_count_" + cutoffnames[i] + '_' + schednames[j] + ".tsv")
        fca_metric(poutfile, foutfile, popfile, facfile, cutoff=cutoffs[i], piecewise=piecewise[i], popnbrfile=popnbrfile, facnbrfile=facnbrfile, crowding=False, distfile=distfile, schedfile=schedules[j])

#------------------------------------------------------------------------------

# Polk County pharmacy access

popfile = os.path.join("..", "processed", "polk", "polk_pop.tsv")
popnbrfile = os.path.join("..", "processed", "polk", "polk_pop_nbr.tsv")

facfile = os.path.join("..", "processed", "polk", "polk_pharmacy.tsv")
facnbrfile = os.path.join("..", "processed", "polk", "polk_pharmacy_nbr.tsv")
distfile = os.path.join("..", "processed", "polk", "polk_dist_pharmacy.tsv")
schedfilefull = os.path.join("..", "processed", "polk", "polk_schedule_pharmacy.tsv")
schedfileabbrv = os.path.join("..", "processed", "polk", "polk_schedule_abbrv_pharmacy.tsv")

#poutfile = os.path.join("..", "results", "polk", "polk_pop_pharm_count_15-cutoff_noschedule.tsv")
#foutfile = os.path.join("..", "results", "polk", "polk_fac_pharm_count_15-cutoff_noschedule.tsv")
#fca_metric(poutfile, foutfile, popfile, facfile, cutoff=15.0, piecewise=None, popnbrfile=popnbrfile, facnbrfile=facnbrfile, crowding=False, distfile=distfile, schedfile=None)

#poutfile = os.path.join("..", "results", "polk", "polk_pop_pharm_count_15-cutoff_fullschedule.tsv")
#foutfile = os.path.join("..", "results", "polk", "polk_fac_pharm_count_15-cutoff_fullschedule.tsv")
#fca_metric(poutfile, foutfile, popfile, facfile, cutoff=15.0, piecewise=None, popnbrfile=popnbrfile, facnbrfile=facnbrfile, crowding=False, distfile=distfile, schedfile=schedfilefull)

#poutfile = os.path.join("..", "results", "polk", "polk_pop_pharm_count_15-cutoff_abbrvschedule.tsv")
#foutfile = os.path.join("..", "results", "polk", "polk_fac_pharm_count_15-cutoff_abbrvschedule.tsv")
#fca_metric(poutfile, foutfile, popfile, facfile, cutoff=15.0, piecewise=None, popnbrfile=popnbrfile, facnbrfile=facnbrfile, crowding=False, distfile=distfile, schedfile=schedfileabbrv)

#poutfile = os.path.join("..", "results", "polk", "polk_pop_pharm_count_30-cutoff_noschedule.tsv")
#foutfile = os.path.join("..", "results", "polk", "polk_fac_pharm_count_30-cutoff_noschedule.tsv")
#fca_metric(poutfile, foutfile, popfile, facfile, cutoff=30.0, piecewise=None, popnbrfile=popnbrfile, facnbrfile=facnbrfile, crowding=False, distfile=distfile, schedfile=None)

#poutfile = os.path.join("..", "results", "polk", "polk_pop_pharm_count_30-cutoff_fullschedule.tsv")
#foutfile = os.path.join("..", "results", "polk", "polk_fac_pharm_count_30-cutoff_fullschedule.tsv")
#fca_metric(poutfile, foutfile, popfile, facfile, cutoff=30.0, piecewise=None, popnbrfile=popnbrfile, facnbrfile=facnbrfile, crowding=False, distfile=distfile, schedfile=schedfilefull)

#poutfile = os.path.join("..", "results", "polk", "polk_pop_pharm_count_30-cutoff_abbrvschedule.tsv")
#foutfile = os.path.join("..", "results", "polk", "polk_fac_pharm_count_30-cutoff_abbrvschedule.tsv")
#fca_metric(poutfile, foutfile, popfile, facfile, cutoff=30.0, piecewise=None, popnbrfile=popnbrfile, facnbrfile=facnbrfile, crowding=False, distfile=distfile, schedfile=schedfileabbrv)

# Travel time cutoffs to consider, respectively:
# 15 minute cutoff
# 30 minute cutoff
# 15-30 minute convex combination by urban/rural split
# 15-30 minute cutoff if more than 50% urban

cutoffnames = ["15-cutoff", "30-cutoff", "15-30-urban-convex", "15-30-urban-piecewise"]
cutoffs = [15.0, 30.0, (15.0, 30.0), (15.0, 30.0)]
piecewise = [None, None, None, 0.5]

# Schedules to consideer, respectively:
# None (just a facility count regardless of schedule)
# Abbreviated

schedules = [None, schedfileabbrv]
schednames = ["all-times", "day-type-times"]

for i in range(len(cutoffnames)):
    for j in range(len(schedules)):
        poutfile = os.path.join("..", "results", "polk", "polk_pop_pharm_count_" + cutoffnames[i] + '_' + schednames[j] + ".tsv")
        foutfile = os.path.join("..", "results", "polk", "polk_fac_pharm_count_" + cutoffnames[i] + '_' + schednames[j] + ".tsv")
        fca_metric(poutfile, foutfile, popfile, facfile, cutoff=cutoffs[i], piecewise=piecewise[i], popnbrfile=popnbrfile, facnbrfile=facnbrfile, crowding=False, distfile=distfile, schedfile=schedules[j])
