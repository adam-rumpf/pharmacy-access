"""Pharmacy accessibility project graphics generation scripts.

The code below is meant for generating graphics for use in data analysis and
result display. This includes both maps to display geospatial distributions as
well as plots for statistical results.
"""

import csv
import os.path

import geopandas as gpd
import geopy.distance
import matplotlib.patches
import matplotlib.pyplot as plt
import matplotlib_scalebar.scalebar as scb
import pandas as pd
import pygris as pg
import scipy.stats as stat
import shapely as shp

#==============================================================================
# Global Constants
#==============================================================================

# Shapefile names and paths
SHP_FL_COUNTIES = os.path.join("..", "shapefiles", "fl", "fl_counties.shp")
SHP_POLK_TRACTS = os.path.join("..", "shapefiles", "polk", "polk_tracts.shp")
SHP_POLK_BLOCKS = os.path.join("..", "shapefiles", "polk", "polk_blocks.shp")

# Polk County result file directory
POLK_PROCESSED = os.path.join("..", "processed", "polk")
POLK_RESULTS = os.path.join("..", "results", "polk")

# Polk County bounding boxes
POLK_X_FULL = (-82.155, -81.087)
POLK_Y_FULL = (27.608, 28.397)

# Lakeland bounding boxes
LAKELAND_X = (-82.057, -81.798)
LAKELAND_Y = (27.949, 28.197)

# Distance of 1 degree near Polk County (meters)
POLK_DEGREE = geopy.distance.geodesic((POLK_X_FULL[0], POLK_Y_FULL[1]),
                                      (POLK_X_FULL[0]+1, POLK_Y_FULL[1])).meters

# Florida Polytechnic University coordinates
FL_POLY = (-81.8513, 28.1508)

# Binary indicator color map
white_red = matplotlib.colors.ListedColormap(["white", "red"])

#==============================================================================
# Data Collection and Arrangement
#==============================================================================

def download_shapefiles(state, county, cfile=None, tfile=None, bfile=None):
    """Downloads and saves shapefiles within a given county.
    
    Positional arguments:
        state (str) -- Two-letter state abbreviation.
        county (str) -- Name of county.
    
    Keyword arguments:
        cfile (str) -- File path for saving the county shapefile. Defaults to
            "None", in which case the county shapefile is skipped.
        tfile (str) -- File path for saving the tract shapefile. Defaults to
            "None", in which case the tract shapefile is skipped.
        bfile (str) -- File path for saving the block shapefile. Defaults to
            "None", in which case the county block is skipped.
    """
    
    # Download each of the specified files
    if cfile != None:
        counties = pg.counties(state=state)
        counties.to_file(cfile)
        del counties
    if tfile != None:
        tracts = pg.tracts(state=state, county=county)
        tracts.to_file(tfile)
        del tracts
    if bfile != None:
        blocks = pg.blocks(state=state, county=county)
        blocks.to_file(bfile)
        del blocks

#------------------------------------------------------------------------------

def _make_geodata(shapefile, datafile, field, shapeid="GEOID", dataid="name",
                 missing=None):
    """Adds additional data fields to a GeoDataFrame.
    
    Positional arguments:
        shapefile (str) -- File path for the shapefile (tract, block, etc.) to
            be used as the basis of the GeoDataFrame.
        datafile (str) -- File path for a data file containing fields to be
            added to the shapefile's fields.
        field (str) -- Fields in the data file to be added to the GeoDataFrame.
    
    Keyword arguments:
        shapeid (str) -- Field in the shapefile to use as the unique name to
            distinguish each field. Defaults to "GEOID", which is the full FIPS
            code for the TIGER/Line shapefiles.
        dataid (str) -- Field in the data file to use as the unique name to
            distinguish each field. Should contain exactly the same sets of
            names as appear in the shapefile's "shapeid" field. Defaults to
            "name", which is what we've used for our population result files.
        missing -- Default entry used to fill missing fields. Default None.
    
    Returns:
        (GeoDataFrame) -- A GeoPandas GeoDataFrame object containing the fields
            from the base file as well as any additional provided fields.
    
    To explain, the shapefile and the data file are both expected to have the
    same number of rows, both containing at least one column full of unique IDs
    (e.g. FIPS codes) to fully distinguish each row. This function uses the IDs
    from the shapefile to rearrange the rows from the data file as necessary so
    that all rows correspond correctly, and then it adds any specified columns
    to the shapefile, returning the resulting GeoDataFrame.
    
    The purpose of this function is to add result data to a data frame made from
    a shapefile for use in generating heat maps.
    """
    
    # Get shapefile dataframe and IDs
    shapedf = gpd.read_file(shapefile)
    shapeids = list(shapedf.get(shapeid))
    
    # Get elements of data file field as a dictionary
    data = dict()
    with open(datafile, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t', quotechar='"')
        for row in reader:
            data[row[dataid]] = float(row[field])
    
    # Generate a new column for the dataframe by matching corresponding IDs
    col = [missing for i in range(len(shapeids))]
    for i in range(len(shapeids)):
        if shapeids[i] in data:
            col[i] = data[shapeids[i]]
    
    # Add the column to the dataframe and return the result
    shapedf[field] = col
    return shapedf

#------------------------------------------------------------------------------

def _get_field(datafile, field, cast=2):
    """Gets a specified field from a data file.
    
    Positional arguments:
        datafile (str) -- File path for a data file.
        field (str) -- Field in the data file.
    
    Positional keyword arguments:
        cast (int) -- Option for how to cast the contents of the field. Options
            include: 0 - string, 1 - integer, 2 - float. Default 2 (float).
    
    Returns:
        (list) -- The contents of the specified field of the data file, as a
            list.
    """
    
    # Get field from data file
    col = []
    with open(datafile, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t', quotechar='"')
        for row in reader:
            if cast == 1:
                val = int(row[field])
            elif cast == 2:
                val = float(row[field])
            else:
                val = row[field]
            col.append(val)
    
    return col

#==============================================================================
# Map Generation
#==============================================================================

def map_county(fname, axes, cfips, cfipsname="COUNTYFP", color="black"):
    """Adds a county border in a given file to a given set of axes.
    
    Positional arguments:
        fname (str) -- File path to a data file containing county shapefiles.
        axes (ax) -- Matplotlib axis object to add this plot to.
        cfips (str) -- The 3-digit county FIPS code.
    
    Keyword arguments:
        cfipsname (str) -- Name of the field containing county FIPS codes.
            Defaults to "COUNTYFP".
        color (str) -- Color of county border. Defaults to "black".
    """
    
    # Cast FIPS code as a string
    cfips = str(cfips)
    
    # Get county shapefile dataframe and extract specified county
    counties = gpd.read_file(fname)
    df = pd.DataFrame(counties)
    row = df.loc[df[cfipsname] == cfips]
    
    # Add the border to the axis
    gpd.GeoDataFrame(row).plot(ax=axes, facecolor="none", edgecolor=color)

#------------------------------------------------------------------------------

def map_shapefile(fname, axes, color="black"):
    """Adds shapefile borders in a given file to a given set of axes.
    
    Positional arguments:
        fname (str) -- File path to a data file containing shapefiles (e.g.
            tracts or blocks).
        axes (ax) -- Matplotlib axis object to add this plot to.
    
    Keyword arguments:
        color (str) -- Color of county border. Defaults to "black".
    """
    
    # Open and plot shapefile
    shapes = gpd.read_file(fname)
    shapes.plot(ax=axes, facecolor="none", edgecolor=color)

#------------------------------------------------------------------------------

def map_points(fname, axes, latname="lat", lonname="lon", color="black"):
    """Adds the points defined in a given file to a given set of axes.
    
    Positional arguments:
        fname (str) -- File path to a data file containing latitude and
            longitude fields.
        axes (ax) -- Matplotlib axis object to add this plot to.
    
    Keyword arguments:
        latname (str) -- Name of the field containing latitude values. Defaults
            to "lat".
        lonname (str) -- Name of the field containing longitude values. Defaults
            to "lon".
        color (str) -- Color of points. Defaults to "black".
    """
    
    # Get latitude and longitude values from file
    with open(fname, 'r') as f:
        reader = list(csv.DictReader(f, delimiter='\t', quotechar='"'))
        lat = [float(row[latname]) for row in reader]
        lon = [float(row[lonname]) for row in reader]
    
    # Create point list
    geom = [shp.Point(pt) for pt in zip(lon, lat)]
    df = gpd.GeoDataFrame(geometry=geom)
    
    # Add the points to the axis
    df.plot(ax=axes, color=color)

#------------------------------------------------------------------------------

def map_points_list(points, axes, color="black"):
    """Adds the points defined in a given list to a given set of axes.
    
    Positional arguments:
        points (list) -- List of coordinate tuples for the points.
        axes (ax) -- Matplotlib axis object to add this plot to.
    
    Keyword arguments:
        color (str) -- Color of points. Defaults to "black".
    """
    
    # Create point list
    geom = [shp.Point(pt) for pt in points]
    df = gpd.GeoDataFrame(geometry=geom)
    
    # Add the points to the axis
    df.plot(ax=axes, color=color)

#------------------------------------------------------------------------------

def map_rectangle(xlim, ylim, axes, color="black"):
    """Adds a polygon defined by a list of points to a given set of axes.
    
    Positional arguments:
        xlim (list(float)) -- Pair of points defining the boundaries in the
            x-direction.
        ylim (list(float)) -- Pair of points defining the boundaries in the
            y-direction.
        axes (ax) -- Matplotlib axis object to add this plot to.
    
    Keyword arguments:
        color (str) -- Color of polygon. Defaults to "black".
    """
    
    # Define rectangle lower corner, width, and height
    ll = min(xlim[0], xlim[1])
    bl = min(ylim[0], ylim[1])
    w = abs(xlim[0] - xlim[1])
    h = abs(ylim[0] - ylim[1])
    
    # Add rectangle to axes
    ax.add_patch(matplotlib.patches.Rectangle((ll, bl), w, h, edgecolor=color,
                                              fill=False))

#------------------------------------------------------------------------------

def map_heat(shapefile, datafile, field, axes, color="viridis", shapeid="GEOID",
             dataid="name", missing=None, legend=None, limits=None):
    """Adds a heat map based on a given field to a given set of axes.
    
    Positional arguments:
        shapefile (str) -- Path to the main shapefile.
        datafile (str) -- Path to the data file on which the heat map colors
            will be based.
        field (str) -- Name of field in data file to use for the heat map.
        axes (ax) -- Matplotlib axis object to add this plot to.
    
    Keyword arguments:
        color (str) -- Name of Matplotlib color scheme. Defaults to "viridis".
            The full list of options can be found here:
            https://matplotlib.org/stable/users/explain/colors/colormaps.html
        shapeid (str) -- Field in the shapefile to use as the unique name to
            distinguish each field. Defaults to "GEOID", which is the full FIPS
            code for the TIGER/Line shapefiles.
        dataid (str) -- Field in the data file to use as the unique name to
            distinguish each field. Should contain exactly the same sets of
            names as appear in the shapefile's "shapeid" field. Defaults to
            "name", which is what we've used for our population result files.
        missing -- Default entry used to fill missing fields. Default None.
        legend (str) -- Legend bar name. Default None, in which case no legend
            bar is generated. Simply specifying an empty string (legend="")
            creates the legend bar with no string.
        limits (float, float) -- Color map minimum and maximum values.
            Default None, in which case the range of colors is dictated
            entirely by the range of values in the data file. Including these
            limits expands the effective range of values in the plotted field
            to include the two specified values.
    """
    
    # Create a single unified GeoDataFrame
    frame = _make_geodata(shapefile, datafile, field, shapeid=shapeid,
                          dataid=dataid, missing=missing)
    
    # Add dummy fields in case an expanded range is required
    if limits != None:
        frame.loc[len(frame)] = [None if c is not field else limits[0]
                                 for c in list(frame.columns)]
        frame.loc[len(frame)] = [None if c is not field else limits[1]
                                 for c in list(frame.columns)]
    
    # Determine whether to create a legend
    make_legend = True
    if legend == None:
        make_legend = False
    
    # Add heat map to axes
    frame.plot(ax=axes, column=field, cmap=color, legend=make_legend,
               legend_kwds={"label": legend},
               missing_kwds={"color": "lightgray", "edgecolor": "gray",
               "hatch": "//////"})

#==============================================================================
# Statistics
#==============================================================================

def compare_two(file1, field1, file2, field2, names=None):
    """Comparison statistics for two data fields.
    
    Positional arguments:
        file1 (str) -- Path to first data file.
        field1 (str) -- Field of column in first data file
        file2 (str) -- Path to second data file.
        field2 (str) -- Field of column in second data file.
    
    Optional keyword arguments:
        names (tuple) -- Axis names. Default None. Otherwise, a tuple of (x, y)
            axis labels.
    
    Runs a collection of statistical tests and generates a collection of
    statistical graphics for comparing two data fields. This includes scatter
    plots, correlation coefficients, and R^2 values.
    """
    
    # Get both fields as lists
    x = _get_field(file1, field1)
    y = _get_field(file2, field2)
    
    print(f"\nComparing:\n{field1} ({file1})\n{field2} ({file2})")
    
    # Linear regression statistics
    m, b, corr, _, _ = stat.linregress(x, y)
    print("Linear regression statistics:")
    print(f"\tSlope:     {m}")
    print(f"\t1/Slope:   {1/m}")
    print(f"\tIntercept: {b}")
    print(f"\tCorr:      {corr}")
    
    # Scatter plot
    fig, ax = plt.subplots()
    plt.scatter(x, y)
    xleft, xright = ax.get_xlim()
    ybottom, ytop = ax.get_ylim()
    ax.set_aspect(abs((xright-xleft)/(ybottom-ytop)))
    #plt.xticks([], [])
    #plt.yticks([], [])
    if names != None:
        ax.set_xlabel(names[0])
        ax.set_ylabel(names[1])
    plt.show()

#==============================================================================

# Download Polk County files
#download_shapefiles("FL", "Polk", cfile=SHP_FL_COUNTIES, tfile=SHP_POLK_TRACTS, bfile=SHP_POLK_BLOCKS)

# Completed result files
RESULTS_PHARM = os.path.join(POLK_RESULTS, "polk_pop_pharm_results.tsv")
RESULTS_UC = os.path.join(POLK_RESULTS, "polk_pop_uc_results.tsv")

## Plot Polk County tracts
#fig, ax = plt.subplots()
#map_shapefile(SHP_POLK_TRACTS, ax)
#plt.xlim(POLK_X_FULL)
#plt.ylim(POLK_Y_FULL)
#ax.add_artist(scb.ScaleBar(POLK_DEGREE))
##plt.xticks([], [])
##plt.yticks([], [])
#ax.set_xlabel("Longitude")
#ax.set_ylabel("Latitude")
#plt.show()

## Plot Polk County blocks
#fig, ax = plt.subplots()
#map_shapefile(SHP_POLK_BLOCKS, ax)
#plt.xlim(POLK_X_FULL)
#plt.ylim(POLK_Y_FULL)
#ax.add_artist(scb.ScaleBar(POLK_DEGREE))
##plt.xticks([], [])
##plt.yticks([], [])
#ax.set_xlabel("Longitude")
#ax.set_ylabel("Latitude")
#plt.show()

## Plot Polk County with pharmacies
#fig, ax = plt.subplots()
##map_county(SHP_FL_COUNTIES, ax, "105")
#map_shapefile(SHP_POLK_TRACTS, ax)
##map_rectangle(LAKELAND_X, LAKELAND_Y, ax)
#map_points(os.path.join(POLK_PROCESSED, "polk_pharmacy.tsv"), ax, color="red")
#map_points(os.path.join(POLK_PROCESSED, "polk_pharmacy_nbr.tsv"), ax, color="red")
#plt.xlim(POLK_X_FULL)
#plt.ylim(POLK_Y_FULL)
#ax.add_artist(scb.ScaleBar(POLK_DEGREE))
##plt.xticks([], [])
##plt.yticks([], [])
#ax.set_xlabel("Longitude")
#ax.set_ylabel("Latitude")
#plt.show()

## Plot Lakeland detail with pharmacies
#fig, ax = plt.subplots()
#map_shapefile(SHP_POLK_TRACTS, ax)
#map_points(os.path.join(POLK_PROCESSED, "polk_pharmacy.tsv"), ax, color="red")
#map_points(os.path.join(POLK_PROCESSED, "polk_pharmacy_nbr.tsv"), ax, color="red")
#plt.xlim(LAKELAND_X)
#plt.ylim(LAKELAND_Y)
#ax.add_artist(scb.ScaleBar(POLK_DEGREE))
##plt.xticks([], [])
##plt.yticks([], [])
#ax.set_xlabel("Longitude")
#ax.set_ylabel("Latitude")
#plt.show()

## Plot Polk County with urgent care
#fig, ax = plt.subplots()
##map_county(SHP_FL_COUNTIES, ax, "105")
#map_shapefile(SHP_POLK_TRACTS, ax)
##map_rectangle(LAKELAND_X, LAKELAND_Y, ax)
#map_points(os.path.join(POLK_PROCESSED, "polk_uc.tsv"), ax, color="blue")
#map_points(os.path.join(POLK_PROCESSED, "polk_uc_nbr.tsv"), ax, color="blue")
#plt.xlim(POLK_X_FULL)
#plt.ylim(POLK_Y_FULL)
#ax.add_artist(scb.ScaleBar(POLK_DEGREE))
##plt.xticks([], [])
##plt.yticks([], [])
#ax.set_xlabel("Longitude")
#ax.set_ylabel("Latitude")
#plt.show()

## Plot Lakeland detail with urgent care
#fig, ax = plt.subplots()
#map_shapefile(SHP_POLK_TRACTS, ax)
#map_points(os.path.join(POLK_PROCESSED, "polk_uc.tsv"), ax, color="blue")
#map_points(os.path.join(POLK_PROCESSED, "polk_uc_nbr.tsv"), ax, color="blue")
#plt.xlim(LAKELAND_X)
#plt.ylim(LAKELAND_Y)
#ax.add_artist(scb.ScaleBar(POLK_DEGREE))
##plt.xticks([], [])
##plt.yticks([], [])
#ax.set_xlabel("Longitude")
#ax.set_ylabel("Latitude")
#plt.show()

## Plot Polk County population by tract
#fig, ax = plt.subplots()
#map_heat(SHP_POLK_TRACTS, RESULTS_PHARM, "pop", ax, color="Greens", legend="Population")
#map_county(SHP_FL_COUNTIES, ax, "105")
#plt.xlim(POLK_X_FULL)
#plt.ylim(POLK_Y_FULL)
#ax.add_artist(scb.ScaleBar(POLK_DEGREE))
##plt.xticks([], [])
##plt.yticks([], [])
#ax.set_xlabel("Longitude")
#ax.set_ylabel("Latitude")
#plt.show()

## Plot Polk County SVI
#fig, ax = plt.subplots()
#map_heat(SHP_POLK_TRACTS, RESULTS_PHARM, "svi", ax, color="Reds", legend="SVI")
#map_county(SHP_FL_COUNTIES, ax, "105")
#plt.xlim(POLK_X_FULL)
#plt.ylim(POLK_Y_FULL)
#ax.add_artist(scb.ScaleBar(POLK_DEGREE))
##plt.xticks([], [])
##plt.yticks([], [])
#ax.set_xlabel("Longitude")
#ax.set_ylabel("Latitude")
#plt.show()

## Plot Polk County pharmacies within 15 minutes
#fig, ax = plt.subplots()
#map_heat(SHP_POLK_TRACTS, RESULTS_PHARM, "fac-count_all-times_cutoff-15", ax, color="Blues", legend="Pharmacies within 15 minutes", limits=(0, 70))
#map_county(SHP_FL_COUNTIES, ax, "105")
#plt.xlim(POLK_X_FULL)
#plt.ylim(POLK_Y_FULL)
#ax.add_artist(scb.ScaleBar(POLK_DEGREE))
##plt.xticks([], [])
##plt.yticks([], [])
#ax.set_xlabel("Longitude")
#ax.set_ylabel("Latitude")
#plt.show()

## Plot Polk County pharmacies within 30 minutes
#fig, ax = plt.subplots()
#map_heat(SHP_POLK_TRACTS, RESULTS_PHARM, "fac-count_all-times_cutoff-30", ax, color="Blues", legend="Pharmacies within 30 minutes", limits=(0, 250))
#map_county(SHP_FL_COUNTIES, ax, "105")
#plt.xlim(POLK_X_FULL)
#plt.ylim(POLK_Y_FULL)
#ax.add_artist(scb.ScaleBar(POLK_DEGREE))
##plt.xticks([], [])
##plt.yticks([], [])
#ax.set_xlabel("Longitude")
#ax.set_ylabel("Latitude")
#plt.show()

## Plot Polk County urgent care within 15 minutes
#fig, ax = plt.subplots()
#map_heat(SHP_POLK_TRACTS, RESULTS_UC, "fac-count_all-times_cutoff-15", ax, color="Blues", legend="Urgent care within 15 minutes", limits=(0, 70))
#map_county(SHP_FL_COUNTIES, ax, "105")
#plt.xlim(POLK_X_FULL)
#plt.ylim(POLK_Y_FULL)
#ax.add_artist(scb.ScaleBar(POLK_DEGREE))
##plt.xticks([], [])
##plt.yticks([], [])
#ax.set_xlabel("Longitude")
#ax.set_ylabel("Latitude")
#plt.show()

## Plot Polk County urgent care within 30 minutes
#fig, ax = plt.subplots()
#map_heat(SHP_POLK_TRACTS, RESULTS_UC, "fac-count_all-times_cutoff-30", ax, color="Blues", legend="Urgent care within 30 minutes", limits=(0, 250))
#map_county(SHP_FL_COUNTIES, ax, "105")
#plt.xlim(POLK_X_FULL)
#plt.ylim(POLK_Y_FULL)
#ax.add_artist(scb.ScaleBar(POLK_DEGREE))
##plt.xticks([], [])
##plt.yticks([], [])
#ax.set_xlabel("Longitude")
#ax.set_ylabel("Latitude")
#plt.show()

## Plot Polk County average 5:00-10:00pm pharmacies within 15 minutes
#fig, ax = plt.subplots()
#map_heat(SHP_POLK_TRACTS, RESULTS_PHARM, "fac-count-avg_Wed_17:00-Wed_22:00_cutoff-15", ax, color="Blues", legend="Average pharmacies within 15 minutes during 5-10pm")
#map_county(SHP_FL_COUNTIES, ax, "105")
#plt.xlim(POLK_X_FULL)
#plt.ylim(POLK_Y_FULL)
#ax.add_artist(scb.ScaleBar(POLK_DEGREE))
##plt.xticks([], [])
##plt.yticks([], [])
#ax.set_xlabel("Longitude")
#ax.set_ylabel("Latitude")
#plt.show()

## Plot Polk County average 5:00-10:00pm pharmacies within 30 minutes
#fig, ax = plt.subplots()
#map_heat(SHP_POLK_TRACTS, RESULTS_PHARM, "fac-count-avg_Wed_17:00-Wed_22:00_cutoff-30", ax, color="Blues", legend="Average pharmacies within 30 minutes during 5-10pm")
#map_county(SHP_FL_COUNTIES, ax, "105")
#plt.xlim(POLK_X_FULL)
#plt.ylim(POLK_Y_FULL)
#ax.add_artist(scb.ScaleBar(POLK_DEGREE))
##plt.xticks([], [])
##plt.yticks([], [])
#ax.set_xlabel("Longitude")
#ax.set_ylabel("Latitude")
#plt.show()

## Plot Polk County average 5:00-10:00pm urgent care within 15 minutes
#fig, ax = plt.subplots()
#map_heat(SHP_POLK_TRACTS, RESULTS_UC, "fac-count-avg_Wed_17:00-Wed_22:00_cutoff-15", ax, color="Blues", legend="Average urgent care within 15 minutes during 5-10pm")
#map_county(SHP_FL_COUNTIES, ax, "105")
#plt.xlim(POLK_X_FULL)
#plt.ylim(POLK_Y_FULL)
#ax.add_artist(scb.ScaleBar(POLK_DEGREE))
##plt.xticks([], [])
##plt.yticks([], [])
#ax.set_xlabel("Longitude")
#ax.set_ylabel("Latitude")
#plt.show()

## Plot Polk County average 5:00-10:00pm urgent care within 30 minutes
#fig, ax = plt.subplots()
#map_heat(SHP_POLK_TRACTS, RESULTS_UC, "fac-count-avg_Wed_17:00-Wed_22:00_cutoff-30", ax, color="Blues", legend="Average urgent care within 30 minutes during 5-10pm")
#map_county(SHP_FL_COUNTIES, ax, "105")
#plt.xlim(POLK_X_FULL)
#plt.ylim(POLK_Y_FULL)
#ax.add_artist(scb.ScaleBar(POLK_DEGREE))
##plt.xticks([], [])
##plt.yticks([], [])
#ax.set_xlabel("Longitude")
#ax.set_ylabel("Latitude")
#plt.show()

## Plot Polk County indicator of no pharmacies within 15 minutes
#fig, ax = plt.subplots()
#map_heat(SHP_POLK_TRACTS, RESULTS_PHARM, "fac-count-below-1_all-times_cutoff-15", ax, color=white_red)
#map_county(SHP_FL_COUNTIES, ax, "105")
#map_shapefile(SHP_POLK_TRACTS, ax)
#plt.xlim(POLK_X_FULL)
#plt.ylim(POLK_Y_FULL)
#ax.add_artist(scb.ScaleBar(POLK_DEGREE))
##plt.xticks([], [])
##plt.yticks([], [])
#ax.set_xlabel("Longitude")
#ax.set_ylabel("Latitude")
#plt.show()

## Plot Polk County indicator of no pharmacies within 30 minutes
#fig, ax = plt.subplots()
#map_heat(SHP_POLK_TRACTS, RESULTS_PHARM, "fac-count-below-1_all-times_cutoff-30", ax, color=white_red)
#map_county(SHP_FL_COUNTIES, ax, "105")
#map_shapefile(SHP_POLK_TRACTS, ax)
#plt.xlim(POLK_X_FULL)
#plt.ylim(POLK_Y_FULL)
#ax.add_artist(scb.ScaleBar(POLK_DEGREE))
##plt.xticks([], [])
##plt.yticks([], [])
#ax.set_xlabel("Longitude")
#ax.set_ylabel("Latitude")
#plt.show()

## Plot Polk County indicator of no urgent cares within 15 minutes
#fig, ax = plt.subplots()
#map_heat(SHP_POLK_TRACTS, RESULTS_UC, "fac-count-below-1_all-times_cutoff-15", ax, color=white_red)
#map_county(SHP_FL_COUNTIES, ax, "105")
#map_shapefile(SHP_POLK_TRACTS, ax)
#plt.xlim(POLK_X_FULL)
#plt.ylim(POLK_Y_FULL)
#ax.add_artist(scb.ScaleBar(POLK_DEGREE))
##plt.xticks([], [])
##plt.yticks([], [])
#ax.set_xlabel("Longitude")
#ax.set_ylabel("Latitude")
#plt.show()

## Plot Polk County indicator of no urgent cares within 30 minutes
#fig, ax = plt.subplots()
#map_heat(SHP_POLK_TRACTS, RESULTS_UC, "fac-count-below-1_all-times_cutoff-30", ax, color=white_red)
#map_county(SHP_FL_COUNTIES, ax, "105")
#map_shapefile(SHP_POLK_TRACTS, ax)
#plt.xlim(POLK_X_FULL)
#plt.ylim(POLK_Y_FULL)
#ax.add_artist(scb.ScaleBar(POLK_DEGREE))
##plt.xticks([], [])
##plt.yticks([], [])
#ax.set_xlabel("Longitude")
#ax.set_ylabel("Latitude")
#plt.show()

## Plot Polk County fraction of 5-10pm with access to at least 1 pharmacy in 15 minutes
#fig, ax = plt.subplots()
#map_heat(SHP_POLK_TRACTS, RESULTS_PHARM, "frac-time-above-1_Wed_17:00-Wed_22:00_cutoff-15", ax, color="Blues", legend="Frac of 5-10pm w/ access to 1 pharmacy in 15 min", limits=(0.0, 1.0))
#map_county(SHP_FL_COUNTIES, ax, "105")
#plt.xlim(POLK_X_FULL)
#plt.ylim(POLK_Y_FULL)
#ax.add_artist(scb.ScaleBar(POLK_DEGREE))
##plt.xticks([], [])
##plt.yticks([], [])
#ax.set_xlabel("Longitude")
#ax.set_ylabel("Latitude")
#plt.show()

## Plot Polk County fraction of 5-10pm with access to at least 1 pharmacy in 30 minutes
#fig, ax = plt.subplots()
#map_heat(SHP_POLK_TRACTS, RESULTS_PHARM, "frac-time-above-1_Wed_17:00-Wed_22:00_cutoff-30", ax, color="Blues", legend="Frac of 5-10pm w/ access to 1 pharmacy in 30 min", limits=(0.0, 1.0))
#map_county(SHP_FL_COUNTIES, ax, "105")
#plt.xlim(POLK_X_FULL)
#plt.ylim(POLK_Y_FULL)
#ax.add_artist(scb.ScaleBar(POLK_DEGREE))
##plt.xticks([], [])
##plt.yticks([], [])
#ax.set_xlabel("Longitude")
#ax.set_ylabel("Latitude")
#plt.show()

## Plot Polk County fraction of 5-10pm with access to at least 1 urgent care in 15 minutes
#fig, ax = plt.subplots()
#map_heat(SHP_POLK_TRACTS, RESULTS_UC, "frac-time-above-1_Wed_17:00-Wed_22:00_cutoff-15", ax, color="Blues", legend="Frac of 5-10pm w/ access to 1 urgent care in 15 min", limits=(0.0, 1.0))
#map_county(SHP_FL_COUNTIES, ax, "105")
#plt.xlim(POLK_X_FULL)
#plt.ylim(POLK_Y_FULL)
#ax.add_artist(scb.ScaleBar(POLK_DEGREE))
##plt.xticks([], [])
##plt.yticks([], [])
#ax.set_xlabel("Longitude")
#ax.set_ylabel("Latitude")
#plt.show()

## Plot Polk County fraction of 5-10pm with access to at least 1 urgent care in 30 minutes
#fig, ax = plt.subplots()
#map_heat(SHP_POLK_TRACTS, RESULTS_UC, "frac-time-above-1_Wed_17:00-Wed_22:00_cutoff-30", ax, color="Blues", legend="Frac of 5-10pm w/ access to 1 urgent care in 30 min", limits=(0.0, 1.0))
#map_county(SHP_FL_COUNTIES, ax, "105")
#plt.xlim(POLK_X_FULL)
#plt.ylim(POLK_Y_FULL)
#ax.add_artist(scb.ScaleBar(POLK_DEGREE))
##plt.xticks([], [])
##plt.yticks([], [])
#ax.set_xlabel("Longitude")
#ax.set_ylabel("Latitude")
#plt.show()

## Plot Polk County fraction of 5-10pm with access to at least 3 pharmacies in 15 minutes
#fig, ax = plt.subplots()
#map_heat(SHP_POLK_TRACTS, RESULTS_PHARM, "frac-time-above-3_Wed_17:00-Wed_22:00_cutoff-15", ax, color="Blues", legend="Frac of 5-10pm w/ access to 3 pharmacies in 15 min", limits=(0.0, 1.0))
#map_county(SHP_FL_COUNTIES, ax, "105")
#plt.xlim(POLK_X_FULL)
#plt.ylim(POLK_Y_FULL)
#ax.add_artist(scb.ScaleBar(POLK_DEGREE))
##plt.xticks([], [])
##plt.yticks([], [])
#ax.set_xlabel("Longitude")
#ax.set_ylabel("Latitude")
#plt.show()

## Plot Polk County fraction of 5-10pm with access to at least 3 pharmacies in 30 minutes
#fig, ax = plt.subplots()
#map_heat(SHP_POLK_TRACTS, RESULTS_PHARM, "frac-time-above-3_Wed_17:00-Wed_22:00_cutoff-30", ax, color="Blues", legend="Frac of 5-10pm w/ access to 3 pharmacies in 30 min", limits=(0.0, 1.0))
#map_county(SHP_FL_COUNTIES, ax, "105")
#plt.xlim(POLK_X_FULL)
#plt.ylim(POLK_Y_FULL)
#ax.add_artist(scb.ScaleBar(POLK_DEGREE))
##plt.xticks([], [])
##plt.yticks([], [])
#ax.set_xlabel("Longitude")
#ax.set_ylabel("Latitude")
#plt.show()

## Plot Polk County fraction of 5-10pm with access to at least 3 urgent cares in 15 minutes
#fig, ax = plt.subplots()
#map_heat(SHP_POLK_TRACTS, RESULTS_UC, "frac-time-above-3_Wed_17:00-Wed_22:00_cutoff-15", ax, color="Blues", legend="Frac of 5-10pm w/ access to 3 urgent cares in 15 min", limits=(0.0, 1.0))
#map_county(SHP_FL_COUNTIES, ax, "105")
#plt.xlim(POLK_X_FULL)
#plt.ylim(POLK_Y_FULL)
#ax.add_artist(scb.ScaleBar(POLK_DEGREE))
##plt.xticks([], [])
##plt.yticks([], [])
#ax.set_xlabel("Longitude")
#ax.set_ylabel("Latitude")
#plt.show()

## Plot Polk County fraction of 5-10pm with access to at least 3 urgent cares in 30 minutes
#fig, ax = plt.subplots()
#map_heat(SHP_POLK_TRACTS, RESULTS_UC, "frac-time-above-3_Wed_17:00-Wed_22:00_cutoff-30", ax, color="Blues", legend="Frac of 5-10pm w/ access to 3 urgent cares in 30 min", limits=(0.0, 1.0))
#map_county(SHP_FL_COUNTIES, ax, "105")
#plt.xlim(POLK_X_FULL)
#plt.ylim(POLK_Y_FULL)
#ax.add_artist(scb.ScaleBar(POLK_DEGREE))
##plt.xticks([], [])
##plt.yticks([], [])
#ax.set_xlabel("Longitude")
#ax.set_ylabel("Latitude")
#plt.show()

## Plot Polk County fraction of 5-10pm with access to at least 5 pharmacies in 15 minutes
#fig, ax = plt.subplots()
#map_heat(SHP_POLK_TRACTS, RESULTS_PHARM, "frac-time-above-5_Wed_17:00-Wed_22:00_cutoff-15", ax, color="Blues", legend="Frac of 5-10pm w/ access to 5 pharmacies in 15 min", limits=(0.0, 1.0))
#map_county(SHP_FL_COUNTIES, ax, "105")
#plt.xlim(POLK_X_FULL)
#plt.ylim(POLK_Y_FULL)
#ax.add_artist(scb.ScaleBar(POLK_DEGREE))
##plt.xticks([], [])
##plt.yticks([], [])
#ax.set_xlabel("Longitude")
#ax.set_ylabel("Latitude")
#plt.show()

## Plot Polk County fraction of 5-10pm with access to at least 5 pharmacies in 30 minutes
#fig, ax = plt.subplots()
#map_heat(SHP_POLK_TRACTS, RESULTS_PHARM, "frac-time-above-5_Wed_17:00-Wed_22:00_cutoff-30", ax, color="Blues", legend="Frac of 5-10pm w/ access to 5 pharmacies in 30 min", limits=(0.0, 1.0))
#map_county(SHP_FL_COUNTIES, ax, "105")
#plt.xlim(POLK_X_FULL)
#plt.ylim(POLK_Y_FULL)
#ax.add_artist(scb.ScaleBar(POLK_DEGREE))
##plt.xticks([], [])
##plt.yticks([], [])
#ax.set_xlabel("Longitude")
#ax.set_ylabel("Latitude")
#plt.show()

## Plot Polk County fraction of 5-10pm with access to at least 5 urgent cares in 15 minutes
#fig, ax = plt.subplots()
#map_heat(SHP_POLK_TRACTS, RESULTS_UC, "frac-time-above-5_Wed_17:00-Wed_22:00_cutoff-15", ax, color="Blues", legend="Frac of 5-10pm w/ access to 5 urgent cares in 15 min", limits=(0.0, 1.0))
#map_county(SHP_FL_COUNTIES, ax, "105")
#plt.xlim(POLK_X_FULL)
#plt.ylim(POLK_Y_FULL)
#ax.add_artist(scb.ScaleBar(POLK_DEGREE))
##plt.xticks([], [])
##plt.yticks([], [])
#ax.set_xlabel("Longitude")
#ax.set_ylabel("Latitude")
#plt.show()

## Plot Polk County fraction of 5-10pm with access to at least 5 urgent cares in 30 minutes
#fig, ax = plt.subplots()
#map_heat(SHP_POLK_TRACTS, RESULTS_UC, "frac-time-above-5_Wed_17:00-Wed_22:00_cutoff-30", ax, color="Blues", legend="Frac of 5-10pm w/ access to 5 urgent cares in 30 min", limits=(0.0, 1.0))
#map_county(SHP_FL_COUNTIES, ax, "105")
#plt.xlim(POLK_X_FULL)
#plt.ylim(POLK_Y_FULL)
#ax.add_artist(scb.ScaleBar(POLK_DEGREE))
##plt.xticks([], [])
##plt.yticks([], [])
#ax.set_xlabel("Longitude")
#ax.set_ylabel("Latitude")
#plt.show()

## Pharmacy and urgent care comparisons
#compare_two(RESULTS_PHARM, "fac-count_all-times_cutoff-15", RESULTS_UC, "fac-count_all-times_cutoff-15", names=("Pharmacies within 15 minutes", "Urgent care within 15 minutes"))
#compare_two(RESULTS_PHARM, "fac-count_all-times_cutoff-30", RESULTS_UC, "fac-count_all-times_cutoff-30", names=("Pharmacies within 30 minutes", "Urgent care within 30 minutes"))
#compare_two(RESULTS_PHARM, "fac-count-avg_Wed_17:00-Wed_22:00_cutoff-15", RESULTS_UC, "fac-count-avg_Wed_17:00-Wed_22:00_cutoff-15", names=("Average pharmacies within 15 minutes during 5-10pm", "Average urgent cares within 15 minutes during 5-10pm"))
#compare_two(RESULTS_PHARM, "fac-count-avg_Wed_17:00-Wed_22:00_cutoff-30", RESULTS_UC, "fac-count-avg_Wed_17:00-Wed_22:00_cutoff-30", names=("Average pharmacies within 30 minutes during 5-10pm", "Average urgent cares within 30 minutes during 5-10pm"))
##compare_two(RESULTS_PHARM, "frac-time-above-1_Wed_17:00-Wed_22:00_cutoff-15", RESULTS_UC, "frac-time-above-1_Wed_17:00-Wed_22:00_cutoff-15", names=("Frac of 5-10pm w/ access to 1 pharmacy within 15 min", "Frac of 5-10pm w/ access to 1 urgent care within 15 min"))
##compare_two(RESULTS_PHARM, "frac-time-above-1_Wed_17:00-Wed_22:00_cutoff-30", RESULTS_UC, "frac-time-above-1_Wed_17:00-Wed_22:00_cutoff-30", names=("Frac of 5-10pm w/ access to 1 pharmacy within 30 min", "Frac of 5-10pm w/ access to 1 urgent care within 30 min"))
##compare_two(RESULTS_PHARM, "frac-time-above-3_Wed_17:00-Wed_22:00_cutoff-15", RESULTS_UC, "frac-time-above-3_Wed_17:00-Wed_22:00_cutoff-15", names=("Frac of 5-10pm w/ access to 3 pharmacies within 15 min", "Frac of 5-10pm w/ access to 3 urgent cares within 15 min"))
##compare_two(RESULTS_PHARM, "frac-time-above-3_Wed_17:00-Wed_22:00_cutoff-30", RESULTS_UC, "frac-time-above-3_Wed_17:00-Wed_22:00_cutoff-30", names=("Frac of 5-10pm w/ access to 3 pharmacies within 30 min", "Frac of 5-10pm w/ access to 3 urgent cares within 30 min"))
##compare_two(RESULTS_PHARM, "frac-time-above-5_Wed_17:00-Wed_22:00_cutoff-15", RESULTS_UC, "frac-time-above-5_Wed_17:00-Wed_22:00_cutoff-15", names=("Frac of 5-10pm w/ access to 5 pharmacies within 15 min", "Frac of 5-10pm w/ access to 5 urgent cares within 15 min"))
##compare_two(RESULTS_PHARM, "frac-time-above-5_Wed_17:00-Wed_22:00_cutoff-30", RESULTS_UC, "frac-time-above-5_Wed_17:00-Wed_22:00_cutoff-30", names=("Frac of 5-10pm w/ access to 5 pharmacies within 30 min", "Frac of 5-10pm w/ access to 5 urgent cares within 30 min"))

## SVI comparisons
#compare_two(RESULTS_PHARM, "svi", RESULTS_PHARM, "fac-count_all-times_cutoff-15", names=("SVI", "Pharmacies within 15 minutes"))
#compare_two(RESULTS_PHARM, "svi", RESULTS_PHARM, "fac-count_all-times_cutoff-30", names=("SVI", "Pharmacies within 30 minutes"))
#compare_two(RESULTS_PHARM, "svi", RESULTS_UC, "fac-count_all-times_cutoff-15", names=("SVI", "Urgent care within 15 minutes"))
#compare_two(RESULTS_PHARM, "svi", RESULTS_UC, "fac-count_all-times_cutoff-30", names=("SVI", "Urgent care within 30 minutes"))
#compare_two(RESULTS_PHARM, "svi", RESULTS_PHARM, "fac-count-avg_Wed_17:00-Wed_22:00_cutoff-15", names=("SVI", "Average pharmacies within 15 minutes during 5-10pm"))
#compare_two(RESULTS_PHARM, "svi", RESULTS_PHARM, "fac-count-avg_Wed_17:00-Wed_22:00_cutoff-30", names=("SVI", "Average pharmacies within 30 minutes during 5-10pm"))
#compare_two(RESULTS_PHARM, "svi", RESULTS_UC, "fac-count-avg_Wed_17:00-Wed_22:00_cutoff-15", names=("SVI", "Average urgent cares within 15 minutes during 5-10pm"))
#compare_two(RESULTS_PHARM, "svi", RESULTS_UC, "fac-count-avg_Wed_17:00-Wed_22:00_cutoff-30", names=("SVI", "Average urgent cares within 30 minutes during 5-10pm"))
##compare_two(RESULTS_PHARM, "svi", RESULTS_PHARM, "frac-time-above-1_Wed_17:00-Wed_22:00_cutoff-15", names=("SVI", "Frac of 5-10pm w/ access to 1 pharmacy within 15 min"))
##compare_two(RESULTS_PHARM, "svi", RESULTS_PHARM, "frac-time-above-1_Wed_17:00-Wed_22:00_cutoff-30", names=("SVI", "Frac of 5-10pm w/ access to 1 pharmacy within 30 min"))
##compare_two(RESULTS_PHARM, "svi", RESULTS_PHARM, "frac-time-above-3_Wed_17:00-Wed_22:00_cutoff-15", names=("SVI", "Frac of 5-10pm w/ access to 3 pharmacies within 15 min"))
##compare_two(RESULTS_PHARM, "svi", RESULTS_PHARM, "frac-time-above-3_Wed_17:00-Wed_22:00_cutoff-30", names=("SVI", "Frac of 5-10pm w/ access to 3 pharmacies within 30 min"))
##compare_two(RESULTS_PHARM, "svi", RESULTS_PHARM, "frac-time-above-5_Wed_17:00-Wed_22:00_cutoff-15", names=("SVI", "Frac of 5-10pm w/ access to 5 pharmacies within 15 min"))
##compare_two(RESULTS_PHARM, "svi", RESULTS_PHARM, "frac-time-above-5_Wed_17:00-Wed_22:00_cutoff-30", names=("SVI", "Frac of 5-10pm w/ access to 5 pharmacies within 30 min"))
##compare_two(RESULTS_PHARM, "svi", RESULTS_UC, "frac-time-above-1_Wed_17:00-Wed_22:00_cutoff-15", names=("SVI", "Frac of 5-10pm w/ access to 1 urgent care within 15 min"))
##compare_two(RESULTS_PHARM, "svi", RESULTS_UC, "frac-time-above-1_Wed_17:00-Wed_22:00_cutoff-30", names=("SVI", "Frac of 5-10pm w/ access to 1 urgent care within 30 min"))
##compare_two(RESULTS_PHARM, "svi", RESULTS_UC, "frac-time-above-3_Wed_17:00-Wed_22:00_cutoff-15", names=("SVI", "Frac of 5-10pm w/ access to 3 urgent cares within 15 min"))
##compare_two(RESULTS_PHARM, "svi", RESULTS_UC, "frac-time-above-3_Wed_17:00-Wed_22:00_cutoff-30", names=("SVI", "Frac of 5-10pm w/ access to 3 urgent cares within 30 min"))
##compare_two(RESULTS_PHARM, "svi", RESULTS_UC, "frac-time-above-5_Wed_17:00-Wed_22:00_cutoff-15", names=("SVI", "Frac of 5-10pm w/ access to 5 urgent cares within 15 min"))
##compare_two(RESULTS_PHARM, "svi", RESULTS_UC, "frac-time-above-5_Wed_17:00-Wed_22:00_cutoff-30", names=("SVI", "Frac of 5-10pm w/ access to 5 urgent cares within 30 min"))
##compare_two(RESULTS_PHARM, "svi", RESULTS_UC, "pov150", names=("SVI", "Fraction of population below 150% of poverty line"))
##compare_two(RESULTS_PHARM, "svi", RESULTS_UC, "noveh", names=("SVI", "Fraction of population with no vehicle access"))

## Population/accessibility comparisons
#compare_two(RESULTS_PHARM, "pop", RESULTS_PHARM, "fac-count_all-times_cutoff-15", names=("Population", "Pharmacies within 15 minutes"))
#compare_two(RESULTS_PHARM, "pop", RESULTS_PHARM, "fac-count_all-times_cutoff-30", names=("Population", "Pharmacies within 30 minutes"))
#compare_two(RESULTS_PHARM, "pop", RESULTS_UC, "fac-count_all-times_cutoff-15", names=("Population", "Urgent care within 15 minutes"))
#compare_two(RESULTS_PHARM, "pop", RESULTS_UC, "fac-count_all-times_cutoff-30", names=("Population", "Urgent care within 30 minutes"))
