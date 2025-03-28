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
             dataid="name", missing=None, legend=None):
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
    """
    
    # Create a single unified GeoDataFrame
    frame = _make_geodata(shapefile, datafile, field, shapeid=shapeid,
                          dataid=dataid, missing=missing)
    
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

# Download Polk County files
#download_shapefiles("FL", "Polk", cfile=SHP_FL_COUNTIES, tfile=SHP_POLK_TRACTS, bfile=SHP_POLK_BLOCKS)

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

# Plot Polk County pharmacies within 15 minutes
fig, ax = plt.subplots()
map_heat(SHP_POLK_TRACTS, os.path.join(POLK_RESULTS, "polk_pop_pharm_count_15-cutoff_all-times.tsv"), "access", ax, color="Blues", legend="Pharmacies within 15 minutes")
map_county(SHP_FL_COUNTIES, ax, "105")
plt.xlim(POLK_X_FULL)
plt.ylim(POLK_Y_FULL)
ax.add_artist(scb.ScaleBar(POLK_DEGREE))
#plt.xticks([], [])
#plt.yticks([], [])
ax.set_xlabel("Longitude")
ax.set_ylabel("Latitude")
plt.show()

# Plot Polk County pharmacies within 30 minutes
fig, ax = plt.subplots()
map_heat(SHP_POLK_TRACTS, os.path.join(POLK_RESULTS, "polk_pop_pharm_count_30-cutoff_all-times.tsv"), "access", ax, color="Blues", legend="Pharmacies within 30 minutes")
map_county(SHP_FL_COUNTIES, ax, "105")
plt.xlim(POLK_X_FULL)
plt.ylim(POLK_Y_FULL)
ax.add_artist(scb.ScaleBar(POLK_DEGREE))
#plt.xticks([], [])
#plt.yticks([], [])
ax.set_xlabel("Longitude")
ax.set_ylabel("Latitude")
plt.show()

# Plot Polk County urgent care within 15 minutes
fig, ax = plt.subplots()
map_heat(SHP_POLK_TRACTS, os.path.join(POLK_RESULTS, "polk_pop_uc_count_15-cutoff_all-times.tsv"), "access", ax, color="Blues", legend="Urgent care within 15 minutes")
map_county(SHP_FL_COUNTIES, ax, "105")
plt.xlim(POLK_X_FULL)
plt.ylim(POLK_Y_FULL)
ax.add_artist(scb.ScaleBar(POLK_DEGREE))
#plt.xticks([], [])
#plt.yticks([], [])
ax.set_xlabel("Longitude")
ax.set_ylabel("Latitude")
plt.show()

# Plot Polk County urgent care within 30 minutes
fig, ax = plt.subplots()
map_heat(SHP_POLK_TRACTS, os.path.join(POLK_RESULTS, "polk_pop_uc_count_30-cutoff_all-times.tsv"), "access", ax, color="Blues", legend="Urgent care within 30 minutes")
map_county(SHP_FL_COUNTIES, ax, "105")
plt.xlim(POLK_X_FULL)
plt.ylim(POLK_Y_FULL)
ax.add_artist(scb.ScaleBar(POLK_DEGREE))
#plt.xticks([], [])
#plt.yticks([], [])
ax.set_xlabel("Longitude")
ax.set_ylabel("Latitude")
plt.show()
