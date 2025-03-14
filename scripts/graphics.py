"""Pharmacy accessibility project graphics generation scripts.

The code below is meant for generating graphics for use in data analysis and
result display. This includes both maps to display geospatial distributions as
well as plots for statistical results.
"""

import csv
import os.path

import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib_scalebar.scalebar as sb
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
POLK_X_FULL = [-82.155, -81.087]
POLK_Y_FULL = [27.608, 28.397]

# Distance of 1 degree near Polk County (meters)
POLK_DEGREE = 98358.86851331996

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

def make_geodata(shapefile, datafile, field, shapeid="GEOID", dataid="name",
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
    
    pass

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
    gpd.GeoDataFrame(row).plot(ax=axes, color="white", edgecolor=color)

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

def map_heat(shapefile, datafile, field, color="viridis", shapeid="GEOID",
             dataid="name", missing=None):
    """Adds a heat map based on a given field to a given set of axes.
    
    Positional arguments:
        shapefile (str) -- Path to the main shapefile.
        datafile (str) -- Path to the data file on which the heat map colors
            will be based.
        field (str) -- Name of field in data file to use for the heat map.
    
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
    """
    
    # Create a single unified GeoDataFrame
    frame = make_geodata(shapefile, datafile, field, shapeid=shapeid,
                         dataid=dataid, missing=missing)
    
    ###

#==============================================================================

# Download Polk County files
#download_shapefiles("FL", "Polk", cfile=SHP_FL_COUNTIES, tfile=SHP_POLK_TRACTS, bfile=SHP_POLK_BLOCKS)

# Plot Polk County with pharmacies
fig, ax = plt.subplots()
map_county(SHP_FL_COUNTIES, ax, "105")
map_points(os.path.join(POLK_PROCESSED, "polk_pharmacy.tsv"), ax, color="red")
map_points(os.path.join(POLK_PROCESSED, "polk_pharmacy_nbr.tsv"), ax, color="red")
plt.xlim(POLK_X_FULL)
plt.ylim(POLK_Y_FULL)
ax.add_artist(sb.ScaleBar(POLK_DEGREE))
plt.show()

# Add urgent care
fig, ax = plt.subplots()
map_county(SHP_FL_COUNTIES, ax, "105")
map_points(os.path.join(POLK_PROCESSED, "polk_uc.tsv"), ax, color="blue")
map_points(os.path.join(POLK_PROCESSED, "polk_uc_nbr.tsv"), ax, color="blue")
plt.xlim(POLK_X_FULL)
plt.ylim(POLK_Y_FULL)
ax.add_artist(sb.ScaleBar(POLK_DEGREE))
plt.show()
