"""Pharmacy accessibility project graphics generation scripts.

The code below is meant for generating graphics for use in data analysis and
result display. This includes both maps to display geospatial distributions as
well as plots for statistical results.
"""

import csv
import os.path

import matplotlib.pyplot as plt
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
POLK_RESULTS = os.path.join("..", "results", "polk")

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

def make_geodata(shapefile, datafile, fields, shapeid="GEOID", dataid="name",
                 missing=None):
    """Adds additional data fields to a GeoDataFrame.
    
    Positional arguments:
        shapefile (str) -- File path for the shapefile (tract, block, etc.) to
            be used as the basis of the GeoDataFrame.
        datafile (str) -- File path for a data file containing fields to be
            added to the shapefile's fields.
        fields (list(str)) -- List of fields in the data file to be added to
            the GeoDataFrame.
    
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
    """
    
    pass

#==============================================================================
# Map Generation
#==============================================================================

### Useful colormap reference: https://matplotlib.org/stable/users/explain/colors/colormaps.html

### To-dos:
    # Plain tract map with population centroids
    # Plain tract map with pharmacy locations
    # Plain tract map with urgent care locations
    # Tract heatmap with average number of facilities available during [[[time slot]]]

### Implement each of these as a map layer generator, then make a master script
### that combines several of them.

#==============================================================================

# Download Polk County files
#download_shapefiles("FL", "Polk", cfile=SHP_FL_COUNTIES, tfile=SHP_POLK_TRACTS, bfile=SHP_POLK_BLOCKS)
