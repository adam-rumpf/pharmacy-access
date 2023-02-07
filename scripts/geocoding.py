"""COVID-19 Vaccine accessibility project geocoding scripts.

The code below includes various scripts for geocoding the addresses from
certain regions. This is meant to speed up the process of getting geographic
coordinates from the pharmacy addresses available online.
"""

import time

from tqdm import tqdm
from geopy.geocoders import Nominatim

#==============================================================================
# Common Functions
#==============================================================================

def address_to_coords(geocoder, address):
    """Computes the geographic coordinates of a given address.
    
    Positional arguments:
        geocoder (Nominatim) -- geocoder object
        address (str) -- address string
    
    Returns:
        (tuple(float)) -- (latitude, longitude) tuple
    """
    
    try:
        location = geocoder.geocode(address)
        return (location.latitude, location.longitude)
    except geopy.exc.GeocoderTimedOut:
        print("geocoder timed out on query: " + address)
        return (None, None)

#==============================================================================
# Location-Specific Processing Scripts
#==============================================================================

###
### Use tqdm for progress bars
### Santa Clara pharmacy address geocoding
### note: Nominatim requires no more than 1 request per second, time.sleep(1.0)
### geocoder = Nominatim(user_agent="geopy/2.3.0")

#==============================================================================
# Execution
#==============================================================================

# Comment or uncomment the function calls below to process each location.
###
