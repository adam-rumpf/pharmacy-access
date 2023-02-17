# Scripts

This subdirectory includes a variety of scripts for use in processing data files and generating statistical results. The specific scripts are described below.

* `preprocessing.py` is the main preprocessing script, used to collect data from the raw data files in the `data/` subdirectory and assemble them into standardized data tables.
* `metrics.py` is the main accessibility metric computing script. It uses preprocessed data from the `processed/` subdirectory.
* `travel_times.py` is used to generate travel time matrices between population centers and vaccine providers.

## Files

The following files should be included in this directory.

* `email.txt`: A text file containing an email address. This will be read and used as the user_agent field for the Nominatim geocoder API. A placeholder file, `_email.txt`, is included in this repo as a reminder.
