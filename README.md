# COVID-19 Vaccine Accessibility Study

A collection of data processing scripts for use in an ongoing research project about measuring accessibility levels to COVID-19 vaccines. These scripts are unlikely to be of use to anyone ourside our research group, but are provided here for anyone interested. Python packages used can be found in the root-level `requirements` file.

## Contents

This repo is organized into several major subdirectories, each with its own README to explain its contents in more detail.

* `data/` is divided into subdirectories by location, and includes raw data files obtained from external sources. The data files, themselves, have been excluded from this repo, but their sources can be found in this subdirectory's README.
* `processed/` is divided into subdirectories by location, and includes data files derived from the raw data, such as collated data tables and calculated travel time matrices.
* `results/` is divied into subdirectories by location, and includes the output files from the statistical analysis.
* `scripts/` contains the programs actually used to preprocess the data and to conduct a statistical analysis.
