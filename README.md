# Pharmacy Accessibility Study

![Python 3.8](https://img.shields.io/badge/python-3%2E8-blue) ![MIT License](https://img.shields.io/github/license/adam-rumpf/pharmacy-access)

A collection of data processing scripts for use in an ongoing research project about measuring accessibility levels to pharmacies. These scripts are unlikely to be of use to anyone outside our research group, but are provided here for anyone interested. Python packages used can be found in the root-level `requirements` file.

## Contents

This repo is organized into several major subdirectories, each with its own README to explain its contents in more detail.

* `data/` is divided into subdirectories by location, and includes raw data files obtained from external sources. The data files, themselves, have been excluded from this repo, but their sources can be found in this subdirectory's README.
* `graphs/` is divided into subdirectories by location, and includes simplified graph representations of raw map files for use in distance computations.
* `maps/` is divided into subdirectories by location, and includes raw map files downloaded for use in travel time calculation.
* `processed/` is divided into subdirectories by location, and includes data files derived from the raw data, such as collated data tables and calculated travel time matrices.
* `results/` is divided into subdirectories by location, and includes the output files from the statistical analysis.
* `scripts/` contains the programs actually used to preprocess the data and to conduct a statistical analysis.
* `shapefiles/` contains shapefiles for use in generating graphics.

## Files

The raw data files, themselves, are not included in this repo, but their sources are indicated in the various subdirectory READMEs. Some of the scripts also require the use of a user email address or API token. These are meant to be stored in local files called `email.txt` or `token.txt`, which are also not included in this repo.

Most of the directories contain scripts and subdirectories that mention Santa Clara County, CA. These are left over from a previous iteration of the project, before the focus was shifted to Polk County, FL.
