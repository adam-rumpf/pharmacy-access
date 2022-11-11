# Preprocessing Scripts

The scripts in this folder are mostly meant for processing the raw data from assorted files into single, unified data sets for use elsewhere. The raw data files are not included in this repo, but are explained and cited below. The preprocessing scripts included in this folder refer to the specific files named below.

## Population File Format

A population file is generated for each location, summarizing population center-specific data in a standardized format. Each population file is a TSV file with the following columns:

* `id`: Index of the population center within this location. All population centers are given ascending integer indices starting with 0.
* `name`: The unique identifier given to this population center within its source file (usually FIPS or ZIP code).
* `lat`: Latitude of centroid.
* `lon`: Longitude of centroid.
* `pop`: Population, as of the 2020 census.
* `vacc`: Vaccination rate, as of the most recent data.
* `adi`: 2022 [Area Deprivation Index (ADI)](https://www.nejm.org/doi/full/10.1056/NEJMp1802313) ranking.

## Facility File Format

A facility file is generated for each location, summarizing vaccination facility-specific data in a standardized format. Each facility file is a TSV file with the following columns:

* `id`: Index of the facilities center within this location. All facilities are given ascending integer indices starting with 0.
* `name`: An identifying string based on its source file (address, coordinates, etc.).
* `lat`: Latitude of facility.
* `lon`: Longitude of facility.
* `cap`: Capacity of facility (for now just 1 for everything, but might include something to account for resource limitations in the future).

## Raw Data

The raw data files for each location, along with sources, are included below.

### Chicago

Data mostly collected from the [Chicago Data Portal](https://data.cityofchicago.org/) on November 7, 2022.

* [`IL_2020_ADI_9 Digit Zip Code_v3.2.csv`](https://www.neighborhoodatlas.medicine.wisc.edu/): ZIP code-level [Area Deprivation Index (ADI)](https://www.nejm.org/doi/full/10.1056/NEJMp1802313) rankings for Illinois, indexed by ZIP code.
* [`COVID-19_Cases__Test_and_Deaths_by_ZIP_Code.csv`](https://data.cityofchicago.org/Health-Human-Services/COVID-19-Cases-Tests-and-Deaths-by-ZIP-Code/yhhz-zm2v): Historical data for Chicago COVID-19 cases and deaths by ZIP code. Also includes centroids of ZIP codes and populations by ZIP code.
* [`COVID-19_Vaccination_Locations.csv`](https://data.cityofchicago.org/Health-Human-Services/COVID-19-Vaccination-Locations/6q3z-9maq): Geographic data for COVID-19 vaccine locations in Chicago.
* [`COVID-19_Vaccinations_by_ZIP_Code.csv`](https://data.cityofchicago.org/Health-Human-Services/COVID-19-Vaccinations-by-ZIP-Code/553k-3xzc): Historical data for Chicago COVID-19 vaccination rates by ZIP code. Also includes centroids of ZIP codes.

### Santa Clara County

Data mostly collected from the [County of Santa Clara Open Data Portal](https://data.sccgov.org/) on November 3, 2022.

* [`CA_2020_ADI_Census_Block_Group_v3.2.zip`](https://www.neighborhoodatlas.medicine.wisc.edu/): Census tract-level [Area Deprivation Index (ADI)](https://www.nejm.org/doi/full/10.1056/NEJMp1802313) rankings for California, indexed by FIPS.
* [`CensusTract2020.csv`](https://data.sccgov.org/Government/CensusTract2020/4z77-invd): Santa Clara County census tract information from the 2020 census.
* [`COVID-19_Vaccination_among_County_Residents_by_Census_Tract.csv`](https://data.sccgov.org/COVID-19/COVID-19-Vaccination-among-County-Residents-by-Cen/qx2e-7jz2): Santa Clara County COVID-19 vaccination rates by census tract.
