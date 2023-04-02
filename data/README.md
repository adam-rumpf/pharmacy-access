# Raw Data

The raw data files for each location should be included in subdirectories with the following structure:
```
chicago/
    COVID-19_Cases__Test_and_Deaths_by_ZIP_Code.csv
    COVID-19_Vaccination_Locations.csv
    COVID-19_Vaccinations_by_ZIP_Code.csv
    IL_2020_ADI_9 Digit Zip Code_v3.2.csv
santa_clara/
    2022_gaz_tracts_06.txt
    CA_2020_ADI_Census_Block_Group_v3.2.csv
    COVID-19_Vaccination_among_County_Residents_by_Census_Tract.csv
    Santa_Clara_County_Pharmacies.csv
```
The files for each location, along with their sources, are explained below.

## Chicago

Data mostly collected from the [Chicago Data Portal](https://data.cityofchicago.org/) on November 7, 2022.

* [`COVID-19_Cases__Test_and_Deaths_by_ZIP_Code.csv`](https://data.cityofchicago.org/Health-Human-Services/COVID-19-Cases-Tests-and-Deaths-by-ZIP-Code/yhhz-zm2v): Historical data for Chicago COVID-19 cases and deaths by ZIP code. Also includes centroids of ZIP codes and populations by ZIP code.
* [`COVID-19_Vaccination_Locations.csv`](https://data.cityofchicago.org/Health-Human-Services/COVID-19-Vaccination-Locations/6q3z-9maq): Geographic data for COVID-19 vaccine locations in Chicago. A couple of geographic coordinates missing from the data file have been filled in by hand using the given street address.
* [`COVID-19_Vaccinations_by_ZIP_Code.csv`](https://data.cityofchicago.org/Health-Human-Services/COVID-19-Vaccinations-by-ZIP-Code/553k-3xzc): Historical data for Chicago COVID-19 vaccination rates by ZIP code. Also includes centroids of ZIP codes.
* [`IL_2020_ADI_9 Digit Zip Code_v3.2.csv`](https://www.neighborhoodatlas.medicine.wisc.edu/): ZIP code-level [Area Deprivation Index (ADI)](https://www.nejm.org/doi/full/10.1056/NEJMp1802313) rankings for Illinois, indexed by ZIP code.

## Santa Clara County

Data mostly collected from the [County of Santa Clara Open Data Portal](https://data.sccgov.org/) on November 3, 2022.

* [`2022_gaz_tracts_06.txt`](https://www.census.gov/geographies/reference-files/time-series/geo/gazetteer-files.html): California gazetteer files from the 2020 census.
* [`CA_2020_ADI_Census_Block_Group_v3.2.csv`](https://www.neighborhoodatlas.medicine.wisc.edu/): Census tract-level [Area Deprivation Index (ADI)](https://www.nejm.org/doi/full/10.1056/NEJMp1802313) rankings for California, indexed by FIPS.
* [`COVID-19_Vaccination_among_County_Residents_by_Census_Tract.csv`](https://data.sccgov.org/COVID-19/COVID-19-Vaccination-among-County-Residents-by-Cen/qx2e-7jz2): Santa Clara County COVID-19 vaccination rates by census tract.
* [`Santa_Clara_County_Pharmacies.csv`](https://www.vaccines.gov/): A data table of Santa Clara County pharmacy information found by manual search of vaccines.gov during Spring 2023. Primarily used to generate the facility location file by geocoding street addresses to find geographic coordinates.
