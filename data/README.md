# Raw Data

The raw data files for each location should be included in subdirectories with the following structure:
```
_general/
    Vaccines.gov__COVID-19_vaccinating_provider_locations.csv
ca/
    cageo2020.pl
chicago/
    COVID-19_Cases__Test_and_Deaths_by_ZIP_Code.csv
    COVID-19_Vaccination_Locations.csv
    COVID-19_Vaccinations_by_ZIP_Code.csv
    IL_2020_ADI_9 Digit Zip Code_v3.2.csv
santa_clara/
    CA_2020_ADI_Census_Block_Group_v3.2.csv
    COVID-19_Vaccination_among_County_Residents_by_Census_Tract.csv
    Santa_Clara_County_Pharmacies.csv
```
The files for each location, along with their sources, are explained below.

## General

* [`Vaccines.gov__COVID-19_vaccinating_provider_locations.csv`](https://data.cdc.gov/Vaccinations/Vaccines-gov-COVID-19-vaccinating-provider-locatio/5jp2-pgaw): The vaccines.gov master list of vaccine providers gathered March 21, 2023.

## CA (California)

* [`cageo2020.pl`](https://www2.census.gov/programs-surveys/decennial/2020/data/01-Redistricting_File--PL_94-171/): Geographic header file for the California portion of the 2020 Census National Redistricting Data Summary File, last updated August 12, 2021. See the header file's documentation [here](https://www2.census.gov/programs-surveys/decennial/2020/technical-documentation/complete-tech-docs/summary-file/2020Census_PL94_171Redistricting_NationalTechDoc.pdf).

## Chicago

Data mostly collected from the [Chicago Data Portal](https://data.cityofchicago.org/) on November 7, 2022.

* [`COVID-19_Cases__Test_and_Deaths_by_ZIP_Code.csv`](https://data.cityofchicago.org/Health-Human-Services/COVID-19-Cases-Tests-and-Deaths-by-ZIP-Code/yhhz-zm2v): Historical data for Chicago COVID-19 cases and deaths by ZIP code. Also includes centroids of ZIP codes and populations by ZIP code.
* [`COVID-19_Vaccination_Locations.csv`](https://data.cityofchicago.org/Health-Human-Services/COVID-19-Vaccination-Locations/6q3z-9maq): Geographic data for COVID-19 vaccine locations in Chicago. A couple of geographic coordinates missing from the data file have been filled in by hand using the given street address.
* [`COVID-19_Vaccinations_by_ZIP_Code.csv`](https://data.cityofchicago.org/Health-Human-Services/COVID-19-Vaccinations-by-ZIP-Code/553k-3xzc): Historical data for Chicago COVID-19 vaccination rates by ZIP code. Also includes centroids of ZIP codes.
* [`IL_2020_ADI_9 Digit Zip Code_v3.2.csv`](https://www.neighborhoodatlas.medicine.wisc.edu/): ZIP code-level [Area Deprivation Index (ADI)](https://www.nejm.org/doi/full/10.1056/NEJMp1802313) rankings for Illinois, indexed by ZIP code.

## Santa Clara County

Data mostly collected from the [County of Santa Clara Open Data Portal](https://data.sccgov.org/) on November 3, 2022.

* [`CA_2020_ADI_Census_Block_Group_v3.2.csv`](https://www.neighborhoodatlas.medicine.wisc.edu/): Census tract-level [Area Deprivation Index (ADI)](https://www.nejm.org/doi/full/10.1056/NEJMp1802313) rankings for California, indexed by FIPS.
* [`COVID-19_Vaccination_among_County_Residents_by_Census_Tract.csv`](https://data.sccgov.org/COVID-19/COVID-19-Vaccination-among-County-Residents-by-Cen/qx2e-7jz2): Santa Clara County COVID-19 vaccination rates by census tract.
* [`Santa_Clara_County_Pharmacies.csv`](https://data.cdc.gov/Vaccinations/Vaccines-gov-COVID-19-vaccinating-provider-locatio/5jp2-pgaw): A data table of Santa Clara County pharmacy information extracted from a the vaccines.gov master list during March-April 2023. Primarily used to generate the facility location file by geocoding street addresses to find geographic coordinates.
