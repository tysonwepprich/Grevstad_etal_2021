# Analysis for Grevstad et al. manuscript

This repository walks through our analysis of *Aphalara itadori* developmental traits, photoperiod-based diapause, and predicted voltinism within North America and Europe for 2010-2019.

**In the code directory**, follow the numbered R scripts to reproduce the analysis. 
* 01_development_rates.R uses the results of three experiments to estimate the development rate at different temperatures and the total degree-days to complete a lifecycle.
* 02_photoperiod_response.R estimates the diapause proportion at different day hours from a growth chamber experiment.
* 03_spring_emergence.R estimates the degree-days needed for overwintering adults to break diapause and lay eggs from a growth chamber experiment.
* 04_DDRP_phenology_model.R uses the developmental traits to predict phenology and voltinism with the Degree-Days, Risk, and Phenological event mapping (DDRP) platform by [Barker et al. (2020)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0244005). This code is abbreviated from [their repository](https://github.com/bbarker505/ddrp_v2) with photoperiod features added for this manuscript. Only 2019 is predicted from the temperature data here.
* 05_map_DDRP_results.R uses 2010-2019 results from DDRP and makes maps for the manuscript.

Unnumbered R scripts are for utility functions and processing temperature rasters. 
* DDRP_cohorts_v1_funcs.R are functions for DDRP (not all used in this manuscript)
* Process_gridded_temperatures.R shows how to convert downloaded gridded temperature data into daily rasters for DDRP. You do not have to run this, as we include 2019 as sample temperature data.
* APH.params is a parameters file for our species using the DDRP format

**In the data directory**, we include raw data and intermediate products described with a metadata file. Daily temperature rasters are too large to include all data in this repository, but instructions on downloading and processing temperature data are included in Process_gridded_temperatures.R. We include 2019 data from Europe as an example.

**In the figures directory**, we include the output from the analysis that we include in the manuscript.