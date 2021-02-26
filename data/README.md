# Metadata
Temperature data and DDRP results are large files that will be included in the final Zenodo archive, but are not uploaded currently in the repository.

**Directories**
* *daymet*: Not included, daily temperature min/max rasters for North America. Instructions on how to process these from downloaded netcdf's are in Process_gridded_temperature.R
* *eobs*: Not included, daily temperature min/max rasters for Europe. Instructions on how to process these from downloaded netcdf's are in Process_gridded_temperature.R
* *DDRP_results*: Not included, Intermediate results (.tif) from DDRP for 2010-2019 that are used to make maps in the manuscript in 05_map_DDRP_results.R. These results come from analysis performed in 04_DDRP_phenology model.R.

**Experimental data**
* *dev_exp1.csv*: Development time at constant temperatures
* replicate = unique ID for replicate and temperature treatment
* time1 = days to complete development
* temperature = constant temperature treatment in degrees Celsius.
* *dev_exp2.csv*: Development time after transfer from low temperature treatments
* plant = experimental unit, individuals developing on same plant
* run = treatment replicates
* time1 = 30 days at low temperature treatment
* time2 = days to complete development once transferred to 15C
* temperature = initial temperature treatment for 30 days in degrees Celsius.
* temp2 = all finished development at 15C
* replicate = unique ID for replicate by plant, run, and temperature treatment
* *dev_myint.csv*: Development times reported from publication
* sex = male/female (we grouped together to compare with our experiments)
* temperature = constant temperature treatments
* n = sample size, number of individuals
* mean and se for multiple lifestage development times and total development time
* *photo_diapause.csv*: Photoperiod-based diapause initiation experiment
* Pop = biotype collected in Japan, Kyushu (southern) or Hokkaido (northern)
* Treat = photoperiod treatment in hours of daylight
* Repro = number of females in treatment reproductive
* Diap = number of females in treatment that entered diapause
* fraction = proportion reproductive for plotting
* *springtime_oviposition.csv*: Photoperiod- and temperature-based diapause termination * experiment
* Population = biotype collected in Japan, Kyushu (southern) or Hokkaido (northern)
* Temperature = constant temperature treatment, 14 or 21C
* Photoperiod = photoperiod treatment in hours of daylight
* Date = date of observation
* AccumDD = accumulated degree-days from start of experiment (LDT 6.9C)
* Feeding = observed feeding from individuals in petri dish, lerp on leaf
* Mating = observed mating within petri dish,
* Dead = number of adults that had died
* Eggs_total = number of eggs in petri dish (top and bottom of leaf)
* Eggs_top = eggs visible on top of leaf
* Eggs_bottom = eggs hidden, counted less frequently
* EggAccumStartDate = date when oviposition started
* StartCount = number of adults in dish at beginning of experiment
* Notes = other observations

**Other**
* *na_lines.rds*: Boundaries for North America for maps, including both states and provinces.

