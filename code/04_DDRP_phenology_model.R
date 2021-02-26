## Header ----
## Script name: 04_DDRP_phenology_model.R
## Purpose of script:
## Authors: Brittany S. Barker, Len Coop, Gericke Cook, Dan Upper, Tyson Wepprich
## Date Created: 2021-02-04
## License: GNU GPLv3
## Email: tyson.wepprich@gmail.com
## ---
## Notes:
## Modified by Tyson Wepprich from 
## Barker et al. (2020). DDRP: Real-time phenology and climatic suitability 
## modeling of invasive insects. PLOS ONE. https://doi.org/10.1371/journal.pone.0244005
## Archived code from Barker et al. (2020) https://doi.org/10.5281/zenodo.4294801
## This script adds photoperiod-based diapause/phenology mismatches and removes
## other features not included in this manuscript, such as phenology event maps 
## and lifestage proportions across the season.

## For analysis, this script was run from terminal on linux server. 
## It is reproduced here without the need for shell (40 min. on laptop for North America, 6 min. for Europe).
## This script demonstrates the DDRP model for 2019. 
## The annual maps produced in this script are not in the manuscript.
## All results 2010-2019 were produced with this script, and are saved in /DDRP_results
## Annual results (*_all_weighted.tif) are then used in 05_map_DDRP.results.R to make maps for manuscript
## ---

# Load the required packages
pkgs <- c("colorspace", "doParallel", "plyr", "dplyr", "foreach", "ggplot2", "ggthemes", 
          "lubridate", "mapdata", "mgsub", "optparse", "parallel",
          "purrr", "RColorBrewer", "rgdal", "raster", "readr", "sp", "stringr", 
          "tidyr", "tictoc", "tools", "viridis")

lapply(pkgs, library, character.only = TRUE)

# Load collection of functions for this model
source("code/DDRP_cohorts_v1_funcs.R")

# Start timing the model run
tic("Total run time")

# (1). PARAM HANDLING -----

#### * Values for params, if not provided in command line ####
spp           <- "APH" # Default species to use
forecast_data <- "EOBS" # Forecast data to use (DAYMET or EOBS)
start_year    <- 2019 # Year to use
start_doy     <- 1 # Start day of year          
end_doy       <- 365 # End day of year - need 365 if voltinism map 
keep_leap     <- 0 # Should leap year be kept?
region_param  <- "EUROPE" # Default REGION to use
exclusions_stressunits    <- 0 # Turn on/off climate stress unit exclusions
pems          <- 0 # Turn on/off pest event maps
out_dir       <- "example_EURO_HOKKAIDO" # Output dir
out_option    <- 1 # Output option category
ncohort       <- 6 # Number of cohorts to approximate end of OW stage
odd_gen_map   <- 0 # Create summary plots for odd gens only (gen1, gen3, ..)
do_photo      <- 1 # Use photoperiod diapause modules in daily loop and results
cp_mean       <- 14.7 # Critical photoperiod mean
cp_sd         <- 0.284 # Standard deviation around cp_mean


# (2). DIRECTORY INIT ------

#### * Param inputs - species params; thresholds, weather, etc. ####
params_dir <- "C:/Users/tmwepprich/Desktop/repo/Grevstad_etal_2021/code/"

# Bring in states feature for summary maps (PNG files)
# Requires these libraries: "mapdata" and "maptools"
cat("\n\nDownloading US states feature\n")
if (region_param == "EUROPE"){
  states <- map_data("world")
}else{
  states <- map_data("state")
}

#### * Weather inputs and outputs - climate data w/subdirs 4-digit year ####
if(region_param %in% c("CONUSPLUS", "WESTPLUS")) forecast_data <- "DAYMET"
if(end_doy == 366) end_doy <- 365 # ignore leap years

# if(forecast_data == "PRISM"){
#   base_dir <- "/data/PRISM/"
#   prism_dir <- paste0(base_dir, start_year)
#   raster_crs <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0 "
# }

if(forecast_data == "DAYMET"){
  # North America Daymet data
  base_dir <- "C:/Users/tmwepprich/Desktop/repo/Grevstad_etal_2021/data/daymet/"
  prism_dir <- paste0(base_dir, start_year)
  raster_crs <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0" 
}

# if(forecast_data == "MACA"){
#   base_dir <- "/home/macav2metdata/GFDL-ESM2M/"
#   prism_dir <- paste0(base_dir, start_year)
#   raster_crs <- "+proj=longlat +pm=0 +a=6378137 +rf=298.257223563" 
# }

if(forecast_data == "EOBS"){
  base_dir <- "C:/Users/tmwepprich/Desktop/repo/Grevstad_etal_2021/data/eobs/"
  prism_dir <- paste0(base_dir, start_year)
  raster_crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0" 
}


cat("\nBASE DIR: ", base_dir, "\n")
cat("\nWORKING DIR: ", prism_dir, "\n")

#### * Output directory, log file, and error message file ####
# MUST remove .tif files or script will crash during processing because it will 
# try to analyze previously processed results. 

output_dir <- paste0("C:/Users/tmwepprich/Desktop/repo/Grevstad_etal_2021/data/DDRP_results/", out_dir)

# Remove all files if output_dir exists, or else create output_dir
if (file.exists(output_dir)) {
  unlink(paste0(output_dir, "/*"), recursive = TRUE, force = TRUE)
  cat("\n", str_wrap(paste0("EXISTING OUTPUT DIR: ", output_dir, 
                            "; removing all files\n"), width = 80), sep = "") 
} else {
  dir.create(output_dir)
  cat("NEW OUTPUT DIR:", output_dir, "\n")
}

# Push out a rlogging file with all main messages in model
# Put all log, message, and metadata files in a separate folder
setwd(output_dir)
dir.create("Logs_metadata")
Model_rlogging <- sprintf("%s%s", "./", "/Logs_metadata/Model_rlogging.txt")

# Make header for logging file
cat(paste0(rep("#", 36), collapse = ""), "\n", 
    "### Log file for DDRP cohorts v1 ###\n", 
    paste0(rep("#", 36), collapse = ""), "\n\n", sep = "", 
    file = Model_rlogging)

# Record PRISM and output dir
cat("BASE DIR: ", base_dir, "\n", file = Model_rlogging, append = TRUE)
cat("WORKING DIR: ", prism_dir, "\n", file = Model_rlogging, append = TRUE)
cat(str_wrap(paste0("EXISTING OUTPUT DIR: ", output_dir, 
                    "; removing all files"), width = 80), "\n\n", sep = "", 
    file = Model_rlogging, append = TRUE)

# Push out a message file with all R error messages
msg <- file("Logs_metadata/rmessages.txt", open = "wt")
sink(msg, type = "message")

# (3). PARAMETER AND SETTINGS SETUP ----- 
cat("PARAMETER AND SETTINGS SETUP: getting parameters to use in model\n", 
    file = Model_rlogging, append = TRUE)
cat("\n\nPARAMETER AND SETTINGS SETUP: getting parameters to use in model\n")

# Read from source param files in ./spp_params/SPP.params
param_file <- sprintf("%s%s", spp, ".params")
spp <- gsub(".params", "", param_file) # Get species abbr.
species_params <- sprintf("%s%s", params_dir, param_file) # Location of file

if (file.exists(species_params)) {
  cat("Species params: ",species_params, "\n", file = Model_rlogging, 
      append = TRUE)
  source(species_params) # Read in species parameters
  cat("Reading params for species: ", spp, " Fullname: ", fullname, "\n", 
      file = Model_rlogging, append = TRUE)
  cat("\nReading params for species: ", spp, " Fullname: ", fullname, "\n")
} else {
  cat("Param file: ", species_params, "...not found; exiting program\n", 
      file = Model_rlogging, append = TRUE)
  cat("\nParam file: ", species_params, "...not found; exiting program\n")
  q()  # No reason to keep going without any params
}

# Change year to numeric if it's a specific year
# If using climate normals, there may be letters in folder name
if (!grepl("[A-z]", start_year)) {
  start_year <- as.numeric(start_year)
}

# Create a list of days to use for daily loop
sublist <- start_doy:end_doy

#### * Format threshold and DD params for Daily Loop ####

# Need to match length and order of stgorder, which is species specific
# Remove "F" from stgorder param - this is for an old version of DDRP
# stgorder <- stgorder[-length(stgorder)]

# Upper and lower thresholds
# OW stage will have the same threshold as actual stage (e.g., OWadult = adult)
# Need to match LDT or UDT value to stage, which first requires changing 
# stage abbr to the stage name ("E" = "egg")
stage_ldt_list <- list()
stage_udt_list <- list()
j <- 1

for (i in 1:length(stgorder)) {
  stg_nam <- stgorder[i]
  stg_nam <- mgsub(string = stg_nam, # Requires "mgsub" package
                   pattern = c("OE", "OL", "OP", "OA", "E", "L", "P", "A"), 
                   replacement = c("egg", "larvae", "pupae", "adult", "egg", 
                                   "larvae", "pupae", "adult"))
  stage_ldt_val <- get(paste0(stg_nam, "LDT")) # returns LDT value for stage
  stage_ldt_list[[j]] <- stage_ldt_val
  stage_udt_val <- get(paste0(stg_nam, "UDT")) # returns UDT value for stage
  stage_udt_list[[j]] <- stage_udt_val
  j <- j + 1
}

# DD parameters - OW stage has it's own DD, so change stage abbr. to stage 
# name or "OW" plus stage name
stage_dd_list <- list()
j <- 1

for (i in 1:length(stgorder)) {
  stg_nam <- stgorder[i]
  stg_nam <- mgsub(string = stg_nam, 
                   pattern = c("OE", "OL", "OP", "OA", "E", "L", "P", "A"), 
                   replacement = c("OWegg", "OWlarvae", "OWpupae", "OWadult", 
                                   "egg", "larvae", "pup", "adult"))
  stage_dd_val <- get(paste0(stg_nam, "DD")) # returns DD value for stage
  stage_dd_list[[j]] <- stage_dd_val
  j <- j + 1
}

# Put the values in the list into a numeric vector
stage_ldt <- as.numeric(do.call(cbind, stage_ldt_list))
stage_udt <- as.numeric(do.call(cbind, stage_udt_list))
stage_dd <- as.numeric(do.call(cbind, stage_dd_list))

# Cohort response parameters
# SD is equal the square root of the variance (sigma-squared) parameter
# May consider changing length.out to something else in future versions...
# Also consider making the percent (in cohort_distro) an input parameter
xdist <- seq(xdist1, xdist2, length.out = 1000)
ydist <- dnorm(xdist, mean = distro_mean, sd = sqrt(distro_var))

inputdist <- data.frame(x = xdist, y = ydist) %>% 
  arrange(x) %>% 
  mutate(CDF = cumsum(y/sum(y)))
cohort_distro <- CohortDistrib(dist = inputdist, numstage = ncohort, perc = .99)
relpopsize <- cohort_distro$weights

# Parameters of required degree-days
# Replace OW gen with emergence distro values
if (ncohort == 1) {
  ddpar <- matrix(stage_dd, nrow = 1, byrow = TRUE)
  stage_dd <- ddpar
} else {
  ddpar <- cbind(cohort_distro$means, 
                 matrix(rep(stage_dd[-1], nrow(cohort_distro)), 
                        nrow = nrow(cohort_distro), byrow = TRUE))
  stage_dd <- ddpar
}

# (4). METADATA OUTPUT FILE -----

# Push out a metadata file with all inputs used in model
cat("\nMETADATA: creating metadata file for all inputs used in model\n", 
    file = Model_rlogging, append = TRUE)
cat("\nMETADATA: creating metadata file for all inputs used in model\n")

# Create the metadata file
setwd(output_dir)
metadata <- sprintf("%s%s", "./", "/Logs_metadata/metadata.txt")
cat("### Metadata for DDRP cohorts v1 ###\n", file = metadata)

# Document species information
cat("\n### Model Species Parameters ###\n Species Abbrev:", spp, 
    "\n Full Name:", fullname, 
    "\n Pest of:", pestof,
    "\n Overwintering Stage:", owstage, file = metadata, append = TRUE)

# Document developmental threshold temperatures
cat("\n \n Developmental threshold temperatures",
    "\n Egg Lower Devel Threshold:", eggLDT, 
    "\n Egg Upper Devel Threshold:", eggUDT, 
    "\n Larvae Lower Devel Threshold:", larvaeLDT, 
    "\n Larvae Upper Devel Threshold:", larvaeUDT, 
    "\n Pupae Lower Devel Threshold:", pupaeLDT,
    "\n Pupae Upper Devel Threshold:", pupaeUDT, 
    "\n Adult Lower Devel Threshold:", adultLDT, 
    "\n Adult Upper Devel Threshold:", adultUDT, file =  metadata, append = T)

# Document stage durations
cat("\n\n Stage durations in degree-days (DDs)",
    "\n Egg DDs:", eggDD, 
    "\n Larvae DDs", larvaeDD, 
    "\n Pupae DDs:", pupDD, 
    "\n Adult DDs:", adultDD, "\n ", 
    file = metadata, append = TRUE)

# Document climate stress exclusion parameter values, if applicable
if (exclusions_stressunits) {
  cat("\n Climate stress parameters",
      "\n Lower Chill Threshold:", chillstress_threshold, 
      "\n Upper Heat Threshold:", heatstress_threshold,
      "\n Max Chill Units (lower bound):", chillstress_units_max1, 
      "\n Max Chill Units (upper bound):", chillstress_units_max2,
      "\n Max Heat Stress Units (lower bound):", heatstress_units_max1,
      "\n Max Heat Stress Units (upper bound):", heatstress_units_max2, 
      file = metadata, append = TRUE)
}

# Document Pest Event Map parameter values, if applicable
if (pems) {
  cat("\n \n Pest Event Map parameters",
      "\n Number of generations to make Pest Event Maps (PEMs): ", PEMnumgens,
      "\n Egg Event DDs and Label: ", eggEventDD, " (", eggEventLabel,")", 
      "\n Larvae Event DDs and Label: ", larvaeEventDD, " (", larvaeEventLabel, ")",
      "\n Pupae Event DDs and Label: ", pupaeEventDD, " (", pupaeEventLabel, ")",
      "\n Adult Event DDs and Label: ", adultEventDD, " (", adultEventLabel, ")",
      sep = "", file = metadata, append = TRUE)
}

cat("\n\n### Model Input Parameters ###\n Start Year:", start_year, 
    "\n Weather data for forecasts:", forecast_data, 
    "\n Start day-of-year:", start_doy,
    "\n End day-of-year:", end_doy, 
    "\n Region:", region_param, 
    "\n Climate stress exclusion maps:", exclusions_stressunits, 
    "\n Pest Event Maps:", pems,    
    "\n Output_Dir:", out_dir, 
    "\n Output option:", out_option, 
    "\n No. of cohorts:", ncohort, 
    "\n Mean of end of OW stage (DDs):", distro_mean, 
    "\n Low bound of end of OW stage (DDs), xdist:", xdist1, 
    "\n High bound of end of OW stage (DDs), ydist:", xdist2, 
    "\n Variance in end of OW stage (DDs):", distro_var,
    "\n Shape of distribution of end of OW stage (DDs):", distro_shape, 
    "\n Plot odd gens only:", odd_gen_map, 
    file = metadata, append = TRUE)

# Make a table of stage DDs for each cohort and print to metadata
stage_dd.print <- as.data.frame(stage_dd)
stage_dd.print[1] <- round(stage_dd.print[1], 0)
colnames(stage_dd.print) <- stgorder
stage_dd.print <- cbind("cohort" = as.integer(rownames(stage_dd.print)), 
                        data.frame(stage_dd.print, row.names = NULL))

cat("\n\n###Durations (in degree-days) of stages in each of", 
    ncohort, "cohorts ###\n", file = metadata, append = TRUE) 
suppressWarnings(write.table(stage_dd.print, file = metadata, 
                             row.names = FALSE, 
                             col.names = colnames(stage_dd.print), 
                             append = TRUE))

cat("\nDurations (degree-days) of stages in each of", ncohort, "cohorts:\n\n") 
print(stage_dd.print, row.names = FALSE)

cat("\nDone writing metafile\n\n", forecast_data, " DATA PROCESSING\n", sep = "",
    file = Model_rlogging, append = TRUE)
cat("\nDone writing metafile\n\n", forecast_data, " DATA PROCESSING\n",
    sep = "")

# (5). WEATHER DATA LOADING AND PROCESSING -----

# Weather inputs and outputs - PRISM climate data w/subdirs 4-digit year
# New feature - choose whether to use PRISM or NMME for weather forecasts 
# (forecast_data = PRISM, or forecast_data = NMME)
if(forecast_data == "PRISM"){
  tminfiles <- list.files(path = prism_dir,
                          pattern = glob2rx(paste0("*PRISM_tmin_*",
                                                   start_year, "*.bil$*")),
                          all.files = FALSE, full.names = TRUE, recursive = TRUE)
  tminfiles <- ExtractBestPRISM(tminfiles, forecast_data, 
                                keep_leap)[start_doy:end_doy]
  
  tmaxfiles <- list.files(path = prism_dir,
                          pattern = glob2rx(paste0("*PRISM_tmax_*",
                                                   start_year, "*.bil$*")),
                          all.files = FALSE, full.names = TRUE, recursive = TRUE)
  tmaxfiles <- ExtractBestPRISM(tmaxfiles, forecast_data, 
                                keep_leap) [start_doy:end_doy]
}



if(forecast_data == "MACA"){
  tminfiles <- list.files(path = prism_dir,
                          pattern = glob2rx(paste0("*MACAV2_tmin_*",
                                                   start_year, "*.grd$*")),
                          all.files = FALSE, full.names = TRUE, recursive = TRUE)[start_doy:end_doy]
  
  tmaxfiles <- list.files(path = prism_dir,
                          pattern = glob2rx(paste0("*MACAV2_tmax_*",
                                                   start_year, "*.grd$*")),
                          all.files = FALSE, full.names = TRUE, recursive = TRUE)[start_doy:end_doy]
}

if(forecast_data == "EOBS"){
  tminfiles <- list.files(path = prism_dir,
                          pattern = glob2rx(paste0("*EOBS_tmin_*",
                                                   start_year, "*.grd$*")),
                          all.files = FALSE, full.names = TRUE, recursive = TRUE)[start_doy:end_doy]
  
  tmaxfiles <- list.files(path = prism_dir,
                          pattern = glob2rx(paste0("*EOBS_tmax_*",
                                                   start_year, "*.grd$*")),
                          all.files = FALSE, full.names = TRUE, recursive = TRUE)[start_doy:end_doy]
}

if(forecast_data == "DAYMET"){
  
  tminfiles <- list.files(path = prism_dir,
                          pattern = glob2rx(paste0("*DAYMET_tmin_*",
                                                   start_year, "*.grd$*")),
                          all.files = FALSE, full.names = TRUE, recursive = TRUE)[start_doy:end_doy]
  
  tmaxfiles <- list.files(path = prism_dir,
                          pattern = glob2rx(paste0("*DAYMET_tmax_*",
                                                   start_year, "*.grd$*")),
                          all.files = FALSE, full.names = TRUE, recursive = TRUE)[start_doy:end_doy]
  
}
## Extract date from temperature files using regex pattern matching
dats <- regmatches(tminfiles, regexpr(pattern = "[0-9]{8}", text = tminfiles))

# Specify sampling frequency (how many days until output maps are generated?)
# This feature may be removed in production version
if (out_option == 0) {
  sample_freq <- 61 # Bimonthly maps
} else if (out_option == 1) {
  sample_freq <- 30 # Monthly maps
} else if (out_option == 2) {
  sample_freq <- 14  # Biweekly maps
} else if (out_option == 3) {
  sample_freq <- 10 # Dekad maps
} else if (out_option == 4) {
  sample_freq <- 7  # Weekly maps
} else if (out_option == 5) {
  sample_freq <- 2  # Even day maps
} else if (out_option == 6) {
  sample_freq <- 1 # Daily maps
  dats2 <- dats
} else if (out_option %in% !c(0, 1, 2, 3, 4, 5, 6)) {
  cat("Error: out_option =", out_option, "is unacceptable; exiting program\n", 
      file = Model_rlogging, append = TRUE)
  cat("Error: out_option =", out_option, "is unacceptable; exiting program\n")
  q()  # No reason to keep going if sampling freq is not correctly specified
}

# Make vector of dates to use when processing results - last day is added too
# Using "unique" will only keep date if it doesn't already occur in vector
# This happens if the end day of year is a multiple of the sampling frequency 
# (e.g. 1 to 300, w/ a 30 day sampling frequency)
dats2 <- unique(c(dats[seq(0, length(dats), sample_freq)], last(dats)))

# Create vector of days in the sublist that will be sampled (rasters are saved 
# for those days). Then tack on last day for sampling; using "unique" will only 
# keep day if it doesn't already occur in vector
sample_pts <- sublist[seq(0, length(sublist), sample_freq)]
sample_pts <- unique(c(sample_pts, last(sublist)))

# Log file and terminal messages
cat("Finished loading ", forecast_data, " files for ", 
    length(start_doy:end_doy), " days\nCreating template file for ", 
    region_param, "\n", sep = "", file = Model_rlogging, append = TRUE)
cat("\nFinished loading ", forecast_data, " files for ", 
    length(start_doy:end_doy), " days\n\nCreating template file for ", 
    region_param, "\n", sep = "")

### * Create blank template from a temp file
# This template is used for cropping the temperature (tmin, tmax) rasters
REGION <- Assign_extent(region_param) # Bounding box
template <- crop(raster(tminfiles[1]), REGION) # Template for cropping
template[!is.na(template)] <- 0
dataType(template) <- "INT2U"
# crs <- crs("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0, 0, 0")
crs1 <- crs(raster_crs)
crs(template) <- crs1

#### * If CONUS or EAST, split template into tiles (and run in parallel)
# Benefit of tiles is lost for smaller regions, so these are not split
# SpaDES.tools requires the 'sf' package, which requires a newer version of GDAL
# .inorder must be set to TRUE so that output files are in correct order!

# Register DoParallel
# The "RegCluster" function determines an appropriate # of cores depending on 
# the "region_param" and "ncohort" parameters, so the server doesn't become
# overloaded
RegCluster(ncohort)

if (region_param %in% c("CONUS", "EAST", "CONUSPLUS", "LOCO")) {
  # Split template (2 pieces per side)
  tile_list <- SplitRas(template, ppside = 2, save = FALSE, plot = FALSE) 
  tile_n <- 1:length(tile_list) # How many tiles?
  cat("Splitting template into", length(tile_list), "tiles\n", 
      file = Model_rlogging, append = TRUE)
  
  # Name the 4 tiles (tile1, tile2, tile3, tile4)
  template <- mapply(function(n, t) {
    names(n) <- paste0("tile", t)
    # values(n) <- 0 # Is this necessary? Obscures boundaries by replacing NA's
    return(n)
  }, n = tile_list, t = tile_n )
  
  # Crop temp files by each template tile
  cat("Cropping tmax and tmin tiles for", region_param, "\n", 
      file = Model_rlogging, append = TRUE)
  cat("\nCropping tmax and tmin tiles for", region_param, "\n")
  
  tmax_list <- foreach(tile = template, .packages = "raster", 
                       .inorder = FALSE) %:% 
    foreach(tmax = tmaxfiles, .packages = "raster", .inorder = TRUE) %dopar% { 
      m <- as.matrix(crop(raster(tmax), tile))
    }
  
  tmin_list <- foreach(tile = template, .packages = "raster", 
                       .inorder = FALSE) %:% 
    foreach(tmin = tminfiles, .packages = "raster", .inorder = TRUE) %dopar% { 
      m <- as.matrix(crop(raster(tmin), tile))
    }
  
  # If region is not CONUS or EAST, simply crop temp files by the single template
} else {
  cat("Cropping tmax and tmin tiles for", region_param, "\n", 
      file = Model_rlogging, append = TRUE)
  cat("\nCropping tmax and tmin tiles for", region_param, "\n")
  
  tmax_list <- foreach(t = tmaxfiles, .packages = "raster", 
                       .inorder = TRUE) %dopar% {
                         m <- as.matrix(crop(raster(t), template))
                       }
  
  tmin_list <- foreach(t = tminfiles, .packages = "raster", 
                       .inorder = TRUE) %dopar% {
                         m <- as.matrix(crop(raster(t), template))
                       }
  
  stopCluster(cl)
}

cat("Done processing ", forecast_data, " data\n\nDAILY LOOP\n", sep = "",
    file = Model_rlogging, append = TRUE)
cat("\nDone processing ", forecast_data, " data\n\nDAILY LOOP\n", sep = "")

## (6). RUN THE DAILY LOOP -----
# First set up number of cores to use
# IMPORTANT: using too many cores will result in low memory, killing the daily 
# loop part-way through. For "mclapply" functions, set mc.cores manually.

# For mclapply and mcmapply: Can not use >1 core on Windows, so detect OS
if (grepl("Windows", Sys.info()[1])) {
  mc.cores <- 1
} else {
  mc.cores <- 4 # use 4 here, because there are 4 tiles being run in parallel
}

# Split cohorts into smaller chunks for CONUS and EAST to avoid overloading 
# memory when running in parallel. Three cohorts puts load up to ~13.
# if (region_param %in% c("CONUS", "EAST", "CONUSPLUS", "LOCO")) {
cohort_chunks <- split(1:ncohort, ceiling(1:length(1:ncohort)/3)) 
# } else {
#   cohort_chunks <- split(1:ncohort, ceiling(1:length(1:ncohort)/5))
# }

tic("Daily loop run time") # Start timing the daily loop run-time
# cat("DAILY LOOP: daily loop log files show loop progress and 
# output file info\n", file = Model_rlogging, append = TRUE)
cat("Sampling every", sample_freq, "days between", first(dats), "and", 
    last(dats), "\n", file = Model_rlogging, append = TRUE) 
cat("\nSampling every", sample_freq, "days between", first(dats), "and", 
    last(dats), "\n") 

RegCluster(ncohort)

# Run it! If there is an error the program will stop
tryCatch(
  {
    # If the region is CONUS or EAST, then both cohorts and tiles will be run in
    # parallel. To avoid overloading the server, mc.cores = 4 (for the 4 tiles,
    # which keeps load < 12. The run-time is approx 1.5-2.5 minutes.
    if (region_param %in% c("CONUS", "EAST", "CONUSPLUS", "LOCO")) {
      # Total number of nodes is mc.cores * 2 because use mclapply twice in loop
      for (c in cohort_chunks) {
        cat("Running daily loop for cohorts", as.character(c), "\n", 
            file = Model_rlogging, append = TRUE)
        cat("\nRunning daily loop for cohorts", as.character(c), "\n")
        cohort_vec <- unname(unlist(c)) # change to an unnamed vector
        # For some reason foreach here just doesn't work well for CONUS 
        # (slow, way too much load) - not sure why
        mclapply(cohort_vec, function(cohort) {
          mclapply(1:length(template), function(tile_num) {
            tile <- template[[tile_num]]
            DailyLoop(cohort, tile_num, tile) 
          }, mc.cores = mc.cores)
        }, mc.cores = mc.cores)  
      }
    } else {
      # If the region is not CONUS or EAST, then we don't need to run function 
      # for multiple tiles.
      for (c in cohort_chunks) {
        cat("Running daily loop for cohorts", as.character(c), "\n", 
            file = Model_rlogging, append = TRUE)
        cat("\nRunning daily loop for cohorts", as.character(c), "\n")
        cohort_vec <- unname(unlist(c)) # change to an unnamed vector
        foreach(cohort = cohort_vec, .packages = pkgs, 
                .inorder = FALSE) %dopar% {
                  DailyLoop(cohort, NA, template)
                }
        
      }
    }
  },
  error = function(e) {
    cat("Error in Daily Loop - stopped run - check rmessages file\n", 
        file = Model_rlogging, append = TRUE) 
    cat("\nError in Daily Loop - stopped run - check rmessages file\n")
  })

stopCluster(cl)

# Document daily loop execution time
loop_exectime <- toc(quiet = TRUE)
loop_exectime <- (loop_exectime$toc - loop_exectime$tic) / 60 

cat("Daily loop done (run time = ", round(loop_exectime, digits = 2), " min)",
    "\n\nFINAL ANALYSES AND MAP PRODUCTION\n", sep = "",
    file = Model_rlogging, append = TRUE) 
cat("\nDaily loop done (run time = ", round(loop_exectime, digits = 2), " min)",
    "\n\nFINAL ANALYSES AND MAP PRODUCTION\n", sep = "")

## (7). PROCESS DAILY LOOP RESULTS -----
setwd(output_dir)
tic("Data processing run time") # Start timing for data processing

# Create a directory ("Misc_output") to put secondary outfiles
dir.create("Misc_output")

#### * If CONUS or EAST, merge and delete tiles ####
# If CONUS or EAST, merge the tiles
if (region_param %in% c("CONUS", "EAST", "CONUSPLUS", "LOCO")) {
  cat("\nMerging tiles for", region_param, "\n\n", 
      file = Model_rlogging, append = TRUE)
  cat("\nMerging tiles for", region_param, "\n")
  # Get list of brick files for each tile, the type of file, 
  # its cohort number, and then make a list of the file types
  # File types are split up to avoid overloading sytem when running in parallel
  brick_files <- list.files(pattern = glob2rx("*cohort*tile*tif$"), 
                            recursive = FALSE)
  type_list <- unique(str_split_fixed(brick_files, 
                                      pattern = "_cohort", 4)[,1]) 
  type_list_split <- split(type_list, ceiling(1:length(type_list)/3)) 
  cohorts <- unique(str_split_fixed(brick_files, pattern = "cohort", 4)[,2])
  cohorts <- unique(substring(cohorts, 1, 1)) # Vector of cohorts (1, 2, 3, ...)
  
  # For each file type, merge the tiles for all cohorts
  # File type and cohorts are both run in parallel to increase speed
  # If system is overloaded then consider: 
  # 1) splitting up the cohorts list (as in the "type_list_split"); or 
  # 2) decreasing the number of splits in "type_list_split")
  RegCluster(ncohort)
  
  mrg_by_type <- foreach(type = type_list_split, .packages = pkgs, 
                         .inorder = FALSE) %dopar% {
                           type_vec <- unname(unlist(type)) # Change to an unnamed vector
                           for (t in type_vec) {
                             # If type is exclusion, stressunits, or ddtotal files, 
                             # then just do cohort 1; no other cohorts present
                             if (grepl("Stress_Excl|Stress_Units|DDtotal", t)) {
                               CombineMaps(brick_files, t, "1")
                               cat("Merged", t, "tiles for cohort 1\n", 
                                   file = Model_rlogging, append = TRUE)
                               # If another file type, then merge tiles for all cohorts
                             } else {
                               MrgTiles <- foreach(c = cohorts, .packages = pkgs, 
                                                   .inorder = TRUE) %dopar% {
                                                     CombineMaps(brick_files, t, c)
                                                     cat("Merged", t, "tiles for cohort", c, "\n", 
                                                         file = Model_rlogging, append = TRUE)
                                                   }
                             } 
                           }
                         }
  stopCluster(cl)
  cat("\nDone merging tiles\n", file = Model_rlogging, append = TRUE)
  cat("\nDone merging tiles\n")
}

# If CONUS or EAST, remove tile files, but first check that merged files 
# exist for each type (e.g., Lifestage, NumGen, ...)
if (region_param %in% c("CONUS", "EAST", "CONUSPLUS", "LOCO")) {
  cat("\nDeleting tiles for", region_param, "\n\n", 
      file = Model_rlogging, append = TRUE)
  cat("\nDeleting tiles for", region_param, "\n")
  for (t in type_list) {
    # Should only be a single merged file for these types
    if (grepl("Stress_Excl|Stress_Units|DDtotal", t)) {
      fls <- list.files(pattern = glob2rx(paste0("*", t, "_*all*tif$")))
      if (length(fls == 1)) {
        unlink(list.files(pattern = glob2rx(paste0("*", t, "_*tile*tif$"))))
      }
      # For other types, the number of merged files should equal the number 
      # of cohorts. The exception is if OW stage DD distro parameters for 
      # cohorts results in some cohorts not making DD cutoffs for PEMs - 
      # a warning will be recorded if this is the case
    } else {
      fls <- list.files(pattern = glob2rx(paste0("*", t, "_*all*tif$")))
      if (length(fls) == ncohort) {
        unlink(list.files(pattern = glob2rx(paste0("*", t, "_*tile*tif$"))))
        cat("Deleted tile files for", t, "\n", 
            file = Model_rlogging, append = TRUE)
        # Generate a warning message if some cohorts are missing
      } else {
        unlink(list.files(pattern = glob2rx(paste0("*", t, "_*tile*tif$"))))
        cat("\nWarning: only", length(fls), 
            "cohort files found for", t, "- check params\n", 
            file = Model_rlogging, append = TRUE)
        cat("\nWarning: only", length(fls), 
            "cohort files found for", t, "- check params\n")
      }
    }
  }
  cat("\nDone deleting tile files\n", file = Model_rlogging, append = TRUE)
  cat("\nDone deleting tile files\n")
}

#### * Some ggplot2 settings ####
# Some ggplot2 settings are specified here, but other settings are specified
# in the functions file

# Map production in ggplot requires specifying plot.height and plot.width
# These need to be dynamic because regions have different aspect ratios, 
# which results warped looking maps
# Calculate bounding box (xmin, xmax, ymin, ymax) of REGION 
coord <- coord_quickmap(xlim = c(REGION@xmin, REGION@xmax), 
                        ylim = c(REGION@ymin, REGION@ymax), expand = FALSE)
asp <- coord$aspect(list(x.range = c(REGION@xmin, REGION@xmax), 
                         y.range = c(REGION@ymin, REGION@ymax))) # aspect ratio

# Adjust base_size for ggplot2 (font size) according to aspect ratio
if (asp >= 1.7) {
  base_size <- 10.5
  legend_units <- 1.4
} else if (asp >= 1.5 & asp < 1.7) {
  base_size <- 9.5
  legend_units <- 1.3
} else if (asp >= 1.2 & asp < 1.5) {
  base_size <- 8.5 
  legend_units <- 1.2
} else if (asp >= 1 & asp < 1.2) {
  base_size <- 8
  legend_units <- 1
} else if (asp < 1 & asp >= 0.6) {
  base_size <- 7.5
  legend_units <- asp
} else if (asp < 0.6 & asp >= 0.45) {
  base_size <- 5.5
  legend_units <- asp
} else if (asp < 0.45) {
  base_size <- 5
  legend_units <- asp
}

# Theme to use for plots
mytheme <- theme(legend.text = element_text(size = rel(1)), 
                 legend.title = element_text(size = rel(1.2), face = "bold"),
                 legend.position = "right", 
                 legend.justification = "left",
                 legend.margin = margin(t = 0, r = 0.10, b = 0, l = 0.10, unit = "cm"),
                 legend.key.width = unit(legend_units, "line"), 
                 legend.key.height = unit(legend_units, "line"),
                 plot.title = element_text(size = rel(1.55), face = "bold", hjust = 0.5, 
                                           vjust = -3, lineheight = 1, 
                                           margin = margin(t = 0, r = 0, b = 2, l = 0)), 
                 plot.subtitle = element_text(size = rel(1.25), hjust = 0.5, vjust = -3, 
                                              lineheight = 1, 
                                              margin = margin(t = 5, r = 0, b = 15, l = 0)),
                 plot.margin = margin(t = 0.05, r = 0.25, b = 0.05, l = 0.25, unit = "cm"),
                 panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                 panel.background = element_blank(), panel.border = element_blank(),
                 axis.title.x = element_blank(), 
                 axis.title.y = element_blank(), 
                 axis.ticks = element_blank(),
                 axis.text.x = element_blank(), 
                 axis.text.y = element_blank())

### * DDtotal and climate stress ####

if (region_param %in% c("CONUS", "EAST", "CONUSPLUS", "LOCO")) {
  DDtotal_brick <- brick("DDtotal_cohort1_all.tif")
} else {
  DDtotal_brick <- brick("DDtotal_cohort1.tif")
}

if (exclusions_stressunits) {
  cat("\n", str_wrap("### SUMMARY MAPS: DDTOTAL, CLIMATE STRESS EXCL., AND 
                     CLIMATE STRESS UNITS ###", width = 80), sep = "", 
      file = Model_rlogging, append = TRUE)
  cat("\n", str_wrap("SUMMARY MAPS: DDTOTAL, CLIMATE STRESS EXCL., AND CLIMATE 
               STRESS UNITS", width = 80), "\n", sep = "")
} else {
  cat("\n### SUMMARY MAPS: DDTOTAL ###", file = Model_rlogging, append = TRUE)
  cat("\nSUMMARY MAPS: DDTOTAL\n")
}

# Split up the dates into chunks - this avoids overloading the server w/ 
# running too many dates in parallel
dats_list <- split(dats2, ceiling(1:length(dats2)/4))
last_date <- last(dats2)

# For each date in a date chunk, plot and save summary maps for:
# degree-day accumulation, chill stress unit accumulation, chill stress 
# exclusion, heat stress unit accumulation, heat stress exclusion, and all 
# stress exclusion
RegCluster(ncohort)

#for (dat in dats_list) {
stress_results <- foreach(dat = dats_list, .packages = pkgs, 
                          .inorder = TRUE) %dopar% {
                            dat_vec <- unname(unlist(dat)) # change to an unnamed vector
                            for (d in dat_vec) {
                              # get position (layer) of date in raster brick
                              lyr <- which(dats2 == d)
                              # make the plots
                              PlotMap(DDtotal_brick[[lyr]],d, "Degree day (DD) accumulation", 
                                      "Degree Days", "Misc_output/DDtotal")
                              
                              if (exclusions_stressunits) {
                                # Bring in climate stress bricks for each cohort
                                # Chill stress unit accumulation
                                chillunitsCUM_patrn <- glob2rx("*Chill_Stress_Units*1*.tif$")
                                chillunitsCUM_brick <- brick(list.files(pattern = chillunitsCUM_patrn))
                                PlotMap_stress(chillunitsCUM_brick[[lyr]], d, chillstress_units_max1,
                                               chillstress_units_max2, "Chill stress units", 
                                               "Chill Stress Units", "Misc_output/Chill_Stress_Units")
                                # Chill stress exlusions (-1 = moderate; -2 = severe)
                                chillEXCL_patrn <- glob2rx("*Chill_Stress_Excl*1*.tif$")
                                chillEXCL_brick <- brick(list.files(pattern = chillEXCL_patrn))
                                PlotMap(chillEXCL_brick[[lyr]], d,  "Chill stress exclusion", 
                                        "Exclusion status", "Misc_output/Chill_Stress_Excl")
                                # Heat unit accumulation
                                heatunitsCUM_patrn <- glob2rx("*Heat_Stress_Units*1*.tif$")
                                heatunitsCUM_brick <- brick(list.files(pattern = heatunitsCUM_patrn))          
                                PlotMap_stress(heatunitsCUM_brick[[lyr]], d, heatstress_units_max1,
                                               heatstress_units_max2, "Heat stress units", 
                                               "Heat Stress Units", "Misc_output/Heat_Stress_Units")
                                # Heat stress exclusions (-1 = moderate; -2 = severe)
                                heatEXCL_patrn <- glob2rx("*Heat_Stress_Excl*1*.tif$")
                                heatEXCL_brick <- brick(list.files(pattern = heatEXCL_patrn))
                                PlotMap(heatEXCL_brick[[lyr]], d,  "Heat stress exclusion", 
                                        "Exclusion status", "Misc_output/Heat_Stress_Excl")
                                # All stress exclusions (chill stress + heat stress exclusions)
                                AllEXCL_patrn <- glob2rx("*All_Stress_Excl*1*.tif$")
                                AllEXCL_brick <- brick(list.files(pattern = AllEXCL_patrn))
                                PlotMap(AllEXCL_brick[[lyr]], d,  "All stress exclusion", 
                                        "Exclusion status", "All_Stress_Excl")
                              }
                            }
                          }

stopCluster(cl)

# Log file messages

# If no PEMS and no climate stress exclusions, then moving on to Lifestage 
# analyses 
if (!pems & !exclusions_stressunits) {
  cat("\n\nDone with DDtotal maps\n", file = Model_rlogging, append = TRUE)  
  cat("\n", str_wrap("### SUMMARY MAPS AND WEIGHTED RASTER OUTPUT: LIFESTAGE 
                     ###", width = 80), "\n\nStages: ", 
      paste(stgorder, collapse = ", "), 
      sep = "", file = Model_rlogging, append = TRUE)
  cat("\nSUMMARY MAPS AND WEIGHTED RASTER OUTPUT: LIFESTAGE ###","\n\nStages: ",
      paste(stgorder, collapse = ", "), "\n")
  # If no PEMS but climate stress exclusions, then moving on to Lifestage 
  # analyses that also include climate stress exclusions
} else if (!pems & exclusions_stressunits) {
  cat("\n\n", str_wrap("Done with DDtotal, climate stress exclusions, and 
                       climate stress unit maps", width = 80), "\n", sep = "",
      file = Model_rlogging, append = TRUE)  
  cat("\n", str_wrap("### SUMMARY MAPS AND WEIGHTED RASTER OUTPUT:
      LIFESTAGE W/ CLIM. STRESS EXCL. ###", width = 80), "\n\nStages: ", 
      paste(stgorder, collapse = ", "), sep = "",
      file = Model_rlogging, append = TRUE)
  cat("\n", str_wrap("### SUMMARY MAPS AND WEIGHTED RASTER OUTPUT:
      LIFESTAGE W/ CLIM. STRESS EXCL. ###", width = 80), "\nStages: ", 
      paste(stgorder, collapse = ", "), sep = "")
  # If PEMS but no climate stress exclusions, then conducting PEM analyses
} else if (pems & !exclusions_stressunits) {
  cat("\n\nDone with DDtotal summary maps\n\n", 
      "### SUMMARY MAPS AND RASTER OUTPUT: PEST EVENT MAPS ###", sep = "",
      file = Model_rlogging, append = TRUE)
  cat("\nDone with DDtotal summary maps\n\n", 
      "SUMMARY MAPS AND RASTER OUTPUT: PEST EVENT MAPS\n", sep = "")
  # If PEMS and climate stress exclusions, then conducting PEM analyses that 
  # also include climate stress exclusions
} else if (pems & exclusions_stressunits) {
  cat("\n\n", str_wrap("Done with DDtotal, climate stress exclusions, and 
                       climate stress unit maps", width = 80), "\n\n", 
      str_wrap("### SUMMARY MAPS AND RASTER OUTPUT: PEST EVENT MAPS W/
                       CLIMATE STRESS EXCL. ###", width = 80), sep = "",
      file = Model_rlogging, append = TRUE)
  cat("\n", str_wrap("Done with DDtotal, climate stress exclusions, and climate 
               stress unit maps", width = 80), "\n\n", 
      str_wrap("SUMMARY MAPS AND RASTER OUTPUT: PEST EVENT MAPS W/
                CLIMATE STRESS EXCL.\n", width = 80), sep = "")
}




#### * Weight bricks for Diapause/Mismatch; save rasters and summary plots ####
# beware that AttVolt, Voltinism, and Diapause need to be /1000 first
if(do_photo){
  if (region_param %in% c("CONUS", "EAST", "CONUSPLUS", "LOCO")) {
    template <- merge(template[[1]], template[[2]], template[[3]], template[[4]])
  }
  
  
  FullGen_fls <- list.files(pattern = glob2rx("*FullGen_*.tif$"))
  AttVolt_fls <- list.files(pattern = glob2rx("*AttVolt_*.tif$"))
  # Voltinism_fls <- list.files(pattern = glob2rx("*Voltinism_*.tif$"))
  Diapause_fls <- list.files(pattern = glob2rx("*Diapause_*.tif$"))
  
  # Calculate the highest generation to occur across all FullGen bricks
  maxgens <- max(unlist(lapply(FullGen_fls, 
                               function(x) { maxValue(brick(x)) })))
  maxgens2 <- max(unlist(lapply(AttVolt_fls, 
                                function(x) { maxValue(brick(x)/1000) })))
  maxgens <- ceiling(round(maxgens, maxgens2))
  
  # Create derived rasters (Mismatch) for each cohort before weighting
  RegCluster(ncohort)
  FullGen_wtd_byGen <- foreach(coh = as.list(1:ncohort), 
                               .packages = pkgs, .inorder = TRUE) %dopar% {
                                 Mismatch <- 1000 * (brick(AttVolt_fls[coh])/1000 - brick(FullGen_fls[coh]))
                                 SaveRaster2(Mismatch, paste0("Mismatch_cohort", coh), "INT2S", "Mismatch")
                               }
  stopCluster(cl)
  
  Mismatch_fls <- list.files(pattern = glob2rx("*Mismatch_*.tif$"))
  
  
  # Weight diapause stuff
  to_weight <- list(FullGen_fls, AttVolt_fls, Diapause_fls, Mismatch_fls)
  RegCluster(ncohort)
  Diap_wtd_fls <- foreach(index = 1:4, 
                          .packages = pkgs, 
                          .inorder = FALSE) %dopar% {
                            fls <- to_weight[[index]]
                            if(index == 1){ # FullGen not multiplied by 1000 like others
                              Ras_weighted <- mapply(function(x,y) { 
                                round(brick(x) * y * 1000)
                              }, x=fls, y=relpopsize)
                            }else{
                              Ras_weighted <- mapply(function(x,y) { 
                                round(brick(x) * y)
                              }, x=fls, y=relpopsize)
                            }
                            Ras_weighted_sum <- Reduce("+", Ras_weighted)
                            SaveRaster2(Ras_weighted_sum, 
                                        paste(c("FullGen", "AttVolt", "Diapause", "Mismatch")[index], 
                                              "all", "weighted", sep="_"), 
                                        c("INT2U", "INT2U", "INT2U", "INT2S")[index], " ")
                            return(paste0(paste(c("FullGen", "AttVolt", "Diapause", "Mismatch")[index], 
                                                "all", "weighted", sep="_"), ".tif"))
                          }
  stopCluster(cl)
  
  
  ### * Create summary maps of Diapause results, weighted across cohorts
  cat("\n### SUMMARY MAP OUTPUT - Diapause ###\n\n", file=Model_rlogging, 
      append=TRUE)
  cat("\nSUMMARY MAP OUTPUT - Diapause\n")
  
  # Create list of raster brick files 
  Diap_wtd_fls_nams <- Diap_wtd_fls %>% gsub(".tif","",.)
  names(Diap_wtd_fls) <- Diap_wtd_fls_nams
  
  Diap_wtd_mrgd_brk <- brick() # blank brick to put named raster into
  j <- 1
  for (Diap_wtd_fl in Diap_wtd_fls) {
    Diap_wtd_brk <- lapply(Diap_wtd_fl, function(x) { 
      brk <- brick(x)
      type <- unique(str_split_fixed(names(brk), pattern = "_", 3)[,1])
      names(brk) <- paste(type,dats2,sep="_")
      return(brk)
    })
    Diap_wtd_mrgd_brk <- addLayer(Diap_wtd_mrgd_brk, Diap_wtd_brk)
  }
  
  RegCluster(ncohort)
  Diap_sum_maps <- foreach(index = 1:nlayers(Diap_wtd_mrgd_brk), 
                           .packages = pkgs, .inorder = TRUE) %dopar%{
                             # print(index)
                             brk <- Diap_wtd_mrgd_brk[[index]] 
                             nam <- unique(str_split_fixed(names(brk), pattern = "_", 2)[,1])
                             # print(nam)
                             dat <- unique(str_split_fixed(names(brk), pattern = "_", 2)[,2])
                             # print(dat)
                             brk <- brk + template # adding template ensures NA's are NA's instead of 0
                             df <- ConvDF(brk)
                             if(all(df$value == 0)){
                               nam <- "skip"
                             }  #probably not interesting to plot
                             
                             sp <- paste0(fullname,":")
                             nicedat <- as.character(format(strptime(dat,format="%Y%m%d"), 
                                                            format="%m/%d/%Y"))
                             titl <- case_when(nam == "AttVolt" ~ "Attempted voltinism",
                                               nam == "FullGen" ~ "Potential voltinism",
                                               nam == "Diapause" ~ "Chose to diapause",
                                               nam == "Mismatch" ~ "Voltinism Mismatch")
                             titl <- paste(titl, nicedat, sep=" ")
                             subtitl <- paste("Maps and modeling",format(Sys.Date(),"%m/%d/%Y"),
                                              "by Oregon IPM Center")
                             
                             
                             if(nam %in% c("AttVolt", "FullGen")){
                               
                               plotmax <- max(df$value / 1000)
                               
                               if (plotmax < 4){
                                 df$value <- as.factor(round(df$value / 1000 / .5) * .5)
                               }
                               # if (plotmax >= 5){
                               #   df$value[df$value >= 5000] <- 5000
                               #   df$value <- as.factor(round(df$value / 1000 / .5) * .5)
                               #   levels(df$value)[length(levels(df$value))] <- 
                               #     paste0(levels(df$value)[length(levels(df$value))], " or more")
                               #   plotmax <- 5
                               # }
                               if (plotmax >= 4){
                                 df$value[df$value >= 5000] <- 5000
                                 df$value <- as.factor(round(df$value / 1000))
                                 levels(df$value)[length(levels(df$value))] <- 
                                   paste0(levels(df$value)[length(levels(df$value))], " or more")
                                 plotmax <- 5
                               }
                               
                               p <- Base_map(df) + 
                                 scale_fill_viridis(discrete = TRUE, name = "Generations", 
                                                    begin = 0, end = plotmax/5) +
                                 # geom_point(data = sites, aes(x = x, y = y), color = "black", size = 3) +
                                 labs(title = str_wrap(paste(sp,titl), width = 55), 
                                      subtitle = str_wrap(subtitl, width = 75)) +
                                 theme_map(base_size = base_size) + mytheme
                               ggsave(p,file=paste0(nam,"_", dat, ".png"), height = asp * 7, 
                                      units = c('in'), dpi = 300) 
                             } 
                             
                             if(nam == "Diapause"){
                               # TODO Diapause map doesn't match lost generations because of pre-diap?
                               df$value <- df$value / 10 # Percent
                               
                               p <- Base_map(df) + 
                                 scale_fill_viridis(discrete = FALSE, name = "% in Diapause", 
                                                    begin = 0, end = 1) +
                                 labs(title = str_wrap(paste(sp,titl), width = 55), 
                                      subtitle = str_wrap(subtitl, width = 75)) +
                                 theme_map(base_size = base_size) + mytheme
                               ggsave(p,file=paste0(nam,"_", dat, ".png"), height = asp * 7, 
                                      units = c('in'), dpi = 300) 
                             } 
                             
                             if(nam == "Mismatch"){
                               minmm <- round(min(df$value/1000), 1) - .1
                               df$value <- df$value / 1000
                               if(minmm < -4){
                                 mmbrk <- c(minmm, -4, -3, -2, -1, 0, .25, .5, .75, 1)
                                 df$value <- cut(df$value, mmbrk)
                                 cols <- setNames(c("#8c510a", "#bf812d", "#dfc27d", "#f6e8c3", "#f5f5f5", 
                                                    "#c7eae5", "#80cdc1", "#35978f", "#01665e"), 
                                                  levels(df$value))
                               }else{
                                 mmbrk <- c(-4, -3, -2, -1, 0, .25, .5, .75, 1)
                                 df$value <- cut(df$value, mmbrk)
                                 cols <- setNames(c("#bf812d", "#dfc27d", "#f6e8c3", "#f5f5f5", 
                                                    "#c7eae5", "#80cdc1", "#35978f", "#01665e"), 
                                                  levels(df$value))
                               }
                               
                               p <- Base_map(df) + 
                                 # geom_point(data = sites, aes(x = x, y = y), color = "black", size = 3) +
                                 scale_fill_manual(values = cols, name = "Mismatch") +
                                 labs(title = str_wrap(paste(sp,titl), width = 55), 
                                      subtitle = str_wrap(subtitl, width = 75)) +
                                 theme_map(base_size = base_size) + mytheme
                               ggsave(p,file=paste0(nam,"_", dat, ".png"), height = asp * 7, 
                                      units = c('in'), dpi = 300) 
                             } 
                           }
  stopCluster(cl)
} # close do_photo


#### * Analyses and map production all done - wrap-up ####
processing_exectime <- toc(quiet = TRUE)
processing_exectime <- (processing_exectime$toc - processing_exectime$tic) / 60 

cat("\n### Done w/ final analyses and map production ###\n", 
    "Run time for analyses and map production = ", 
    round(processing_exectime, digits = 2), " min\n", sep = "", 
    "Deleting, renaming, and moving some remaining files\n",
    file = Model_rlogging, append = TRUE)
cat("\nDone w/ final analyses and map production\n\n", 
    "Run time for analyses and mapping run time = ", 
    round(processing_exectime, digits = 2),
    " min\n\n", "Deleting, renaming, and moving some remaining files\n\n", 
    sep = "")

# Delete all output files from daily loop now that they have been processed
if (region_param %in% c("CONUS", "EAST", "CONUSPLUS", "LOCO")) {
  unlink(list.files(pattern = glob2rx(paste0("*1_all.tif$"))))
} else {
  unlink(list.files(pattern = glob2rx(paste0("*_cohort.tif$"))))
}

# Rename files for last day of year
last_dat_fls <- list.files(pattern = glob2rx(paste0("*", last_date, "*.png$")))
new_names <- paste0(spp, "_", last_dat_fls)
invisible(file.rename(last_dat_fls, new_names))
cat("Renamed all files for last day (", last(dats2), ") to include ", spp, 
    " in file name\n", sep = "", file = Model_rlogging, append = TRUE)
cat("\nRenamed all files for last day (", last(dats2), ") to include ", spp, 
    " in file name\n", sep = "")

# Move all leftover files without "LBAM" in file name to "Misc_output" 
# (grep finds the inverse of pattern here)
misc_fls <- grep(list.files(path = output_dir), 
                 pattern = glob2rx(paste0("*", last_date, "*.png$")), 
                 invert = TRUE, value = TRUE)
invisible(file.copy(misc_fls, paste0(output_dir, "/Misc_output/")))
invisible(file.remove(misc_fls))

# Wrap up log file and report run time for entire model
cat("\nMODEL RUN DONE\n", file = Model_rlogging, append = TRUE)
cat("\nMODEL RUN DONE\n")
total_exectime <- toc(quiet = TRUE) # Execution time for entire run
total_exectime <- round((total_exectime$toc - total_exectime$tic) / 60, 
                        digits = 2)
cat("Run time for entire model =", total_exectime, "min", 
    file = Model_rlogging, append = TRUE)
cat("\nRun time for entire model =", total_exectime, "min\n\n")
