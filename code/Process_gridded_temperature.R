## HEADER -----
## Script name: Process_gridded_temperature.R
## Purpose of script: Convert downloaded data to daily min/max temperature grids
## Author: Tyson Wepprich
## Date Created: 2021-02-05
## License: GNU GPLv3
## Email: tyson.wepprich@gmail.com
## ---
## Notes:
## We processed data downloaded from Daymet and E-OBS to make daily temperature
## min/max rasters for easier input into DDRP framework. Alternatively, DDRP could
## be recoded to accept .nc files or raster bricks.
## ---


# Process Daymet ----
# Downloaded 2019 version 3 tmin and tmax .nc files from 
# https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=1328
# Put .nc files into data/daymet directory 
library(raster)
library(dplyr)
library(ncdf4)
library(foreach)
library(doParallel)
library(stringr)

# Steps
# 1. Continent raster brick per year in .nc file
# 2. Process in parallel for each day (takes 40 minutes with 4 cores on laptop)
# 3. projectRaster by day slice and write daily file for DDRP

new_dir <- "data/daymet"

cl <- makePSOCKcluster(4)
registerDoParallel(cl)

agg_res <- 4 # aggregated 1-km resolution data to 4-km for computational speed
fs <- list.files("data/daymet", pattern = '.nc4', full.names = TRUE, recursive = TRUE)

outfiles <- foreach(f = 1:length(fs),
                    .packages= c("raster", "stringr"),
                    .inorder = FALSE)%:%
  foreach(day = 1:365,
          .packages= c("raster", "stringr"),
          .inorder = FALSE)%dopar%{
            
            fl <- fs[f]
            newname <- stringr::str_split_fixed(string = fl, pattern = "_", n = 3)[,3]
            varname <- stringr::str_split_fixed(string = newname, pattern = "_", n = 3)[,1]
            yr <- stringr::str_split_fixed(string = newname, pattern = "_", n = 3)[,2]
            
            ifelse(!dir.exists(file.path(new_dir, yr)), dir.create(file.path(new_dir, yr), recursive = TRUE), FALSE)
            
            
            newname <- stringr::str_split_fixed(string = newname, pattern = "_", n = 3)[,1:2]
            # for CONUSPLUS
            e_2 <- extent(-131, -63, 24, 55)
            e_1 <- extent(-3146358, 4689822, -2035046, 2219754)
            
            fslice <- brick(fl, varname = varname, level = 1)[[day]]
            proj4string(fslice) <- CRS("+proj=lcc +lon_0=-100 +lat_0=42.5 +x_0=0 +y_0=0 +lat_1=25 +lat_2=60 +datum=WGS84 +units=m")
            
            # reproject to prism CRS
            prism_res <- 25.875 / 621 # divided by 4 for Daymet 1 x 1
            newproj <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
            
            agg <- aggregate(crop(fslice, e_1), fact=agg_res, fun=mean, expand=FALSE, na.rm=TRUE)
            
            outr <- projectRaster(agg, res = prism_res, crs = newproj, method="bilinear")
            
            out2 <- crop(outr, e_2)
            
            dat <- gsub(pattern = "X", replacement = "", x = names(fslice), fixed = TRUE) 
            dat <- gsub(pattern = ".", replacement = "", x = dat, fixed = TRUE) 
            
            fname <-  paste0(getwd(), "/", new_dir, "/", newname[2], "/DAYMET_", varname, "_stable_4kmD1_", dat, ".grd")
            writeRaster(out2, filename = fname, overwrite = TRUE)
            
            removeTmpFiles(h=0.1)
          }

if(exists("cl")){
  stopCluster(cl)
}

# check files
fs <- list.files(new_dir, pattern = '.grd', full.names = TRUE, recursive = TRUE)



# Process E-OBS ----
# Downloaded 2011-2020 version 22.0e TN and TX .nc files from 
# https://surfobs.climate.copernicus.eu/dataaccess/access_eobs_chunks.php
# Put .nc files into data/eobs directory 
# For manuscript 1995-2010 data downloaded to include 2010, too.

# longitude in degrees EAST
# time is days since 1950-01-01 00:00:00
# temperature in C
setwd("data/eobs")
targetdir <- getwd()

# Function to split 5-year data into daily tmin/tmax rasters for use in lifecycle model
SplitEOBS <- function(eobsfile, tvar, targetdir, ncores = 4){
  cl <<- makePSOCKcluster(ncores) # export to global environment
  registerDoParallel(cl)
  
  ras <- brick(eobsfile, varname = tvar, lvar = 3, level = 4)
  
  # if you have permission to write, you could make directories for each year here
  yrs <- gregexpr(pattern = "[0-9]{4}", text = eobsfile)
  yr1 <- as.numeric(substr(eobsfile, start = yrs[[1]][1], stop = yrs[[1]][1] + 3))
  yr2 <- as.numeric(substr(eobsfile, start = yrs[[1]][2], stop = yrs[[1]][2] + 3))
  # var <- sub(pattern = "as", replacement = "", x = stringr::str_split_fixed(eobsfile, pattern = "_", n = 3)[, 2], fixed = TRUE)
  
  for (y in yr1:yr2){
    if (dir.exists(paste(targetdir, y, sep = "/")) == FALSE){
      dir.create(paste(targetdir, y, sep = "/"))
    }
  }
  
  fname <- stringr::str_split_fixed(eobsfile, pattern = "_", 9)
  varname <- paste(ifelse(tvar == "tn", "tmin", "tmax"), fname[2], fname[3], fname[4], sep = "_")
  
  indices <- 1:nlayers(ras)
  # loop every day and write new daily rasters
  foreach(d = indices, .packages = c('raster', "ncdf4", "stringr")) %dopar% {
    dras <- ras[[d]]
    date <- gsub(pattern = "[^0-9]", replacement = "", x = names(dras))
    yr <- substr(date, start = 1, stop = 4)
    newname <- paste0("EOBS_", varname, "_", date, ".grd")
    
    writeRaster(x = dras, filename = paste(targetdir, yr, newname, sep = "/"), format = "raster", overwrite = TRUE)
  }
  
  
  stopCluster(cl)
  
}

# Run this function to split .nc into 10 years' daily rasters for DDRP
# I only keep 2019 for this repository's example
# Takes 15 minutes with 12 cores on a laptop
eobsfile <- "tn_ens_mean_0.1deg_reg_2011-2020_v22.0e.nc"
tvar <- "tn"
SplitEOBS(eobsfile = "tn_ens_mean_0.1deg_reg_2011-2020_v22.0e.nc", 
          tvar = "tn", targetdir, ncores = 6)
SplitEOBS(eobsfile = "tx_ens_mean_0.1deg_reg_2011-2020_v22.0e.nc", 
          tvar = "tx", targetdir, ncores = 12)


# check files
fs <- list.files("data/eobs/2019", pattern = '.grd', full.names = TRUE, recursive = TRUE)

