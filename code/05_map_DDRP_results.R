## Header ----
## Script name: 05_map_DDRP_results.R
## Purpose of script: Map results from DDRP annual lifecycle model
## Author: Tyson Wepprich
## Date Created: 2021-02-04
## License: GNU GPLv3
## Email: tyson.wepprich@gmail.com
## ---
## Notes:
## Script makes maps in manuscript from Europe and North America
## Can take multiple years of DDRP results to average (2010-2019) or
## run a single year's results. 
## ---
# Load the required packages

pkgs <- c("colorspace", "dplyr", "ggplot2", "ggthemes", 
          "lubridate", "mapdata", 
          "purrr", "RColorBrewer", "rgdal", "raster", "readr", "sp", "stringr", 
          "tidyr", "tools", "viridis")
lapply(pkgs, library, character.only = TRUE)
source("code/DDRP_cohorts_v1_funcs.R")

theme_set(theme_bw(base_size = 12) +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())) 

# Function to average years or map single year results
# In this manuscript, for parameters:
# region is EUROPE or CONUSPLUS
# years is 2010:2019 or single years of warmest/coolest conditions within decade
# maptype is FullGen, AttVolt, Mismatch or Diapause
# -FullGen is potential voltinism based on degree-days alone
# -AttVolt is attempted voltinism based on degree-days and photoperiod response
# -Mismatch is Attempted - Potential voltinism
# -Diapause is the percent of the population choosing diapause upon photo-sensitivity
# biotype is Kyushu or Hokkaido
# filename matches manuscript figures

MultiyearMaps <- function(region, years, maptype, biotype, filename){
  nam <- maptype
  REGION <- Assign_extent(region) # Bounding box
  if (region == "EUROPE") filereg <- "EURO"
  if (region == "CONUSPLUS") filereg <- "NA"
  
  # Create list of raster brick files 
  Diap_wtd_fls <- list.files(path = "data/DDRP_results", pattern = "all_weighted.tif", full.names = TRUE, recursive = TRUE)
  Diap_wtd_fls <- Diap_wtd_fls[grep(pattern = "APH2", x = Diap_wtd_fls)]
  Diap_wtd_fls <- Diap_wtd_fls[grep(pattern = filereg, x = Diap_wtd_fls)]
  fls_10yr <- Diap_wtd_fls[grep(pattern = nam, x = Diap_wtd_fls, fixed = TRUE)]
  fls_10yr <- fls_10yr[grep(pattern = biotype, x = fls_10yr, fixed = TRUE)]
  fls <- lapply(years, FUN = function(x){
    return(fls_10yr[grep(pattern = x, x = fls_10yr)])
  })
  
  brks <- lapply(X = fls, FUN = function(x) brick(x)[[6]]) # for bimonthly maps, take last result
  if (length(years) == 1){
    brk <- stack(brks[[1]])
    yr <- as.character(years)
  }
  if (length(years) > 1){
    brk <- mean(stack(brks), na.rm = TRUE)
    yr <- paste(years[1], "-", years[length(years)], " (average)")
  }
  
  df <- ConvDF(brk)
  # nam <- "FullGen"
  titl <- paste(case_when(nam == "AttVolt" ~ "Attempted voltinism",
                          nam == "FullGen" ~ "Potential voltinism",
                          nam == "Diapause" ~ "Chose to diapause",
                          nam == "Mismatch" ~ "Voltinism mismatch"), "for", yr, sep = " ")
  if(nam == "FullGen"){
    subtitl <- "No photoperiod response"
  }else{
    subtitl <- paste0(biotype, " photoperiod response")
  }
  
  base_size <- 14
  
  if (region == "EUROPE") reg.df <- map_data("world"); map.width <- 6
  if (region == "CONUSPLUS") reg.df <- readRDS("data/na_lines.rds"); map.width <- 10
  if (region == "CONUS") reg.df <- map_data("states"); map.width <- 10
  
  
  
  if(nam %in% c("AttVolt", "FullGen")){
    
    # too many generations is overwhelming on map 
    plotmax <- max(df$value / 1000)
    df$value[df$value >= 6000] <- 6000
    df$value <- as.factor(round(df$value / 1000 / .5) * .5) # round to nearest 1/2 generation
    # df$value <- as.factor(round(df$value / 1000)) # round to nearest integer generation
    if(plotmax > 6){
      levels(df$value)[length(levels(df$value))] <- 
        paste0(levels(df$value)[length(levels(df$value))], " or more")
    }
    plotmax <- 6
    
    
    p <- ggplot(reg.df, aes(x = long, y = lat)) + 
      geom_raster(data = df, aes(x = x, y = y, fill = value)) + 
      geom_path(aes(group = group), color = "gray20", lwd = 0.4) +
      coord_quickmap(xlim = c(REGION@xmin, REGION@xmax), 
                     ylim = c(REGION@ymin, REGION@ymax), expand = FALSE) +
      # scale_fill_viridis(discrete = TRUE, name = "Generations", begin = 0, end = plotmax/6) +
      scale_fill_manual(name = "Generations",
                        values = Turbo(out.colors = 13))+ ## This is the change
      # theme_map(base_size = base_size) +
      theme(legend.position = "right") +
      ggtitle(titl, subtitle = subtitl) +
      xlab("Longitude (°)") +
      ylab("Latitude (°)") +
      guides(colour=FALSE)
    
    ggsave(p, path = "figures", file = paste0(filename,"_", maptype, ".png"), 
           height = 6, width = map.width,
           units = c('in'), dpi = 300) 
  }
  
  if(nam == "Diapause"){
    df$value <- df$value / 10 # Percent
    
    p <- ggplot(reg.df, aes(x = long, y = lat)) + 
      geom_raster(data = df, aes(x = x, y = y, fill = value)) + 
      geom_path(aes(group = group), color = "gray20", lwd = 0.4) +
      coord_quickmap(xlim = c(REGION@xmin, REGION@xmax), 
                     ylim = c(REGION@ymin, REGION@ymax), expand = FALSE) +
      scale_fill_viridis(discrete = FALSE, name = "Diapause %", begin = 0, end = 1) +
      # theme_map(base_size = base_size) +
      theme(legend.position = "right") +
      ggtitle(titl, subtitle = subtitl) +
      xlab("Longitude (°)") +
      ylab("Latitude (°)") +
      guides(colour=FALSE)
    
    ggsave(path = "figures", file = paste0(filename,"_", maptype, ".png"), 
           height = 6, width = map.width,
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
    
    p <- ggplot(reg.df, aes(x = long, y = lat)) + 
      geom_raster(data = df, aes(x = x, y = y, fill = value)) + 
      geom_path(aes(group = group), color = "gray20", lwd = 0.4) +
      coord_quickmap(xlim = c(REGION@xmin, REGION@xmax), 
                     ylim = c(REGION@ymin, REGION@ymax), expand = FALSE) +
      scale_fill_manual(values = cols, name = "Mismatch") +
      theme(legend.position = "right") +
      ggtitle(titl, subtitle = subtitl) +
      xlab("Longitude (°)") +
      ylab("Latitude (°)") +
      guides(colour=FALSE)
    
    ggsave(p, path = "figures", file = paste0(filename,"_", maptype, ".png"), 
           height = 6, width = map.width,
           units = c('in'), dpi = 300) 
  } 
}

# Reproduce maps in manuscript

# In supplement, annual variation was from two extreme years within the decade
# 2017 coldest, 2018 warmest for Europe by degree-days
# 2013 coldest, 2012 warmest for North America by degree-days

# Figure 4
MultiyearMaps(region = "CONUSPLUS", years = c(2010:2019), maptype = "FullGen", biotype = "Kyushu", filename = "fig4a_NA_mean")
MultiyearMaps(region = "EUROPE", years = c(2010:2019), maptype = "FullGen", biotype = "Kyushu", filename = "fig4b_EURO_mean")

# Figure 5
MultiyearMaps(region = "CONUSPLUS", years = c(2010:2019), maptype = "AttVolt", biotype = "Hokkaido", filename = "fig5a_NA_mean")
MultiyearMaps(region = "CONUSPLUS", years = c(2010:2019), maptype = "AttVolt", biotype = "Kyushu", filename = "fig5c_NA_mean")
MultiyearMaps(region = "EUROPE", years = c(2010:2019), maptype = "AttVolt", biotype = "Hokkaido", filename = "fig5b_EURO_mean")
MultiyearMaps(region = "EUROPE", years = c(2010:2019), maptype = "AttVolt", biotype = "Kyushu", filename = "fig5d_EURO_mean")

# Figure 6
MultiyearMaps(region = "CONUSPLUS", years = c(2010:2019), maptype = "Mismatch", biotype = "Hokkaido", filename = "fig6a_NA_mean")
MultiyearMaps(region = "CONUSPLUS", years = c(2010:2019), maptype = "Mismatch", biotype = "Kyushu", filename = "fig6c_NA_mean")
MultiyearMaps(region = "EUROPE", years = c(2010:2019), maptype = "Mismatch", biotype = "Hokkaido", filename = "fig6b_EURO_mean")
MultiyearMaps(region = "EUROPE", years = c(2010:2019), maptype = "Mismatch", biotype = "Kyushu", filename = "fig6d_EURO_mean")

# Figure S2
MultiyearMaps(region = "CONUSPLUS", years = 2012, maptype = "FullGen", biotype = "Kyushu", filename = "figS2a_NA_2012_warmest")
MultiyearMaps(region = "CONUSPLUS", years = 2013, maptype = "FullGen", biotype = "Kyushu", filename = "figS2c_NA_2013_coolest")
MultiyearMaps(region = "EUROPE", years = 2018, maptype = "FullGen", biotype = "Kyushu", filename = "figS2b_EURO_2018_warmest")
MultiyearMaps(region = "EUROPE", years = 2017, maptype = "FullGen", biotype = "Kyushu", filename = "figS2d_EURO_2017_coolest")

# Figure S3
MultiyearMaps(region = "CONUSPLUS", years = 2012, maptype = "AttVolt", biotype = "Hokkaido", filename = "figS3a_NA_2012_warmest")
MultiyearMaps(region = "CONUSPLUS", years = 2013, maptype = "AttVolt", biotype = "Hokkaido", filename = "figS3c_NA_2013_coolest")
MultiyearMaps(region = "EUROPE", years = 2018, maptype = "AttVolt", biotype = "Hokkaido", filename = "figS3b_EURO_2018_warmest")
MultiyearMaps(region = "EUROPE", years = 2017, maptype = "AttVolt", biotype = "Hokkaido", filename = "figS3d_EURO_2017_coolest")

# Figure S4
MultiyearMaps(region = "CONUSPLUS", years = 2012, maptype = "AttVolt", biotype = "Kyushu", filename = "figS4a_NA_2012_warmest")
MultiyearMaps(region = "CONUSPLUS", years = 2013, maptype = "AttVolt", biotype = "Kyushu", filename = "figS4c_NA_2013_coolest")
MultiyearMaps(region = "EUROPE", years = 2018, maptype = "AttVolt", biotype = "Kyushu", filename = "figS4b_EURO_2018_warmest")
MultiyearMaps(region = "EUROPE", years = 2017, maptype = "AttVolt", biotype = "Kyushu", filename = "figS4d_EURO_2017_coolest")

# Figure S5
MultiyearMaps(region = "CONUSPLUS", years = 2012, maptype = "Mismatch", biotype = "Hokkaido", filename = "figS5a_NA_2012_warmest")
MultiyearMaps(region = "CONUSPLUS", years = 2013, maptype = "Mismatch", biotype = "Hokkaido", filename = "figS5c_NA_2013_coolest")
MultiyearMaps(region = "EUROPE", years = 2018, maptype = "Mismatch", biotype = "Hokkaido", filename = "figS5b_EURO_2018_warmest")
MultiyearMaps(region = "EUROPE", years = 2017, maptype = "Mismatch", biotype = "Hokkaido", filename = "figS5d_EURO_2017_coolest")

# Figure S6
MultiyearMaps(region = "CONUSPLUS", years = 2012, maptype = "Mismatch", biotype = "Kyushu", filename = "figS6a_NA_2012_warmest")
MultiyearMaps(region = "CONUSPLUS", years = 2013, maptype = "Mismatch", biotype = "Kyushu", filename = "figS6c_NA_2013_coolest")
MultiyearMaps(region = "EUROPE", years = 2018, maptype = "Mismatch", biotype = "Kyushu", filename = "figS6b_EURO_2018_warmest")
MultiyearMaps(region = "EUROPE", years = 2017, maptype = "Mismatch", biotype = "Kyushu", filename = "figS6d_EURO_2017_coolest")
