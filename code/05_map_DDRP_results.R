# Multiyear maps

# Load the required packages
pkgs <- c("doParallel", "plyr", "dplyr", "foreach", "ggplot2", "ggthemes", 
          "lubridate", "mapdata", "mgsub", "optparse", "parallel",
          "purrr", "RColorBrewer", "rgdal", "raster", "readr", "sp", "stringr", 
          "tidyr", "tictoc", "tools", "viridis")
lapply(pkgs, library, character.only = TRUE)
source("DDRP_cohorts_v1_funcs_euro.R")
source("turbo.R")
params_dir <- "/home/tyson/REPO/ddrp-cohorts-v1/spp_params/" # tyson's GRUB
REGION <- Assign_extent("EUROPE") # Bounding box
theme_set(theme_bw(base_size = 12) +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())) 

# Create list of raster brick files 
Diap_wtd_fls <- list.files(path = "DDRP_results", pattern = "all_weighted.tif", full.names = TRUE, recursive = TRUE)
Diap_wtd_fls <- Diap_wtd_fls[grep(pattern = "APH2_EURO", x = Diap_wtd_fls, fixed = TRUE)]


# FullGen plots are the same regardless of photoperiod response
maptype <- nam <-  "FullGen"
index <- 1
fls_10yr <- Diap_wtd_fls[grep(pattern = "FullGen", x = Diap_wtd_fls, fixed = TRUE)]
fls_10yr <- fls_10yr[grep(pattern = "Kyushu", x = fls_10yr, fixed = TRUE)]

brks <- lapply(X = fls_10yr, FUN = function(x) brick(x)[[6]]) # for bimonthly maps, take last result
rankbrk <- lapply(brks, FUN = function(x) mean(getValues(x), na.rm = TRUE))
# 
# # rank by degree-days? same years (2017 coldest, 2018 warmest for Europe)
# gdd_fls <- list.files(path = "DDRP_results", pattern = "DDtotal_cohort1.tif", full.names = TRUE, recursive = TRUE)
# gdd_fls <- gdd_fls[grep(pattern = "APH2_EURO", x = gdd_fls, fixed = TRUE)]
# brks <- lapply(X = gdd_fls[seq(1, 19, by = 2)], FUN = function(x) brick(x)[[6]]) # for bimonthly maps, take last result
# rankbrk <- lapply(brks, FUN = function(x) mean(getValues(x), na.rm = TRUE))


for (quant in c("mean", "2017", "2018")){
  # for (quant in c("2010", "2018")){
  if (quant == "2017"){
    brk <- stack(brks[[8]])
    yr <- paste("2017 (coolest)")
  }
  if (quant == "2018"){
    brk <- stack(brks[[9]])
    yr <- paste("2018 (warmest)")
  }
  # if (quant == "min"){
  #   brk <- min(stack(brks), na.rm = TRUE)
  #   yr <- paste("2010-2019 (minimum)")
  # }
  if (quant == "mean"){
    brk <- mean(stack(brks), na.rm = TRUE)
    yr <- paste("2010-2019 (average)")
  }
  # if (quant == "max"){
  #   brk <- max(stack(brks), na.rm = TRUE)
  #   yr <- paste("2010-2019 (maximum)")
  # }
  
  df <- ConvDF(brk)
  # nam <- "FullGen"
  titl <- paste(case_when(nam == "AttVolt" ~ "Attempted voltinism",
                          nam == "FullGen" ~ "Potential voltinism",
                          nam == "Diapause" ~ "Chose to diapause",
                          nam == "Mismatch" ~ "Voltinism mismatch"), "for", yr, sep = " ")
  subtitl <- "No photoperiod response"
  base_size <- 14
  
  reg.df <- map_data("world")
  # sites$shp <- "Other"
  
  if(nam %in% c("AttVolt", "FullGen")){
    # too many generations is overwhelming on map 
    # TODO: add maxgen+ as highest discrete generation value
    
    plotmax <- max(df$value / 1000)
    df$value[df$value >= 6000] <- 6000
    df$value <- as.factor(round(df$value / 1000 / .5) * .5)
    # df$value <- as.factor(round(df$value / 1000))
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
      # scale_fill_viridis(discrete = TRUE, name = "Generations", 
      # begin = 0, end = plotmax/6) +
      scale_fill_manual(name = "Generations",
                        values = Turbo(out.colors = 13))+ ## This is the change
      # theme_map(base_size = base_size) +
      theme(legend.position = "right") +
      ggtitle(titl, subtitle = subtitl) +
      xlab("Longitude (°)") +
      ylab("Latitude (°)") +
      guides(colour=FALSE)
    
    ggsave(p,file=paste0(nam,"_", quant, ".png"), height = 6, width = 6,
           units = c('in'), dpi = 300) 
  }
  
  if(nam == "Diapause"){
    # TODO Diapause map doesn't match lost generations because of pre-diap?
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
    
    ggsave(p,file=paste0(nam,"_", quant, ".png"), height = 6, width = 6,
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
    
    ggsave(p,file=paste0(nam,"_", quant, ".png"), height = 6, width = 6,
           units = c('in'), dpi = 300) 
  } 
}



fls_all <- expand.grid(maptype = c("AttVolt", "Diapause", "Mismatch"),
                       site = c("Kyushu", "Hokkaido"))

# redo Attvolt
fls_all <- expand.grid(maptype = "AttVolt",
                       site = c("Kyushu", "Hokkaido"))

for (index in 1:nrow(fls_all)){
  
  fls_10yr <- Diap_wtd_fls[grep(pattern = fls_all$maptype[index], x = Diap_wtd_fls, fixed = TRUE)]
  fls_10yr <- fls_10yr[grep(pattern = fls_all$site[index], x = fls_10yr, fixed = TRUE)]
  sitename <- fls_all$site[index]
  brks <- lapply(X = fls_10yr, FUN = function(x) brick(x)[[6]]) # for bimonthly maps, take last result
  
  for (quant in c("mean", "2017", "2018")){
    if (quant == "2017"){
      brk <- stack(brks[[8]])
      yr <- paste("2017 (coolest)")
    }
    if (quant == "2018"){
      brk <- stack(brks[[9]])
      yr <- paste("2018 (warmest)")
    }
    if (quant == "mean"){
      brk <- mean(stack(brks), na.rm = TRUE)
      yr <- paste("2010-2019 (average)")
    }
    
    df <- ConvDF(brk)
    nam <- fls_all$maptype[index]
    titl <- paste(case_when(nam == "AttVolt" ~ "Attempted voltinism",
                            nam == "FullGen" ~ "Potential voltinism",
                            nam == "Diapause" ~ "Chose to diapause",
                            nam == "Mismatch" ~ "Voltinism mismatch"), "for", yr, sep = " ")
    subtitl <- paste0(as.character(sitename), " photoperiod response")
    base_size <- 14
    
    reg.df <- map_data("world")
    # sites$shp <- "Other"
    
    if(nam %in% c("AttVolt", "FullGen")){
      # too many generations is overwhelming on map 
      # TODO: add maxgen+ as highest discrete generation value
      
      plotmax <- max(df$value / 1000)
      df$value[df$value >= 6000] <- 6000
      df$value <- as.factor(round(df$value / 1000 / .5) * .5)
      # df$value <- as.factor(round(df$value / 1000))
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
        # scale_fill_viridis(discrete = TRUE, name = "Generations", 
        # begin = 0, end = plotmax/6) +
        scale_fill_manual(name = "Generations",
                          values = Turbo(out.colors = 13))+ ## This is the change
        # theme_map(base_size = base_size) +
        theme(legend.position = "right") +
        ggtitle(titl, subtitle = subtitl) +
        xlab("Longitude (°)") +
        ylab("Latitude (°)") +
        guides(colour=FALSE)
      
      ggsave(p,file=paste0(sitename, "_", nam,"_", quant, ".png"), height = 6, width = 6,
             units = c('in'), dpi = 300) 
    }
    
    if(nam == "Diapause"){
      # TODO Diapause map doesn't match lost generations because of pre-diap?
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
      
      ggsave(p,file=paste0(sitename, "_", nam,"_", quant, ".png"), height = 6, width = 6,
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
      
      ggsave(p,file=paste0(sitename, "_", nam,"_", quant, ".png"), height = 6, width = 6,
             units = c('in'), dpi = 300) 
    } 
  }
}
