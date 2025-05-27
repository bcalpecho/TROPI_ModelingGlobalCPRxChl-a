# ---
# title: Extracting OC-CCI data
# author: Bryan Alpecho
# date: 2025-
# output: RData for extracted chl-a
# ---

#Aims
#1. aggregate raster file of OC-CCI (based on 'wrangle-netcdf-imos' script from 2024 UQ MME Lab R Workshop)
#2. extract raster data at CPR sampling points
#3. export chl-a data

setwd("C:/Users/power/OneDrive - Université Libre de Bruxelles/4 UQ/CurrentWork/db")
# Load libraries

  library(tidyverse)
  library(stars)

## Load data
# read in raster (chl-a)
  
  #file.list <- list.files(path = "C:/Users/power/Desktop/db/OC-CCI/AusCPR - unedited", pattern = "*\\.nc4", full.names = TRUE)
  file.list <- list.files(path = "C:/Users/power/Desktop/OC-CCI/Global", pattern = "*\\.nc4", full.names = TRUE)
  
  #set study area for extraction 
    #auscpr
    study_area <- read_csv("Input/CPR on-process/cpr_aus_complete.csv", col_names=TRUE) %>% 
      select(TripCode, Sample_ID, Region, Latitude, Longitude, SampleTime_UTC) %>% 
      st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)
    
    #npacific
    study_area <- read_csv("Input/CPR on-process/cpr_npacific_complete.csv", col_names=TRUE) %>% 
      select(Latitude, Longitude, Date, Time) %>% 
      st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)
    
    #natlantic
    study_area <- read_csv("Input/CPR on-process/cpr_natlantic_coords.csv", col_names=TRUE) %>% 
      select(Latitude, Longitude, Year, Month, Day) %>% 
      st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)
    
    #socpr
    study_area <- read_csv("Input/CPR on-process/cpr_SO_complete.csv", col_names=TRUE) %>% 
      select(Tow_Number, Latitude, Longitude, Date, Time) %>% 
      st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)
    
#extract from netcdf function
  extract_ncdf <- function(ncdf_list) {
    filenames <- basename(ncdf_list)
    time.l <- list()
    for(j in 1:length(ncdf_list)){
      
      dat <- read_ncdf(ncdf_list[j], var = c('chlor_a'), make_time = T, proxy = T)
       
      #extract time dimension from ncdf
      time_values <- data.frame(st_get_dimension_values(dat, "time"))
      time.l <- c(time.l, list(time_values))  
      #print(raster_files[i])
      #dat <- raster(raster_files[i])
    
      # make a grid to chunk for processing 
      chunk <- st_make_grid(st_bbox(dat), cellsize = c(10, 10)) %>% st_as_sf()
      
      # loop through chunks and make coarser res
      tmp <- list()
      for(i in 1:nrow(chunk)){ #is not atomic error
        tmp[[i]] <- dat[st_bbox(chunk[i,])] %>% 
          st_as_stars(proxy = F) %>% 
          st_warp(cellsize = c(0.2, 0.2), crs = st_crs(4326))
      }
      # bind processed raster cubes back together
      dat_all <- do.call(st_mosaic, tmp)
        #turns off use of spherical geometry
        sf_use_s2(FALSE)
        #assumes planar coordinates
      
      dat_all2 <- dat_all[st_bbox(study_area)] %>%
        st_set_dimensions("band", values = time_values[,1], names = "time")
      # sequence generation of bands from 1 to 17 by 1 step
      # dat_all2 #shows summary of stars object bydimensions and attribute
     
      #save
      print(paste("file #:",j))
      print(paste(sub("*\\.nc4","",filenames[j]),".grd",sep=""))
      write_stars(dat_all2, paste("Output/raster/",sub("*\\.nc4","",filenames[j]),".grd",sep=""))
      #write_stars(dat_all2, 'Output/raster/')
      
    }
    saveRDS(time.l, file="Output/raster/raster_time_list.RData")
  }
  
  extract_ncdf(file.list)
  
# read in cpr spatiotemporal coordinates 
    #auscpr
    cpr_coords <- read_csv("Input/CPR on-process/cpr_aus_complete.csv", col_names=T) %>% 
      select(TripCode, Sample_ID, Region, Latitude, Longitude, SampleTime_UTC) %>% 
      st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)
      
    #npacific
    cpr_coords <- read_csv("Input/CPR on-process/cpr_npacific_complete.csv", col_names=T) %>% 
      select(Latitude, Longitude, Date, Time) %>% 
      st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)
    
    #natlantic
    cpr_coords <- read_rds( "Output/data/natlantic/cpr_natlantic_complete.rds") %>% 
      select(Sample_ID, SampleTime_UTC, Longitude, Latitude) %>% 
      st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)
    
    #socpr
    cpr_coords <- read_csv("Input/CPR on-process/cpr_SO_complete.csv", col_names=T) %>% 
      select(Tow_Number, Latitude, Longitude, Date, Time) %>% 
      st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)
    
    #global
    auscpr <- read_csv("Input/CPR on-process/cpr_aus_complete.csv", col_names=T) %>% 
      select(TripCode, Sample_ID, Region, Latitude, Longitude, SampleTime_UTC) %>% 
      st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)
    
    npacific <- read_csv("Input/CPR on-process/cpr_npacific_complete.csv", col_names=T) %>% 
      select(Latitude, Longitude, Date, Time) %>% 
      st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

    natlantic <- read_rds( "Output/data/natlantic/cpr_natlantic_complete.rds") %>% 
      select(Sample_ID, SampleTime_UTC, Longitude, Latitude) %>% 
      st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

    socpr <- read_rds()
    
    
####Extract chl-a at sampling point and time####
  raster_files <- list.files(path = "Output/raster/auscpr/global8day/", pattern = "*\\.grd$", full.names = TRUE)  
  
  #read time list for raster files
  time.l <- readRDS("Output/raster/auscpr/global8day/raster_time_list.RData")
  head(time.l)
    #
    extract_chla <- function(rast_list){
      chla_df <- cpr_coords
      for(i in 1:length(rast_list)){
        
        #insert time coordinates
        chla <- read_stars(paste("Output/raster/auscpr/global8day/",basename(rast_list[i]),sep="")) %>% 
          st_set_dimensions("band", values = time.l[[i]][,1], names = "time")
        
        #extract
        #attr(aus_chla, "dimensions")$time$point = FALSE
        file.name <- paste("Output/raster/auscpr/",basename(rast_list[i]),sep="")
        print(paste("file #:",i))
        print(file.name)
        
        chla_ext <- st_extract(chla, cpr_coords, time_column = "SampleTime_UTC", interpolate_time = F)
        
        chla_df <- c(chla_df, chla_ext[1])
      
      }
      chla_df <- cbind(chla_df, cpr_coords %>% select("SampleTime_UTC"))
      saveRDS(chla_df, file="Output/data/auscpr/chla_extracted_noInterpolate.rds")
    }
    
    extract_chla(raster_files)
    
    #merge into a single column for non-empty columns
    chla <- read_rds("Output/data/auscpr/chla_extracted_NOInterpolate.rds")
    
    #exporting
    #chla - true interpolation
    #chla2 - false interpolation
    chla <- chla[1:59]
    chla$geometry = NULL
    
    chla_long <- chla %>%
      select(-c(TripCode, Region, SampleTime_UTC)) %>%
      pivot_longer(-Sample_ID) %>% 
      filter(value != "")
    
    saveRDS(chla_long2, file="Output/data/chla_extracted_long.rds")
      