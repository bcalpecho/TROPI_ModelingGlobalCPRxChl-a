# ---
# title: Extracting OC-CCI data
# author: Bryan Alpecho
# date: 2025-
# output: RData for extracted chl-a
# ---

#Aims
#1. aggregate raster file of OC-CCI (based on 'wrangle-netcdf-imos' script from 2024 UQ MME Lab R Workshop)
#2. extract raster data at CPR sampling points

# Load libraries
  library(tidyverse)
  library(stars)

#1. aggregate raster file of OC-CCI (based on 'wrangle-netcdf-imos' script from 2024 UQ MME Lab R Workshop)
  # repeat the process per survey
  ## Load data
  # read in list of netcdf files
    file.list <- list.files(path = "data_input/OC-CCI/Global", pattern = "*\\.nc4", full.names = TRUE)
    
  #set study area for extraction 
    #auscpr
    study_area <- read_csv("data_input/CPR_on-process/cpr_auscpr_metadata.csv", col_names=TRUE) %>% 
      select(Latitude, Longitude) %>% 
      st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)
    
    #npacific
    study_area <- read_csv("data_input/CPR_on-process/cpr_npacific_metadata.csv", col_names=TRUE) %>% 
      select(Latitude, Longitude) %>% 
      st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)
    
    #natlantic
    study_area <- read_csv("data_input/CPR_on-process/cpr_natlantic_metadata.csv", col_names=TRUE) %>% 
      select(Latitude, Longitude) %>% 
      st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)
    
    #socpr
    study_area <- read_csv("data_input/CPR_on-process/cpr_socpr_metadata.csv", col_names=TRUE) %>% 
      select(Latitude, Longitude) %>% 
      st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)
    

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

#2. extract raster data at CPR sampling points
# repeat the process per survey
# set date for version control
  date <- "12062025"
  
# read in cpr spatiotemporal coordinates 
    #auscpr
    cpr_coords <- read_csv("data_input/CPR_on-process/cpr_auscpr_metadata.csv", col_names=T) %>% 
      select(TripCode, Sample_ID, Region, Latitude, Longitude, SampleTime_UTC) %>% 
      st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)
    survey <- "auscpr"
      
    #npacific
    cpr_coords <- read_csv("data_input/CPR_on-process/cpr_npacific_metadata.csv", col_names=T) %>% 
      select(Latitude, Longitude, Date, Time) %>% 
      st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)
    survey <- "npacific"
    
    #natlantic
    cpr_coords <- read_csv( "data_input/CPR_on-process/cpr_natlantic_metadata.csv") %>% 
      mutate(SampleTime_UTC = as.POSIXct(paste(Day, Month, Year, Hour, sep="/"), format = "%d/%m/%Y/%H", tz = "UTC")) %>% 
      select(Latitude, Longitude, SampleTime_UTC) %>% 
      st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)
    survey <- "natlantic"
    
    #socpr
    cpr_coords <- read_csv("data_input/CPR_on-process/cpr_socpr_metadata.csv", col_names=T) %>% 
      select(Tow_Number, Latitude, Longitude, Date, Time) %>% 
      st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)
    survey <- "socpr"

####Extract chl-a at sampling point and time####
  raster_files <- list.files(path = paste("Output/raster/", survey, sep=""), pattern = "*\\.grd$", full.names = TRUE)  
  
  #read time list for raster files
  time.l <- readRDS(paste("Output/raster/",survey,"/raster_time_list.RData",sep=""))

    #
    extract_chla <- function(rast_list){
      chla_df <- cpr_coords
      for(i in 1:length(rast_list)){
        
        #insert time coordinates
        chla <- read_stars(paste("Output/raster/",survey,"/",basename(rast_list[i]),sep="")) %>% 
          st_set_dimensions("band", values = time.l[[i]][,1], names = "time")
        
        #extract
        #attr(aus_chla, "dimensions")$time$point = FALSE
        file.name <- paste("Output/raster/",survey,"/",basename(rast_list[i]),sep="")
        print(paste("file #:",i))
        print(file.name)
        
        chla_ext <- st_extract(chla, cpr_coords, time_column = "SampleTime_UTC", interpolate_time = T)
        
        chla_df <- c(chla_df, chla_ext[1])
      
      }
      chla_df <- cbind(chla_df, cpr_coords %>% select("SampleTime_UTC"))
      saveRDS(chla_df, file=paste("Output/data/",survey,"/chla_extracted_",date,".rds",sep=""))
    }
  
  extract_chla(raster_files)

  #12.06.2025 - revision
      