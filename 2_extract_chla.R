# ---
# title: Extracting OC-CCI data
# author: Bryan Alpecho
# date: 2025-
# output: RData for extracted chl-a
# ---

#Aims
#1. aggregate raster file of 8-day OC-CCI (based on 'wrangle-netcdf-imos' script from 2024 UQ MME Lab R Workshop)
#2. extract 8-day OC-CCI (raster data) at CPR sampling points
#3. fill-up gaps of 8-day by monthly OC-CCI values

# Load libraries
  library(tidyverse)
  library(stars)
  library(raster)

#1. aggregate raster file of OC-CCI (based on 'wrangle-netcdf-imos' script from 2024 UQ MME Lab R Workshop)
  #set date for version control
  date <- "02082025"
  
  # repeat the process per survey
  ## Load data
  # read in list of netcdf files
    file.list <- list.files(path = "data_input/OC-CCI/Global", pattern = "*\\.nc4", full.names = TRUE)
    file.list <- list.files(path = "C:/Users/power/Desktop/OC-CCI/Global", pattern = "*\\.nc4", full.names = TRUE) #temporary file source due to large file size
    file.list <- list.files(path = "C:/Users/power/Desktop/OC-CCI/Global monthly", pattern = "*\\.nc4", full.names = TRUE) #temporary file source due to large file size

  #set study area for extraction 
    #auscpr
    study_area <- read_csv("data_input/CPR_on-process/cpr_auscpr_metadata.csv", col_names=TRUE) %>% 
      dplyr::select(latitude, longitude) %>% 
      st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
    
    #npacific
    study_area <- read_csv("data_input/CPR_on-process/cpr_npacific_metadata.csv", col_names=TRUE) %>% 
      dplyr::select(latitude, longitude) %>% 
      st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
    
    #natlantic
    study_area <- read_csv("data_input/CPR_on-process/cpr_natlantic_metadata.csv", col_names=TRUE) %>% 
      dplyr::select(latitude, longitude) %>% 
      st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
    
    #socpr
    study_area <- read_csv("data_input/CPR_on-process/cpr_socpr_metadata.csv", col_names=TRUE) %>% 
      dplyr::select(latitude, longitude) %>% 
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
        print(paste(sub("*\\.nc4","",filenames[j]),"_",date,".grd",sep=""))
        write_stars(dat_all2, paste("Output/raster/monthly/",sub("*\\.nc4","",filenames[j]),"_",date,".grd",sep=""))
        #write_stars(dat_all2, 'Output/raster/')
        
      }
      saveRDS(time.l, file=paste("Output/raster/monthly/raster_time_list_",date,".RData",sep=""))
    }
  extract_ncdf(file.list)

#2. extract raster data at CPR sampling points
# repeat the process per survey
# set date for version control
  date <- "03082025"
  
   #read in cpr spatiotemporal coordinates 
    #auscpr
    cpr_coords <- read_csv("data_input/CPR_on-process/cpr_auscpr_metadata.csv", col_names=T) %>% 
      select(sample_id, latitude, longitude, sampleTime_utc) %>% 
      st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
    survey <- "auscpr"
      
    #npacific
    cpr_coords <- read_csv("data_input/CPR_on-process/cpr_npacific_metadata.csv", col_names=T) %>% 
      st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
    survey <- "npacific"
    
    #natlantic
    cpr_coords <- read_csv( "data_input/CPR_on-process/cpr_natlantic_metadata.csv") %>% 
      select(sample_id, latitude, longitude, sampleTime_utc) %>% 
      st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
    survey <- "natlantic"
    
    #socpr
    cpr_coords <- read_csv("data_input/CPR_on-process/cpr_socpr_metadata.csv", col_names=T) %>% 
      select(sample_id, latitude, longitude, sampleTime_utc) %>% 
      st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
    survey <- "socpr"

  #Extract chl-a at sampling point and time####
    #set temporal resolution
    #temp_res <- "monthly"
    raster_files <- list.files(path = paste("Output/raster/", survey, sep=""), pattern = "*\\.grd$", full.names = TRUE)  
    
    #read time list for raster files
    time.l <- read_rds(paste("Output/raster/",survey,"/raster_time_list_02082025.RData",sep=""))

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
        
        chla_ext <- st_extract(chla, cpr_coords, time_column = "sampleTime_utc", interpolate_time = T)
        
        chla_df <- c(chla_df, chla_ext[1])
      
      }
      chla_df <- cbind(chla_df, cpr_coords %>% select("sample_id", "sampleTime_utc"))
      saveRDS(chla_df, file=paste("Output/data/",survey,"/chla_extracted_",date,".rds",sep=""))
      print(paste("Extracted Chl-a Output: ","Output/data/",survey,"/chla_extracted_",date,".rds",sep=""))
    }
  
  extract_chla(raster_files)

#3. Fill-up gaps by monthly OC-CCI values
  
  #Aims
  #1. To average 8-day chl-a into monthly chl-a
  #2. To fill up gaps of 8-day with monthly chl-a
  
  # set for version control
  date <- "03082025"
  survey <- "npacific"
  
  #read raster files (aggregated to 0.2 x 0.2 deg)
  #raster_files <- list.files(path = paste("Output/raster/", survey, sep=""), pattern = "*\\.grd$", full.names = TRUE)  
  raster_files <- list.files(path = paste("Output/raster/npacific", sep=""), pattern = "*\\.grd$", full.names = TRUE)  
  
  #read time list for raster files
  #time.l <- read_rds(paste("Output/raster/",survey,"/raster_time_list.RData",sep=""))
  time.l <- read_rds(paste("Output/raster/npacific/raster_time_list_02082025.RData",sep=""))
  
  fill_up_gaps <- function(rast_list){
    chla_df <- cpr_coords
    for(i in 1:length(rast_list)){
      #insert time coordinates
        chla <- read_stars(paste("Output/raster/",survey,"/",basename(rast_list[i]),sep="")) %>% 
          st_set_dimensions("band", values = time.l[[i]][,1], names = "time")
     
      #average 8-day into monthly raster
        chla_monthly <- stars::st_apply(chla, "time", function(x) raster::movingFun(x, 3, mean, na.rm=TRUE)) # (3) 8-weeks -> monthly
      
      #extract from monthly raster 
        file.name <- paste("Output/raster/",survey,"/",basename(rast_list[i]),sep="")
        print(paste("file #:",i)) #progress update
        print(file.name)
        
        chla_ext <- stars::st_extract(chla_monthly, cpr_coords, time_column = "SampleTime_UTC", interpolate_time = T)
        
        chla_df <- c(chla_df, chla_ext[1])
        print(paste("Monthly average computed:",file.name,sep=""))  
    }
    #export 
    chla_df <- cbind(chla_df, cpr_coords %>% select("SampleTime_UTC"))
    write_rds(chla_df, file=paste("Output/data/",survey,"/chla_extracted_monthly_",date,".rds",sep=""))
    
      #combine monthly and 8-day chl-a
      chla_monthly <- read_rds("Output/data/chla_extracted_monthly.rds") %>% 
        rename(chla_monthly = value)
      chla_eightday <- read_rds("Output/data/chla_extracted_long_09052025.rds") %>% 
        rename(chla_eightday = value)
      
      chla_df2 <- chla_df %>% 
        select(Sample_ID) %>% 
        left_join(chla_monthly, by="Sample_ID") %>%
        left_join(chla_eightday, by="Sample_ID") %>% 
        select(-c(name.x, name.y))
      
      #fill-up gaps of 8-day by monthly values
      chla_df2 <- chla_df2 %>% 
        mutate(chla_f = ifelse(is.na(chla_eightday), chla_monthly, chla_eightday))
      
      chla_df_final <- chla_df2 %>% 
        select(Sample_ID, SampleTime_UTC, chla_f) %>% 
        rename(chla = chla_f)
    
    file.name <- paste("Output/data/",survey,"/chla_merged_",date,".rds",sep="")
    write_rds(chla_df_final, file=file.name)
    print(paste("Output",file.name,sep=""))  
    
  }

  fill_up_gaps(rast_trial)
  
#27.06.2025 - latest revision: added "fill_up_gaps" function
      