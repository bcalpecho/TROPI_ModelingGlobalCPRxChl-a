# ---
# title: generate_CompleteDF
# author: Bryan Alpecho
# date: 2025-
# output: csv file
# ---

#Aims
##1. Compute for proportions of feeding groups
##2. Generate values of random effects and combine variables altogether into a dataframe
##3. Update df with proportions of feeding groups
##4. generate complete DF (global) inc. compute for ratios

# Load libraries
library(tidyverse)
library(sf)
library(stars)

##1. Compute for proportions of feeding groups
    #00 get abundance list 
    file.list <- list.files(path = "data_input/CPR_on-process/", pattern = "*\\complete_10062025.csv", full.names = TRUE)
    
    compute_sums <- function(abundance_list){
      #01 get trait list
      traits <- read_csv("Output/traits/TraitTable_06-06-2025.csv", col_names=T) %>% 
        rename(FG_complete = FG_latest)
      
      #aphiaIDs
      aphia_FF <- traits %>% 
        filter(FG_complete == "filter-feeder")
      
      aphia_omni <- traits %>% 
        filter(FG_complete == "omnivore")
      
      aphia_carni <- traits %>% 
        filter(FG_complete == "carnivore")
      
      aphia_na <- traits %>% 
        filter(FG_complete == "na")
      
      aphia_copepods <- traits %>% 
        filter(class == "Copepoda") %>% 
        filter(!is.na(aphiaID))
      
      aphia_cope_omni <- aphia_copepods %>% 
        filter(FG_complete == "omnivore")
      
      aphia_cope_carni <- aphia_copepods %>% 
        filter(FG_complete == "carnivore")
      
      aphia_copepods_noHarp <- traits %>% 
        filter(class =="Copepoda") %>% 
        filter(order != "Harpacticoida" & order != "Siphonostomatoida")
      
      aphia_calanoids <- traits %>% 
        filter(order == "Calanoida")
      
      aphia_calanoid_omni <- aphia_calanoids %>% 
        filter(FG_complete == "omnivore")
      
      aphia_calanoid_carni <- aphia_calanoids %>% 
        filter(FG_complete == "carnivore")
      
      aphia_cyclopoids <- traits %>% 
        filter(order == "Cyclopoida")
      
      aphia_cyclopoid_omni <- aphia_cyclopoids %>% 
        filter(FG_complete == "omnivore")
      
      aphia_cyclopoid_carni <- aphia_cyclopoids %>% 
        filter(FG_complete == "carnivore")
      
      aphia_siphonostomatoida <- traits %>% 
        filter(order == "Siphonostomatoida")
      
      aphia_harpaticoida <- traits %>% 
        filter(order == "Harpacticoida" | order =="Siphonostomatoida")
      
      aphia_harpacticoid_omni <- aphia_harpaticoida %>% 
        filter(FG_complete == "omnivore")
      
      aphia_harpacticoid_carni <- aphia_harpaticoida %>% 
        filter(FG_complete == "carnivore")
      
      #
      for(i in 1:length(abundance_list)){
        #get metadata
        filenames <- basename(abundance_list)
        survey <- str_extract(filenames, "(?<=_)[^_]+") 
        date <- str_extract(filenames, "[^_]+$")
        
        #loop progress
        print(survey[i])
        print(date[i])
        
        #read in cpr abundance data
        cpr <- read_csv(abundance_list[i], col_names = T, name_repair = "minimal")
        cpr_dat <- cpr %>%  select(-c(Sample_ID))
        
        #compute for sums
        cpr_dat_noNAaphia <- cpr_dat %>% 
          select(-matches(as.character(aphia_na$aphiaID)))
        
        cpr_dat_copepods <- cpr_dat %>% 
          select(matches(as.character(aphia_copepods_noHarp$aphiaID)))
        # compute for total abundance of omnivores
        Cope_omni <- cpr_dat_copepods %>%  select(.,matches(as.character(aphia_cope_omni$aphiaID)))
        Cope_omni_sum <- rowSums(Cope_omni)
        # compute for total abundance of carnivores
        Cope_carni <- cpr_dat_copepods %>% select(.,matches(as.character(aphia_cope_carni$aphiaID)))
        Cope_carni_sum <- rowSums(Cope_carni)
        
        #978 natlantic taxa with known feeding group
        #118 natlantic taxa without feeding group (Phytoplankton, Too general)
        #noFG <- natlantic_dat %>% select(.,matches(as.character(aphia_na$aphiaID))) 
        
        # compute for total abundance of filter-feeder
        FF <- cpr_dat_noNAaphia %>% select(.,matches(as.character(aphia_FF$aphiaID)))
        FF_sum <- rowSums(FF)
        # compute for total abundance of omnivores
        Omni <- cpr_dat_noNAaphia %>%  select(.,matches(as.character(aphia_omni$aphiaID)))
        Omni_sum <- rowSums(Omni)
        # compute for total abundance of carnivores
        Carni <- cpr_dat_noNAaphia %>% select(.,matches(as.character(aphia_carni$aphiaID)))
        Carni_sum <- rowSums(Carni)
        # compute for total abundance of non filter-feeder
        NonFF <- cpr_dat_noNAaphia %>% select(.,-matches(as.character(aphia_FF$aphiaID)))
        NonFF_sum <- rowSums(NonFF)
        # compute for total abundance of non carnivores
        NonCarni <- cpr_dat_copepods %>% select(.,-matches(as.character(aphia_carni$aphiaID)))
        NonCarni_sum <- rowSums(NonCarni)
        # compute for total abundance of non omnivores
        NonOmni <- cpr_dat_copepods %>% select(.,-matches(as.character(aphia_omni$aphiaID)))
        NonOmni_sum <- rowSums(NonOmni)
        
        # compute for total abundance of calanoids
        #omnivore calanoids
        Cal_omni <- cpr_dat_copepods %>%  select(.,matches(as.character(aphia_calanoid_omni$aphiaID)))
        Cal_omni_sum <- rowSums(Cal_omni)
        #carnivore calanoids
        Cal_carni <- cpr_dat_copepods %>%  select(.,matches(as.character(aphia_calanoid_carni$aphiaID)))
        Cal_carni_sum <- rowSums(Cal_carni)
        # compute for total abundance of cyclopoids
        #omnivore cyclopoids
        Cyc_omni <- cpr_dat_copepods %>%  select(.,matches(as.character(aphia_cyclopoid_omni$aphiaID)))
        Cyc_omni_sum <- rowSums(Cyc_omni)
        #carnivore cyclopoids
        Cyc_carni <- cpr_dat_copepods %>%  select(.,matches(as.character(aphia_cyclopoid_carni$aphiaID)))
        Cyc_carni_sum <- rowSums(Cyc_carni)
        # compute for total abundance of harpaticoids
        #omnivore harpac
        Har_omni <- cpr_dat_copepods %>%  select(.,matches(as.character(aphia_harpacticoid_omni$aphiaID)))
        Har_omni_sum <- rowSums(Har_omni)
        
        # bind all sum
        trait_sums <- data.frame(rbind(cpr$Sample_ID,FF_sum,Omni_sum, Cope_omni_sum, Carni_sum, Cope_carni_sum, Cal_omni_sum, Cal_carni_sum, Cyc_omni_sum, Cyc_carni_sum, Har_omni_sum, NonFF_sum, NonCarni_sum, NonOmni_sum)) %>% 
          t() %>% 
          data.frame() %>% 
          rename(Sample_ID = V1)
        
        print(paste("Output/data/trial/df_sums_",survey[i],"_",date[i], sep=""))
        write_csv(trait_sums, paste("Output/data/trial/df_sums_",survey[i],"_",date[i], sep=""))
      }
    }
    compute_sums(file.list)
    

##2. Combining variables altogether into a dataframe (shown for NAtlantic)
    #set date for version control 
    date <- "12062025"
    
    #load cpr and keep metadata
    df <- read_rds("Output/data/natlantic/cpr_natlantic_complete.rds") %>% 
      select(Sample_ID, Latitude, Longitude, SampleTime_UTC)
    
    #01 assign Tow_No
    df <- df %>% 
      mutate(Tow_No = str_extract(Sample_ID, "^.*?(?=-)")) 
    
    #02 compute Tow_Days
    df <- df %>% 
      mutate(SampleTime_UTC = as.POSIXct(paste(Day, Month, Year, Hour, sep="/"), format = "%d/%m/%Y/%H", tz = "UTC")) %>% 
      group_by(Tow_No) %>% 
      mutate(Tow_Days = as.double(difftime(SampleTime_UTC, min(SampleTime_UTC), units = "days"))) %>% 
      ungroup()
    
    #03 intersect Longhurst Provinces
    #load longhurst provinces
    LH <- st_read(dsn = "Input/LonghurstProvinces/", layer = "Longhurst_world_v4_2010")
    #downloaded from https://marineregions.org/downloads.php
    
    st_is_valid(LH)
    LH <- st_make_valid(LH)
    
    df <- df %>% 
      st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>% 
      st_transform(crs = st_crs(LH)) 
    
    coords <- read_rds("Output/data/natlantic/cpr_natlantic_complete.rds") %>% 
      select(Sample_ID, Latitude, Longitude)
    df <- st_intersection(df, LH) %>% 
      select(-c(ProvDescr, geometry)) %>% 
      left_join(coords, by="Sample_ID") %>% 
      relocate(c("Longitude","Latitude"), .after=(Sample_ID)) %>% 
      rename("Longhurst" = "ProvCode")
    
    #load chlorophyll-a
    ## Extracted chlorophyll matched by time and space
    chl <- read_rds("Output/data/natlantic/chla_extracted_combined_26052025.rds")
    
    df <- df %>% 
      left_join(chl, by = "Sample_ID") %>% 
      rename(Chla = chla) 
    
    #load proportions
    ## Relative abundances of filter-feeders, omnivores, and carnivores 
    #read df1 below first
    traits <- read_csv("Output/data/natlantic/df_sums_05062025.csv", col_names =T) 
    traits <- traits %>% 
      rename(Sample_ID = V1)
    df <- df %>% 
      left_join(traits, by = "Sample_ID")
    df$geometry = NULL
    
    write_rds(df, paste("Output/data/natlantic/df_complete_",date,".rds",sep=""))


##3. Update df with proportions of feeding groups
    #set date for version control
    date <- "12062025"
    
    #
    survey <- c("auscpr","natlantic","npacific","socpr")
    file.list <- NA
    for(i in 1:length(survey)){
      file.list[i] <- list.files(path = paste("Output/data/",survey[i],"/",sep=""), pattern = "*\\.rds", full.names = TRUE)
    }
    remove(i)
    
    update_df_traits <- function(df_list){
      for(i in 1:length(df_list)){
        #get metadata
        filenames <- basename(file.list)
        survey <- str_extract(df_list, "(?<=/data/)[^/]+") 
        date <- str_extract(filenames, "[^_]+$")
        
        #progress
        print(survey[i])
        
        #read in previous df
        df_prev <- paste("Output/data/",survey[i],"/",filenames[i],sep="")
        df <- read_rds(df_prev)
        
        #remove previous traits
        df <- df %>% select(-c(FF_sum, Omni_sum, Cope_omni_sum, Carni_sum, Cope_carni_sum, Cal_omni_sum, Cal_carni_sum, Cyc_omni_sum, Cyc_carni_sum, Har_omni_sum, NonFF_sum, NonOmni_sum, NonCarni_sum))
        
        #load updated traits
        traits_file <- paste("Output/data/",survey,"/df_sums_09062025.csv", sep="")
        traits <- read_csv(traits_file)
        
        #update df
        df_new <- df %>% 
          left_join(traits, by = "Sample_ID")
        
        #save updated df
        write_rds(df_new, paste("Output/data/",survey[i],"/df_complete_",date,".rds",sep=""))
      } 
    }

update_df_traits(file.list)

##4. generate complete DF (global) inc. compute for ratios
    #set date for version control
    current_version <- "11062025"
    date <- "12062025"

    #Combining the survey-specific CPR datasets
    survey <- c("auscpr","natlantic","npacific","socpr")
      for(i in 1:length(survey)){
        assign(paste(survey[i]), read_rds(paste("Output/data/",survey[i],"/df_complete_",current_version,".rds",sep="")) %>%  mutate(Survey = survey[i]))
      }
    
    df <- auscpr %>% 
      rows_insert(socpr, by = "Sample_ID") %>% 
      rows_insert(npacific, by ="Sample_ID") %>% 
      rows_insert(natlantic, by = "Sample_ID")
    
    #transform chl-a and compute for feeding group ratios
    df <- df %>% 
      mutate(Chla_sqrt = sqrt(Chla)) %>%
      mutate(ROC = Cope_omni_sum/(Cope_omni_sum + Cope_carni_sum)) %>% 
      mutate(RCO = Cope_carni_sum/(Cope_omni_sum + Cope_carni_sum)) %>% 
      mutate(RFF = FF_sum/(FF_sum+NonFF_sum)) %>% 
      mutate(ROC_Cal = Cal_omni_sum/(Cal_omni_sum + Cal_carni_sum)) %>% 
      mutate(RCO_Cal = Cal_carni_sum/(Cal_omni_sum + Cal_carni_sum)) %>% 
      mutate(ROC_Cyc = Cyc_omni_sum/(Cyc_omni_sum + Cyc_carni_sum)) %>% 
      mutate(RCO_Cyc = Cyc_carni_sum/(Cyc_omni_sum + Cyc_carni_sum)) 
    
    df <- df %>% 
      mutate(RFF_SVT = case_when(
        RFF == 0 | RFF == 1 ~ ((RFF*325075)+0.5)/325076,
        .default = RFF)) %>% 
      mutate(ROC_SVT = case_when(
        ROC == 0 | ROC == 1 ~ ((ROC*325075)+0.5)/325076,
        .default = ROC)) %>% 
      mutate(RCO_SVT = case_when(
        RCO == 0 | RCO == 1 ~ ((RCO*325075)+0.5)/325076,
        .default = RCO)) %>%  
      mutate(ROC_Cal_SVT = case_when(
        ROC_Cal == 0 | ROC_Cal == 1 ~ ((ROC_Cal*325075)+0.5)/325076,
        .default = ROC_Cal)) %>% 
      mutate(RCO_Cal_SVT = case_when(
        RCO_Cal == 0 | RCO_Cal == 1 ~ ((RCO_Cal*325075)+0.5)/325076,
        .default = RCO_Cal)) %>%
      mutate(ROC_Cyc_SVT = case_when(
        ROC_Cyc == 0 | ROC_Cyc == 1 ~ ((ROC_Cyc*325075)+0.5)/325076,
        .default = ROC_Cyc)) %>% 
      mutate(RCO_Cyc_SVT = case_when(
        RCO_Cyc == 0 | RCO_Cyc == 1 ~ ((RCO_Cyc*325075)+0.5)/325076,
        .default = RCO_Cyc)) %>%  
      mutate(RFF_SVT_zib = case_when(
        RFF == 1 ~ ((RFF*325075)+0.5)/325076,
        .default = RFF)) %>% 
      mutate(ROC_SVT_zib = case_when(
        ROC == 1 ~ ((ROC*325075)+0.5)/325076,
        .default = ROC)) %>% 
      mutate(RCO_SVT_zib = case_when(
        RCO == 1 ~ ((RCO*325075)+0.5)/325076,
        .default = RCO)) %>% 
      mutate(ROC_Cal_SVT_zib = case_when(
        ROC_Cal == 1 ~ ((ROC_Cal*325075)+0.5)/325076,
        .default = ROC_Cal)) %>% 
      mutate(RCO_Cal_SVT_zib = case_when(
        RCO_Cal == 1 ~ ((RCO_Cal*325075)+0.5)/325076,
        .default = RCO_Cal)) %>% 
      mutate(ROC_Cyc_SVT_zib = case_when(
        ROC_Cyc == 1 ~ ((ROC_Cyc*325075)+0.5)/325076,
        .default = ROC_Cyc)) %>% 
      mutate(RCO_Cyc_SVT_zib = case_when(
        RCO_Cyc == 1 ~ ((RCO_Cyc*325075)+0.5)/325076,
        .default = RCO_Cyc)) 
    
    #finalize df without Group
    df <- df %>% 
      relocate(Survey, .before= Sample_ID) %>% 
      relocate(Chla_sqrt, .after=Chla)
    
    write_rds(df, paste("Output/data/global/df_complete_",date,".rds",sep=""))   
    
    #subdivide data into cyclopoid and calanoid
    
    df_cal <- df %>% 
      select(c(Survey,Sample_ID,Longitude, Latitude, SampleTime_UTC, Tow_No, Tow_Days, Longhurst, Chla, Chla_sqrt, FF_sum, Omni_sum, Cope_omni_sum, Carni_sum, Cope_carni_sum, ROC_Cal,ROC_Cal_SVT,ROC_Cal_SVT_zib,RCO_Cal,RCO_Cal_SVT,RCO_Cal_SVT_zib)) %>% 
      mutate(Group = "Calanoids") %>% 
      rename(ROC = ROC_Cal) %>% 
      rename(ROC_SVT = ROC_Cal_SVT) %>% 
      rename(ROC_SVT_zib = ROC_Cal_SVT_zib) %>% 
      rename(RCO = RCO_Cal) %>% 
      rename(RCO_SVT = RCO_Cal_SVT) %>% 
      rename(RCO_SVT_zib = RCO_Cal_SVT_zib)
    
    df_cyc <- df %>% 
      select(c(Survey,Sample_ID,Longitude, Latitude, SampleTime_UTC, Tow_No, Tow_Days, Longhurst, Chla, Chla_sqrt, FF_sum, Omni_sum, Cope_omni_sum, Carni_sum, Cope_carni_sum, ROC_Cyc,ROC_Cyc_SVT,ROC_Cyc_SVT_zib,RCO_Cyc,RCO_Cyc_SVT,RCO_Cyc_SVT_zib)) %>% 
      mutate(Group = "Cyclopoids") %>% 
      rename(ROC = ROC_Cyc) %>% 
      rename(ROC_SVT = ROC_Cyc_SVT) %>% 
      rename(ROC_SVT_zib = ROC_Cyc_SVT_zib) %>% 
      rename(RCO = RCO_Cyc) %>% 
      rename(RCO_SVT = RCO_Cyc_SVT) %>% 
      rename(RCO_SVT_zib = RCO_Cyc_SVT_zib)
    
    df_withGroup <- bind_rows(df_cal, df_cyc)    
    
    write_rds(df_withGroup, paste("Output/data/global/df_complete_withGroup_",date,".rds",sep=""))

#12.06.2025 latest revision