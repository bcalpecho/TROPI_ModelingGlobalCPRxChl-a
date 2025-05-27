# ---
# title: Completing the dataframe (CPR + chla + trait + Random effects)
# author: Bryan Alpecho
# date: 2025-
# output: rds file for AusCPR and global data
# ---

#Aims
#1. Generate values of random effects
#2. Combine all variables into a dataframe
#3. Export the dataframe as a rds file

setwd("C:/Users/power/OneDrive - Université Libre de Bruxelles/4 UQ/CurrentWork/db")

# Load libraries
library(tidyverse)
library(sf)
library(stars)
# for reference
  campbell <- readRDS("Codes/Campbell et al (2023)/CopeData.RDS")


## Random Effects ##

    #load AusCPR data
      df_large <- read_csv("Input/CPR on-process/NorthAtlantic - latest/cpr_natlantic_largezooplankton.csv", col_names=T) 
      df_small <- read_csv("Input/CPR on-process/NorthAtlantic - latest/cpr_natlantic_smallzooplankton.csv", col_names=T)
      meta_large <- read_csv("Input/CPR on-process/NorthAtlantic - latest/cpr_natlantic_largezooplankton_meta.csv", col_names=T)
      meta_small <- read_csv("Input/CPR on-process/NorthAtlantic - latest/cpr_natlantic_smallzooplankton_meta.csv", col_names=T)
      day <- read_csv("Input/CPR on-process/NorthAtlantic - latest/cpr_natlantic_day.csv", col_names=T)
      
      #combine data 
      df <- df_large %>% 
        left_join(df_small, by=c("SampleId","Latitude","Longitude","Year","Month")) %>% 
        left_join(day, by="SampleId") %>% 
        relocate(Day, .after=Month)
      
      #update 'accepted_id' to 'aphia_id' 
      names(df)[match(meta_large$accepted_id, names(df))] <- meta_large$aphia_id
      names(df)[match(meta_small$accepted_id, names(df))] <- meta_small$aphia_id
        
      
    #01 assign Tow_No
      df <- df %>% 
        mutate(Tow_No = str_extract(SampleId, "([0-9]{3})"))
      
      df <- df %>% 
        group_by(TripCode) %>% 
        mutate(Tow_No = paste(Region, cur_group_id(), sep = "")) #Instead of Region, update it to 'SurveyName' (AUS)
      
    #02 compute Tow_Days
      str(cpr)
      df <- df %>% 
        mutate(time = as.POSIXct(SampleTime_UTC, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")) %>%  
        group_by(TripCode) %>% 
        mutate(Tow_Days = as.double(difftime(time, min(time), units = "days"))) %>% 
        ungroup() %>% 
        select(-c(time))
        
      #export
      write.csv(cpr, "Output/data/cpr_aus_TowNo-TowDays.csv", row.names = FALSE)
      
    #03 intersect Longhurst Provinces
      #load longhurst provinces
      LH <- st_read(dsn = "Input/LonghurstProvinces/", layer = "Longhurst_world_v4_2010") #downloaded from https://marineregions.org/downloads.php
      
      df <- df %>% 
        st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>% 
        st_transform(crs = st_crs(LH)) 
        
      st_is_valid(LH)
      LH <- st_make_valid(LH)
      
      df <- st_intersection(df, LH) %>% 
        select(-c(ProvDescr, TripCode)) %>% 
        relocate(c("geometry"), .after=(Region))
        
      write.csv(inter, "Output/data/Longhurst_intersection.csv", row.names = FALSE) #note there are 5 points left to be determined for its longhurst province

## Compute for proportions of feeding groups ##
    
    #get trait list
    #traits <- read.csv("C:/Users/power/OneDrive - Université Libre de Bruxelles/4 UQ/CurrentWork/db/Input/trait-table/2025.21.04_species_list_with_traits.csv", sep=",", header=T) %>% 
    #  select(aphiaID, taxonID, taxon, traitValue, FG_order, filter_or_not)
    
      traits <- read_csv("Input/trait-table/trait-table-final_26052025.csv", col_names=T) %>% 
        select(aphiaID, taxonID, class, order, scientificName, taxon_cpr_list, FG_mod, filter_or_not) %>% 
        rename(FG_complete = FG_mod)
      
    #%>% 
    #  rename(FG_complete = FG_mod)
    
    #campbell (Sep 2018) updated last 22.05.2025
      campbell <- read_csv("Input/trait-table/Campbelletal_copepods_traits.csv", col_names=T)
      
      auscpr_aphia <- read_csv("Input/CPR on-process/cpr_aus_aphia.csv", col_names=T)
      
      campbell <- campbell %>% 
        left_join(auscpr_aphia, by=c("taxon_name"="taxon"))
      
      campbell <- campbell %>% 
        left_join(traits, by="aphiaID")
      
      write_csv(campbell,"Output/data/traits/campbell-comparison.csv")
      
      campbellVer2 <- read_csv("Output/data/traits/campbell-comparison.csv") %>% 
        select(taxon_no, taxon_name, aphiaID_mod, diet, Size_resolved)
      
      campbellVer2 <- campbellVer2 %>% 
        left_join(traits, by=c("aphiaID_mod"="aphiaID"))  %>% 
        rename(campbell_diet = diet) %>% 
        relocate(campbell_diet,.after=FG_complete)
        #rw 43 and 700 
      
      campbellVer2 <- campbellVer2  %>% 
        mutate(same_or_not = (str_extract(FG_complete,"([a-z]{4})") == campbell_diet)) 
      
      summary(campbellVer2$same_or_not) #22.05.2025 72F 635T
      
      write_csv(campbellVer2, "Output/data/traits/campbell-comparison-ver2.csv")
        #24.05.2025 to know if calanoids are the "unsure traits'
        campbellVer2 <- read_csv("Output/data/traits/campbell-comparison-ver2.csv")
        
        campbell_diet <- campbellVer2 %>% 
          select(-c(taxon_no, Size_resolved, taxon_name)) %>% 
          distinct(aphiaID_mod, .keep_all=T) 
      
        #how to extract portions of data frame based on a logical column
        
      campbell_diet_aphia <- campbell_diet %>%     
          select(c("aphiaID_mod","campbell_diet","same_or_not"))
     
      traits_mod <- traits %>% 
        left_join(campbell_diet_aphia,by=c("aphiaID"="aphiaID_mod"))
      
        #26.05.2025 modified for campbell_diet
        aphia_different <- traits_mod %>% 
          filter(same_or_not == FALSE) %>%
          select(c("aphiaID", "class", "order", "scientificName", "taxon_cpr_list", "FG_complete", "campbell_diet", "same_or_not")) %>% 
          group_by(aphiaID)
        
        print(aphia_different, n= 37)
        summary(as.factor(aphia_different$order))
        
      summary(as.factor(campbellVer2$campbell_diet))
      
      traits_mod <- traits_mod %>% 
        mutate(FG_mod = case_when(
          campbell_diet == "omni" ~ "omnivore",
          campbell_diet == "carn" ~ "carnivore",
          campbell_diet == "herb" ~ "omnivore",
          .default = FG_complete
        )) %>% 
        filter(aphiaID != "1080") #remove damaged copepods
      
      write_csv(traits_mod, "Input/trait-table/trait-table-final_26052025.csv")
      
    #
    
    
    #get abundance list 
      #auscpr <- read_csv("Input/CPR on-process/cpr_aus_complete.csv", col_names=T) %>% 
      #    select(-c(TripCode, Region, Latitude, Longitude, SampleTime_UTC, SampleTime_Local, Year_Local, Month_Local, Day_Local, Time_Local24hr, SatSST_degC, SatChlaSurf_mgm3, PCI, BiomassIndex_mgm3, SampleVolume_m3))
      
      cpr <- read.csv("Input/CPR on-process/cpr_aus_complete_aphia.csv", sep=",", header=T, check.names = F) 
      cpr_dat <- cpr[,2:1097] * 1.5
      cpr_dat_r <- round(cpr_dat) #default decimal place is 0

      
    #    #02 group trait table by feeding group
      aphia_FF <- traits %>% 
        filter(FG_complete == "filter-feeder")
      
      aphia_omni <- traits %>% 
        filter(FG_complete == "omnivore")
      
      aphia_carni <- traits %>% 
        filter(FG_complete == "carnivore")
      
      aphia_na <- traits %>% 
        filter(FG_complete == "na")
      
      aphia_copepods <- traits %>% 
        filter(class == "Copepoda")
      
      aphia_cope_omni <- aphia_copepods %>% 
        filter(FG_complete == "omnivore")
      
      aphia_cope_carni <- aphia_copepods %>% 
        filter(FG_complete == "carnivore")
      
      summary(as.factor(aphia_copepods$order))
      
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
      
      aphia_harpaticoida <- traits %>% 
        filter(order == "Harpacticoida" | order =="Siphonostomatoida")
      
      aphia_harpacticoid_omni <- aphia_harpaticoida %>% 
        filter(FG_complete == "omnivore")
      
      aphia_harpacticoid_carni <- aphia_harpaticoida %>% 
        filter(FG_complete == "carnivore")
      
      #
      
      summary(as.factor(aphia_copepods$FG_complete)) #246 carni, 751 omni copepods; 130 carni, 267 omni 26052025
      
    #modify dat without NA for FeedingGroup #reupdate to rounded values of abundance 
    auscpr_dat <- auscpr_dat_r %>% 
      select(-matches(as.character(aphia_na$aphiaID)))
    
    auscpr_dat_copepods <- auscpr_dat_r %>% 
      select(matches(as.character(aphia_copepods$aphiaID)))
        #978 AusCPR taxa with known feeding group
        #118 AusCPR taxa without feeding group (Phytoplankton, Too general)
        #noFG <- auscpr_dat %>% select(.,matches(as.character(aphia_na$aphiaID))) 
    
    # compute for total abundance
      #filter-feeder
      FF <- auscpr_dat %>% select(.,matches(as.character(aphia_FF$aphiaID)))
      FF_sum <- rowSums(FF)
      # omnivores
      Omni <- auscpr_dat_copepods %>% select(.,matches(as.character(aphia_omni$aphiaID)))
      Omni_sum <- rowSums(Omni)
      # carnivores
      Carni <- auscpr_dat_copepods %>% select(.,matches(as.character(aphia_carni$aphiaID)))
      Carni_sum <- rowSums(Carni)
      # non filter-feeder
      NonFF <- auscpr_dat %>% select(.,-matches(as.character(aphia_FF$aphiaID)))
      NonFF_sum <- rowSums(NonFF)
      # non carnivores
      NonCarni <- auscpr_dat_copepods %>% select(.,-matches(as.character(aphia_carni$aphiaID)))
      NonCarni_sum <- rowSums(NonCarni)
      # non omnivores
      NonOmni <- auscpr_dat_copepods %>% select(.,-matches(as.character(aphia_omni$aphiaID)))
      NonOmni_sum <- rowSums(NonOmni)
    
    # bind all sum and export
      t <- data.frame(rbind(auscpr$Sample_ID,FF_sum,Omni_sum,Carni_sum,NonFF_sum, NonCarni_sum, NonOmni_sum)) %>% 
        t() %>% 
        data.frame() %>% 
        rename(Sample_ID=V1)
      write_csv(t, "Output/data/df_sums_22052025.csv",row.names = FALSE)
      
    ##
    
      
    # bind all proportions
    # t <- data.frame(rbind(auscpr$Sample_ID,FF_prop,Omni_prop,Carni_prop,NonFF_prop)) %>% 
    #   t() %>% 
    #   data.frame() %>% 
    #   rename(c(Sample_ID=V1,Filter_prop = FF_prop,NonFilter_prop = NonFF_prop))
    
    # write.csv(t, "Output/data/df_proportions.csv",row.names = FALSE)
    #to double check the proportion (should be < 2 as total = (Carni + Omni + FF) = 1)
    t.tot <- colSums(t[2:5,])
    y <- t[2:5,]
    colSums(y)
    
    
    colSums(t[,2:5])
    t <- t %>% 
      rbind(t.tot) 
    
      mutate(prop_filter = rowSums(select(.,matches(aphia_FF$aphiaID)))/rowSums(auscpr_dat),
             prop_omni = rowSums(select(.,matches(aphia_omni$aphiaID)))/rowSums(auscpr_dat),
             prop_carn = rowSums(select(.,matches(aphia_carni$aphiaID)))/rowSums(auscpr_dat))
    
             
    df <- auscpr %>% #need to update the taxon names based on 'taxon alignment sheet'
      t() %>% 
      data.frame() %>% 
      rownames_to_column(var = "taxon") %>% 
      left_join(aphia_aus, by="taxon") %>% 
      relocate(aphiaID, .after=taxon) %>% 
      left_join(traits, by="aphiaID") %>% #occurence of many-to-many relationships
      t() %>% 
      data.frame()
    
    #21.04

    ####mine####
      df.t <- auscpr %>% 
        column_to_rownames("Sample_ID") %>% 
        t() %>% 
        data.frame()
      
      df.c <- df.t %>% 
        rownames_to_column(var = "taxon") %>% 
        left_join(traits, by = ("taxon")) #note resulted to an increase in objects
      #sum of filter-feeders 
      
      
      FF.total <- function(zoop_abund){
        zoop_abund %>% 
          filter(filter_or_not %in% "1") %>% 
          select(-c(taxonID, taxon, aphiaID, traitValue, FG_order)) %>% 
          colSums()
      
        zoop_abund %>%
          filter(filter_or_not %in% "0") %>%
          
        }  
          # summarise(across(everything(),~(.x)[colSums()]),
          #           rowname = 'Filter_prop') %>% 
          # bind_rows(zoop_abund %>% rownames_to_column()) %>% 
          # column_to_rownames("rowname")
      
      
        df.c2 <- FF.total(df.c) %>% 
          data.frame()
      
        
    ####second mine####
      #  next step is to try it in with the plankton taxon as columns
      
        df <- auscpr %>% #need to update the taxon names based on 'taxon alignment sheet'
          t() %>% 
          data.frame() %>% 
          rownames_to_column(var = "taxon") %>% 
          left_join(aphia_aus, by="taxon") %>% 
          relocate(aphiaID, .after=taxon) %>% 
          left_join(traits, by="aphiaID") %>% #occurence of many-to-many relationships
          
          
      #
        
      FF.total <- df.t %>% 
        filter(traitName == "trophicGroup") %>% 
        select(taxonID, stageID, scientificName, order, family, genus, traitValue) %>% 
        distinct() %>% 
        mutate(traitValue = as.character(traitValue))
      
      df.t2 <- df.t %>% 
        map(colnames(.), mutate(filter_tot = )
            %>% mutate(prop_filter = compute(filter_tot/(filter_tot+non_filter_tot))))
    
         #20.04
    # aphia.t <- aphia_aus %>% 
    #   t() %>% 
    #   data.frame() %>% 
    #   add_column(a = "a") %>% 
    #   relocate(a, .before = X1) 
    # 
    # rownames(aphia.t) <- 
    #df <- rbind(aphia.t, auscpr)
    
    # auscpr_header <- read.csv("CPR on-process/cpr_aus_complete.csv", sep=",", header=F, nrows=1) %>% 
    #   select(-c(V1,V3:V16)) %>% 
    #   data.frame()
    
    df <- rbind(auscpr_header, auscpr)
    #
    
    auscpr.t <- t(auscpr) %>% 
      data.frame() %>% 
      add_column(taxon = auscpr_header) %>% 
      relocate(taxon, .before = X1) %>% 
      left_join(aphia_aus, by="taxon") %>% 
      relocate(aphiaID, .after = taxon)
   
    #compute for proportion of filter-feeders per coordinate
    
    
    # auscpr_header <- read.csv("CPR on-process/cpr_aus_complete.csv", sep=",", header=F, nrows=1) %>% 
    #  select(-c(V1:V16)) %>% 
    #   t() %>% 
    #   data.frame()
    
    
    #19.04
    df.t <- auscpr %>% 
      select(-c(TripCode, Region, Latitude, Longitude, SampleTime_UTC, SampleTime_Local, Year_Local, Month_Local, Day_Local, Time_Local24hr, SatSST_degC, SatChlaSurf_mgm3, PCI, BiomassIndex_mgm3, SampleVolume_m3)) %>% 
      t()    
   
      
      df.t2 <- df.t[-1,] %>% 
        data.frame() %>% 
        add_column(taxon = auscpr_header)
        
        relocate(taxon, .before = X534.1) %>% 
        left_join(aphia_aus, by="taxon") %>% 
        relocate(aphiaID, .after = taxon) %>% 
        left_join(traits, by = "aphiaID") 
        
    df.t <- df.t %>% 
      inner_join(aphia_aus, by = "taxon")
    str(df.t)
    str(aphia_aus)
    #combine abundance and aphia
    
    
    #combine abundance data by aphia
    
    
    #match trait list to the abundance-aphia table
    # abundance is by column
    # trait is by row
    
    #compute for relative abundance
    
    # i just want to get percentage of feeding group in each coordinate
    
    #attach relative abundance to the cpr coordinate data frame
    
##






    
##
    
####Combining variables altogether into a dataframe####
    
    #for updating 
      df <- read_rds("Output/data/auscpr/df_complete_22052025.rds") %>% 
        select(-c("FF_sum","Carni_sum","Omni_sum","NonFF_sum","NonCarni_sum","NonOmni_sum"))
      
      #read updated data
      
      #update df_complete
      df <- df %>% 
        left_join(traits, by = "Sample_ID")
      #save updated df_complete
      write_rds(df, "Output/data/auscpr/df_complete_22052025.rds")
      
    #load aus cpr
    ##Keep TripCode, Sample ID, survey, region, latitude, longitude, date, BiomassIndex
      df <- read_csv("Input/CPR on-process/cpr_aus_complete.csv", col_names = TRUE) %>% 
        select(TripCode, Sample_ID, Latitude, Longitude, SampleTime_UTC)
      
    #load random effects
    ## Tow No., Tow Days, Longhurst
    #01 assign Tow_No
      df <- df %>% 
        group_by(TripCode) %>% 
        mutate(Tow_No = paste("AusCPR", cur_group_id(), sep = "")) #Instead of Region, update it to 'SurveyName' (AUS)
      
        
    #02 compute Tow_Days
      df <- df %>% 
        mutate(time = as.POSIXct(SampleTime_UTC, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")) %>% 
        group_by(TripCode) %>% 
        mutate(Tow_Days = as.double(difftime(time, min(time), units = "days"))) %>% 
        ungroup() %>% 
        select(-c(time))
    
    #03 intersect Longhurst Provinces
    #load longhurst provinces
      LH <- st_read(dsn = "Input/LonghurstProvinces/", layer = "Longhurst_world_v4_2010")
      #downloaded from https://marineregions.org/downloads.php
      
      df <- df %>% 
        st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>% 
        st_transform(crs = st_crs(LH)) 
      
        st_is_valid(LH)
        LH <- st_make_valid(LH)
      
      coords <- read_csv("Input/CPR on-process/cpr_aus_complete.csv", col_names = TRUE) %>% 
        select(Sample_ID, Latitude, Longitude)
      df <- st_intersection(df, LH) %>% 
        select(-c(ProvDescr, TripCode, geometry)) %>% 
        left_join(coords, by="Sample_ID") %>% 
        relocate(c("Longitude","Latitude"), .after=(Sample_ID)) %>% 
        rename("Longhurst" = "ProvCode")
      
    #load chlorophyll-a
    ## Extracted chlorophyll matched by time and space
      chl <- read_rds("Output/data/chla_extracted_long_09052025.rds")
              
      df <- df %>% 
        left_join(chl, by = "Sample_ID") %>% 
        rename(Chla = value) %>% 
        select(-name)
      
      
    #load proportions
    ## Relative abundances of filter-feeders, omnivores, and carnivores 
      #read df1 below first
      traits <- read_csv("Output/data/auscpr/df_sums_27052025.csv", col_names =T) %>% 
        rename(Sample_ID = V1)
      
      df <- df %>% 
        left_join(traits, by = "Sample_ID")
      
      df$geometry = NULL
     
      write_rds(df, "Output/data/df_complete_27052025.rds")
      
#Combining the global CPR datasets
      
      auscpr <- read_rds("Output/data/auscpr/df_complete_27052025.rds") %>% 
        mutate(Survey = "AusCPR")
      socpr <- read_rds("Output/data/socpr/df_complete_27052025.rds") %>% 
        mutate(Survey = "SOCPR")
      npacific <- read_rds("Output/data/npacific/df_complete_27052025.rds") %>% 
        mutate(Survey = "NPacific")
      natlantic <- read_rds("Output/data/natlantic/df_complete_27052025.rds") %>% 
        mutate(Survey = "NAtlantic") 
      
      df <- auscpr %>% 
        rows_insert(socpr, by = "Sample_ID") %>% 
        rows_insert(npacific, by ="Sample_ID") %>% 
        rows_insert(natlantic, by = "Sample_ID")
      df <- df %>% 
        mutate(Chla_sqrt = sqrt(Chla)) %>%
        mutate(ROC = Cope_omni_sum/(Cope_omni_sum + Cope_carni_sum)) %>% 
        mutate(RCO = Cope_carni_sum/(Cope_omni_sum + Cope_carni_sum)) %>% 
        mutate(RFF = FF_sum/(FF_sum+NonFF_sum)) %>% 
        mutate(ROC_Cal = Cal_omni_sum/(Cope_omni_sum + Cope_carni_sum)) %>% 
        mutate(RCO_Cal = Cal_carni_sum/(Cope_omni_sum + Cope_carni_sum)) %>% 
        mutate(ROC_Cyc = (Cyc_omni_sum+Har_omni_sum)/(Cope_omni_sum + Cope_carni_sum)) %>% 
        mutate(RCO_Cyc = Cyc_carni_sum/(Cope_omni_sum + Cope_carni_sum)) 
      
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
      
      #not impt as of 27052025  
      df <- df %>% 
        mutate(Chla_withCeilingAt5 = case_when(
          Chla >= 5 ~ 5,
          .default = Chla
        )) %>% 
        
        mutate(Chla_withCeilingAt8 = case_when(
          Chla >= 8 ~ 8,
          .default = Chla
        )) %>% 
        
        mutate(Chla_withCeilingAt12 = case_when(
          Chla >= 12 ~ 12,
          .default = Chla
        )) %>% 
        mutate(Chla_log = log(Chla + 1)) %>% 
        relocate(Survey, .after = Sample_ID) %>% 
        relocate(Chla_withCeilingAt5, .after = Chla) %>% 
        relocate(Chla_withCeilingAt8, .after = Chla_withCeilingAt5) %>% 
        relocate(Chla_withCeilingAt12, .after = Chla_withCeilingAt8) %>% 
        relocate(Chla_log, .after = Chla_withCeilingAt12)
      
      #Import eight day and combine
      auscpr_eightday <- read_rds("Output/data/chla_extracted_long_09052025.rds")
      socpr_eightday <- read_rds("Output/data/socpr/chla_extracted_long_socpr.rds")
      npacific_eightday <- read_rds("Output/data/npacific/chla_extracted_long_npacific.rds")
      natlantic_eightday <- read_rds("Output/data/natlantic/chla_extracted_long_natlantic.rds")
      
      chla_df_eightday <- auscpr_eightday %>% 
        rows_insert(socpr_eightday, by = "Sample_ID") %>% 
        rows_insert(npacific_eightday, by = "Sample_ID") %>% 
        rows_insert(natlantic_eightday, by = "Sample_ID") %>% 
        rename(Chla_8day = value) %>% 
        select(-c(name))
      
      df <- df %>% 
        left_join(chla_df_eightday, by="Sample_ID") 
      
      df <- df%>% 
        mutate(Chla_eightday_log = log(Chla_8day + 1)) %>% 
        relocate(Chla_8day, .after=Chla_log) %>% 
        relocate(Chla_eightday_log, .after = Chla_8day)
      
      #finalize df without Group
      df <- df %>% 
        relocate(Survey, .before= Sample_ID) %>% 
        relocate(Chla_sqrt, .after=Chla)
      
      write_rds(df, "Output/data/global/df_complete_27052025.rds")   

      #subdivide data into cyclopoid and calanoid
      
      df_cal <- df %>% 
        select(c(Survey,Sample_ID,Longitude, Latitude, SampleTime_UTC, Tow_No, Tow_Days, Longhurst, Chla, Chla_sqrt, ROC_Cal,ROC_Cal_SVT,ROC_Cal_SVT_zib,RCO_Cal,RCO_Cal_SVT,RCO_Cal_SVT_zib)) %>% 
        mutate(Group = "Calanoids") %>% 
        rename(ROC = ROC_Cal) %>% 
        rename(ROC_SVT = ROC_Cal_SVT) %>% 
        rename(ROC_SVT_zib = ROC_Cal_SVT_zib) %>% 
        rename(RCO = RCO_Cal) %>% 
        rename(RCO_SVT = RCO_Cal_SVT) %>% 
        rename(RCO_SVT_zib = RCO_Cal_SVT_zib)
      
      df_cyc <- df %>% 
        select(c(Survey,Sample_ID,Longitude, Latitude, SampleTime_UTC, Tow_No, Tow_Days, Longhurst, Chla, Chla_sqrt,ROC_Cyc,ROC_Cyc_SVT,ROC_Cyc_SVT_zib,RCO_Cyc,RCO_Cyc_SVT,RCO_Cyc_SVT_zib)) %>% 
        mutate(Group = "Cyclopoids") %>% 
        rename(ROC = ROC_Cyc) %>% 
        rename(ROC_SVT = ROC_Cyc_SVT) %>% 
        rename(ROC_SVT_zib = ROC_Cyc_SVT_zib) %>% 
        rename(RCO = RCO_Cyc) %>% 
        rename(RCO_SVT = RCO_Cyc_SVT) %>% 
        rename(RCO_SVT_zib = RCO_Cyc_SVT_zib)

      df_withGroup <- bind_rows(df_cal, df_cyc)    
      
      write_rds(df_withGroup, "Output/data/global/df_complete_withGroup_27052025.rds")
      