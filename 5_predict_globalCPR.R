# ---
# title: Predict Global CPR
# author: Bryan Alpecho
# date: 2025-
# output: csv and plot for predictions
# ---

#Aims
#01 to compute for change in trophic group ratio per unit decline of chl-a
#02 to extract data for SSP predictions
#03 to predict trophic group ratios from 1850-2100 using model for globalCPR
#04 to predict % change by 2100 relative to the 1870-1899 baseline value extracted from ACCESS-ESM1-5

# Load libraries
library(tidyverse)
library(glmmTMB)
library(marginaleffects)
library(stars)

#set date for version control
    date <- "28102025"
    df_date <- "16092025" #df of final models

#Import data
    #final models
    load("Output/previousModels/revision/merged_chla_mdls.RData")
    
    #dataframe for final models
    df <- read_rds(paste("data_input/global/df_complete_",df_date,".rds",sep=""))
    df_noNA <- df %>% filter(!is.na(chla)) 
    
#01 to compute for change in trophic group ratio per unit decline of chl-a (mg m-3)
  
  ### Extremes of observed chl-a ###
  min(df$chla, na.rm=TRUE) #0.02847524
  max(df$chla, na.rm=TRUE) #19.26117
  
  extract_chla <- df %>% select(chla, chla_sqrt) %>% 
    summarise(Min_Chla = min(chla, na.rm=T),
      Max_Chla = max(chla, na.rm=T),
      Mean_Chla = mean(chla, na.rm=T),
      Min_Chla_sqrt = sqrt(min(chla, na.rm=T)),
      Max_Chla_sqrt = sqrt(max(chla, na.rm=T)),
      Mean_Chla_sqrt = sqrt(mean(chla, na.rm=T)))
  
  #Gelatinous filter-feeders
      #Predictions for chl-a
      #compute for rate of change across the range of observed chl-a values  
      pop_preds_FF<- predictions(Filter_mdl_zib, 
                             newdata = datagrid(chla_sqrt = seq(floor(min(df_noNA$chla_sqrt)),ceiling(max(df_noNA$chla_sqrt)),0.01)),
                             re.form = NA) # Zeros out random effects
  
      FF_max <- pop_preds_FF %>% filter(chla_sqrt == round(sqrt(max(df_noNA$chla)), 2))
      FF_min <- pop_preds_FF %>% filter(chla_sqrt == round(sqrt(min(df_noNA$chla)), 2))
      FF_rate <- abs(FF_max$estimate - FF_min$estimate)/abs(extract_chla$Max_Chla - extract_chla$Min_Chla) * 100
      FF_range <- (abs(extract_chla$Max_Chla - extract_chla$Min_Chla))*FF_rate
      FF_preds_chla_summary <- tibble(FF_max = FF_max$estimate, FF_min = FF_min$estimate, FF_rate = FF_rate, FF_range = FF_range)
      FF_preds_chla_summary
      
      #Predictions by survey
      #compute for difference in estimates across CPR surveys
      pop_preds_FF_surveys <- predictions(Filter_mdl_zib, 
                                          newdata = datagrid(survey = unique(df_noNA$survey), 
                                                             chla_sqrt = seq(floor(min(df_noNA$chla_sqrt)),ceiling(max(df_noNA$chla_sqrt)),0.05)), 
                                          re.form = NA) # Zeros out random effects
      
      estimates_surveys <- pop_preds_FF_surveys %>% select(survey, estimate) %>% 
        group_by(survey) %>% 
        summarise(mean(estimate)) %>% column_to_rownames("survey") %>% t() %>% data.frame()
      
      survey_summary <- tibble(auscprVSnatlantic = estimates_surveys$natlantic - estimates_surveys$auscpr, 
                               auscprVSnpacific = estimates_surveys$auscpr - estimates_surveys$npacific,
                               auscprVSsocpr = estimates_surveys$auscpr - estimates_surveys$socpr)
      
      unique(pop_preds_FF$survey)

      
  #Omnivorous zooplankton
      #Predictions for chl-a
      pop_preds_omni_chla <- predictions(Omni_mdl_zib, 
                                         newdata = datagrid(chla_sqrt = seq(floor(min(df_noNA$chla_sqrt)),ceiling(max(df_noNA$chla_sqrt)),0.01)), 
                                         re.form = NA) # Zeros out random effects
      
      Omni_max <- pop_preds_omni_chla %>% filter(chla_sqrt == round(sqrt(max(df_noNA$chla)), 2))
      Omni_min <- pop_preds_omni_chla %>% filter(chla_sqrt == round(sqrt(min(df_noNA$chla)), 2))
      Omni_rate <- abs(Omni_max$estimate - Omni_min$estimate)/abs(extract_chla$Max_Chla - extract_chla$Min_Chla) * 100
      Omni_range <- (abs(extract_chla$Max_Chla - extract_chla$Min_Chla))*Omni_rate
      Omni_preds_chla_summary <- tibble(Omni_max = Omni_max$estimate, Omni_min = Omni_min$estimate, Omni_rate = Omni_rate, Omni_range = Omni_range)
      Omni_preds_chla_summary
      
      #Predictions by survey
      pop_preds_omni_surveys <- predictions(Omni_mdl_zib, 
                                            newdata = datagrid(survey = unique(df_noNA$survey), 
                                                               chla_sqrt = seq(floor(min(df_noNA$chla_sqrt)),ceiling(max(df_noNA$chla_sqrt)),0.5)), 
                                            re.form = NA) # Zeros out random effects
      omni_estimates_surveys <- pop_preds_omni_surveys %>% select(survey, estimate) %>% 
        group_by(survey) %>% 
        summarise(mean(estimate)) %>% column_to_rownames("survey") %>% t() %>% data.frame()
      
      omni_survey_summary <- tibble(auscprVSnatlantic = abs(omni_estimates_surveys$natlantic - omni_estimates_surveys$auscpr), 
                               auscprVSnpacific = abs(omni_estimates_surveys$auscpr - omni_estimates_surveys$npacific),
                               auscprVSsocpr = abs(omni_estimates_surveys$auscpr - omni_estimates_surveys$socpr))
      
      
  #Carnivorous zooplankton
      #Predictions for chl-a
      pop_preds_carni_chla <- predictions(Carni_mdl_zib, 
                                    newdata = datagrid(chla_sqrt = seq(floor(min(df_noNA$chla_sqrt)),ceiling(max(df_noNA$chla_sqrt)),0.01)), 
                                    re.form = NA) # Zeros out random effects
      
      Carni_max <- pop_preds_carni_chla %>% filter(chla_sqrt == round(sqrt(max(df_noNA$chla)), 2))
      Carni_min <- pop_preds_carni_chla %>% filter(chla_sqrt == round(sqrt(min(df_noNA$chla)), 2))
      Carni_rate <- abs(Carni_max$estimate - Carni_min$estimate)/abs(extract_chla$Max_Chla - extract_chla$Min_Chla) * 100
      Carni_range <- (abs(extract_chla$Max_Chla - extract_chla$Min_Chla))*Carni_rate
      Carni_preds_chla_summary <- tibble(Carni_max = Carni_max$estimate, Carni_min = Carni_min$estimate, Carni_rate = Carni_rate, Carni_range = Carni_range)
      Carni_preds_chla_summary
      
      #Predictions by survey
      pop_preds_carni_surveys <- predictions(Carni_mdl_zib, 
                                            newdata = datagrid(survey = unique(df_noNA$survey), 
                                                               chla_sqrt = seq(floor(min(df_noNA$chla_sqrt)),ceiling(max(df_noNA$chla_sqrt)),0.5)), 
                                            re.form = NA) # Zeros out random effects
      carni_estimates_surveys <- pop_preds_carni_surveys %>% select(survey, estimate) %>% 
        group_by(survey) %>% 
        summarise(mean(estimate)) %>% column_to_rownames("survey") %>% t() %>% data.frame()
      
      carni_survey_summary <- tibble(auscprVSnatlantic = abs(carni_estimates_surveys$natlantic - carni_estimates_surveys$auscpr), 
                                    auscprVSnpacific = abs(carni_estimates_surveys$auscpr - carni_estimates_surveys$npacific),
                                    auscprVSsocpr = abs(carni_estimates_surveys$auscpr - carni_estimates_surveys$socpr)
                                    )
      
  #### #Campbell et al (2021)####
  #     # Calculate the proportional change
  #     1- y_hat[2]/y_hat[1]
  #     
  #     # Absolute change
  #     y_hat[2] - y_hat[1]
  #     
  #     # length per degree
  #     (y_hat[1] - y_hat[2])/ diff(df$Sqrt_chl)
  #     
  #     # Percentage change in Mass
  #     ww_y_hat <- convert_to_weight(y_hat)
  #     1 - ww_y_hat[2] / ww_y_hat[1]
  #     
  #     
  #     # Chl predictions for climate change --------------------------------------
  #     
  #     # Chlorophyll changes
  #     old_NPP <- 52.1 # As per Bopp et al. 2013
  #     
  #     # As per Bopp et al. 2013 (percent changes)
  #     RCP8.5_per <- -8.6 
  #     RCP2.6_per <- -2.0 
  #     
  #     # Calculate future estimates of NPP
  #     new_NPP_RCP8.5 <- old_NPP + old_NPP * RCP8.5_per/100
  #     new_NPP_RCP2.6 <- old_NPP + old_NPP * RCP2.6_per/100
  #     
  #     # Conversion to Chl from Maranon et al. 2014
  #     convert_to_chl <- function(x) { 10^((log10(x) - 1.58)/1.29) }
  #     
  #     # Absolute changes
  #     abs_chl_change_RCP8.5<- convert_to_chl(new_NPP_RCP8.5) - convert_to_chl(old_NPP)
  #     abs_chl_change_RCP2.6<- convert_to_chl(new_NPP_RCP2.6) - convert_to_chl(old_NPP)
  #     
  #     
  #     ### RCP8.5 predictions ###
  #     
  #     newdat <- with(CopeData, data.frame( Sqrt_chl = c(sqrt(convert_to_chl(old_NPP)),
  #                                                       sqrt(convert_to_chl(new_NPP_RCP8.5))),
  #                                          SST  = mean(SST),
  #                                          Asin_omni = mean(Asin_omni)))
  #     
  #     # Make predictions based on our model
  #     y_hat <- predict(hr_mod, newdata = newdat, re.form = NA, type = "response")
  #     
  #     # Percentage change in Mass
  #     ww <- convert_to_weight(y_hat)
  #     ww[2] / ww[1]
  #     
  #     
  #     ### RCP2.6 predictions ###
  #     
  #     newdat <- with(CopeData, data.frame( Sqrt_chl = c(sqrt(convert_to_chl(old_NPP)),
  #                                                       sqrt(convert_to_chl(new_NPP_RCP2.6))),
  #                                          SST  = mean(SST),
  #                                          Asin_omni = mean(Asin_omni)))
  #     
  #     # Make predictions based on our model
  #     y_hat <- predict(hr_mod, newdata = newdat, re.form = NA, type = "response")
  #     
  #     # Percentage change in Mass
  #     ww <- convert_to_weight(y_hat)
  #     ww[2] / ww[1]
  #     
  #     
  #     # Combined changes
  #     1-1*1.012*0.927 #RCP8.5
  #     1-1*1.003*0.984 #RCP2.6
    #### ####  

#02 to extract data for SSP predictions
  #input netcdf
  ncdf <- "data_input/chla/ETOPO_intpp_Omon_ACCESS-ESM1-5_r1i1p1f1_185001-210012_yearmonths.nc"
  dat <- read_ncdf(ncdf, var = "intpp", proxy = T, ignore_bounds = T) 
  ##data source
  #Cite: Buchanan, P. J. (2022). CMIP6 model vertically-integrated net primary production data (v1.0) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.7226186
  #temporal coverage: 1850-2100
  
  #23.08.2025
    #Context: SSP projections (anomalies) are provided as 2080-2099 mean values relative to the 1870-1899 mean (Kwiatkowski et al., 2010).
    
    #convert from NPP by mol C (C/m^2/s) to NPP by mass C (g C/m^2/year)  
    convert_to_nppPerYr <- function(x) { (x * 12.01)*60*60*24*365.25 }
    ext <- st_apply(dat, 1:4, convert_to_nppPerYr, PROGRESS = T)  
    
    #compute for means per year
    ext_means <- st_apply(ext, 4, mean, na.rm=T) 
    
    
    #get baseline value (1870-1899 mean)
    compute_NPP <- function(baseline) { 
      
      #get annual means 
      means <- baseline$mean %>%  data.frame() %>% 
        rename(chla = ".") %>% 
        mutate(year = 1850:2100) %>% 
        #mutate(NPP_global = convert_to_carbonFixation(NPP)) %>% 
        #mutate(chla = convert_to_chl(NPP_global)) %>% 
        mutate(chla_sqrt = sqrt(chla))
      
      baseline <- means %>% filter(year >= 1870 & year <= 1899) %>% summarise(mean = mean(chla))
      
      #low emission (SSP1-2.6)
      ssp_low <- -0.56
      
      #high emission (SSP5-8.5)
      ssp_high <- -2.99
      
      # Calculate future estimates of NPP given SSP projection
      new_NPP_SSP_high <- baseline + baseline * ssp_high/100
      new_NPP_SSP_low <- baseline + baseline * ssp_low/100
      
      NPP_histAndSSP <- tibble(
        baseline_NPP = baseline,
        highSSP_NPP = new_NPP_SSP_high,
        lowSSP_NPP = new_NPP_SSP_low
      )
      
      return(NPP_histAndSSP)
    }
    NPP_means <- compute_NPP(ext_means)
    #data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABIAAAASCAYAAABWzo5XAAAAbElEQVR4Xs2RQQrAMAgEfZgf7W9LAguybljJpR3wEse5JOL3ZObDb4x1loDhHbBOFU6i2Ddnw2KNiXcdAXygJlwE8OFVBHDgKrLgSInN4WMe9iXiqIVsTMjH7z/GhNTEibOxQswcYIWYOR/zAjBJfiXh3jZ6AAAAAElFTkSuQmCC
    #convert from NPP by mass C (g C/ m^2/ year) to chl-a (), with intermediate conversion to NPP in Maranon et al 2014 (mg C/m^3/day)
    #convert_to_NPPbyMassC <- function(x) { ((x/365.25)/50)*1000 }
    
    convert_to_chl <- function(x) { 10^((log10(x) - 1.58)/1.29) }
    
    #with assumption of 50-m depth for depth-integrated NPP estimates
    convert_to_NPP <- function(x){
      (((x/365.25)/50)*1000)
    }
    NPP_means$baseline_NPP$mean
    convert_to_NPP(76)
    chla_sqrt_baseline <- convert_to_chl(convert_to_NPP(NPP_means$baseline_NPP$mean))
    chla_sqrt_SSP_high <- convert_to_chl(convert_to_NPP(NPP_means$highSSP_NPP$mean))
    chla_sqrt_SSP_low <- convert_to_chl(convert_to_NPP(NPP_means$lowSSP_NPP$mean))
    
    
    
#### #21.08.2025 ####
  #   #convert from NPP of CMIP (mol C/m2/s) to NPP of Maranon (mg C/m3)
  #   convert_to_npp <- function(x) { (x * 12.01) }
  #   
  #   #KWiatkowski et al 2020 unit of NPP is g C m-2 y-1
  # 
  # 
  # #15.08.2025
  # 
  # ext_means_ver2 <- st_apply(dat, 4, mean, na.rm=T) 
  # 
  # convert_to_nppPerYr <- function(x) { (x * 12.01)*60*60*24*365.25 }
  # 
  # means <- ext_means_ver2$mean %>%  data.frame() %>% 
  #   rename(chla = ".") %>% 
  #   mutate(year = 1850:2100) %>% 
  #   #mutate(NPP_global = convert_to_carbonFixation(NPP)) %>% 
  #   #mutate(chla = convert_to_chl(NPP_global)) %>% 
  #   mutate(nppPerYr = convert_to_nppPerYr(chla))
  # 
  # baseline <- means %>% filter(year >= 1870 & year <= 1899) %>% summarise(mean = mean(nppPerYr))
  # 
  # #convert all the way to chla
  # 
  # convert_to_carbonFixation_ver2 <- function(x) { ((x/365.25)/50)*1000 }
  # #12.08.2025
  #   
  #   #12.01 g/mol C for molar mas  
  #   #assuming 50-m as depth of the intpp
  #   #60 * 60  * 24
  #   
  #   convert_to_carbonFixation <- function(x) { (((x * 12.01)/50)*60*60*24)*1000 } #to convert mol C/m2/s to mg C/m3/d
  #   
  #   # Conversion from global NPP (Pg C/y) to Chl from Maranon et al. 2014 (referred to by Campbell et al., 2021) #chl-a in unit of mg C/m3/day
  #   convert_to_chl <- function(x) { 10^((log10(x) - 1.58)/1.29) }
  #   
  # #extract data from 1870-1899 because average of log values is not equal to log of average of values
  # #convert from vertically-integrated NPP to carbon fixation rate
  # ext <- st_apply(dat, 1:4, convert_to_carbonFixation) 
  # 
  # 
  # #convert from carbon fixation rate to chl-a
  # ext_chl <- st_apply(ext, 1:4, convert_to_chl) 
  # #compute for means per year
  # ext_means <- st_apply(ext_chl, 4, mean, na.rm=T) 
  # #summarise by year
  # means <- ext_means$mean %>%  data.frame() %>% 
  #   rename(chla = ".") %>% 
  #   mutate(year = 1850:2100) %>% 
  #   #mutate(NPP_global = convert_to_carbonFixation(NPP)) %>% 
  #   #mutate(chla = convert_to_chl(NPP_global)) %>% 
  #   mutate(chla_sqrt = sqrt(chla))
  # 
  # baseline <- means %>% filter(year >= 1870 & year <= 1899) %>% summarise(mean = mean(chla))
  # convert_frm_chla_to_npp <- function(x) { 10^(1.29*(log10(x))+1.58) }
  # convert_frm_chla_to_npp(baseline)
  # 
  # 
  # 
  # #convert from 
  # #plot means
  # plot_annual <- ggplot(data = means, aes(x = year, y = chla)) + geom_point() + geom_smooth(method = 'lm')
  # 
  # #get the baseline value 1870-1899
  # 
  # 
  # #1870-1899 [1] 2.007322e-07
  # #2100 2.045731e-07
  
# #03 to predict trophic group ratios using model for globalCPR
#   #input final models
#   Filter_mdl <- read_rds("Output/previousModels/mdl_export/Filter_mdl_zib.rds")
#   Carni_mdl <- read_rds("Output/previousModels/mdl_export/Carni_mdl_zib.rds")
#   Omni_mdl <- read_rds("Output/previousModels/mdl_export/Omni_mdl_zib.rds")
#   
#   #predictions
#   pop_preds_FF<- predictions(Filter_mdl, 
#                              newdata = datagrid(Chla_sqrt = means$chla,
#                                                 year = means$year),
#                              re.form = NA, type = "response") 
#   
#   #plot
#   ff_plot <- ggplot(data = pop_preds_FF) + theme_bw(base_size = 18) + 
#     geom_ribbon(aes(x = year, y = estimate, 
#                     ymin = conf.low, ymax = conf.high), alpha = 0.3, colour = "blue", linetype=2) + 
#     geom_line(aes(x = year, y = estimate),colour="blue") + 
#     labs(y = "Filter-feeders proportion", x = "Chl-a")  
####  ####

#03 to predict using the 1870-1899 baseline value extracted from ACCESS-ESM1-5
  # #determine the average 1870-1899
  # baseline <- means %>% filter(year >= 1870 & year <= 1899) %>% summarise(mean = mean(nppPerYr))
  # 
  # #low emission (SSP1-2.6)
  # ssp_low <- -0.56
  # 
  # #high emission (SSP5-8.5)
  # ssp_high <- -2.99
  # 
  # # Calculate future estimates of NPP
  # new_NPP_SSP_high <- baseline + baseline * ssp_high/100
  # new_NPP_SSP_low <- baseline + baseline * ssp_low/100
  # 
  # carbonFixation_baseline <- convert_to_carbonFixation_ver2(baseline)
  # chla_sqrt_baseline <- sqrt(convert_to_chl(carbonFixation_baseline))
  # 
  # carbonFixation_SSP_high <- convert_to_carbonFixation_ver2(new_NPP_SSP_high)
  # chla_sqrt_SSP_high <- sqrt(convert_to_chl(carbonFixation_SSP_high))
  # 
  # carbonFixation_SSP_low <- convert_to_carbonFixation_ver2(new_NPP_SSP_low)
  # chla_sqrt_SSP_low <- sqrt(convert_to_chl(carbonFixation_SSP_low))
  
  #predictions
  ssp_preds_FF <- predictions(Filter_mdl_zib, 
                             newdata = datagrid(Chla_sqrt = c(chla_sqrt_baseline, chla_sqrt_SSP_low, chla_sqrt_SSP_high)),
                             re.form = NA, type = "response") 
  
  ssp_preds_Omni <- predictions(Omni_mdl_zib, 
                             newdata = datagrid(Chla_sqrt = c(chla_sqrt_baseline, chla_sqrt_SSP_low, chla_sqrt_SSP_high)),
                             re.form = NA, type = "response") 
  
  ssp_preds_Carni <- predictions(Carni_md_zib, 
                             newdata = datagrid(Chla_sqrt = c(chla_sqrt_baseline, chla_sqrt_SSP_low, chla_sqrt_SSP_high)),
                             re.form = NA, type = "response") 
  
  #plot predictions from 1850 - 2100
  
  
##
  
  