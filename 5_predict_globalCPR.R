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
library(stars)

#01 to extract data for SSP predictions
  #input netcdf
  ncdf <- "data_input/chla/ETOPO_intpp_Omon_ACCESS-ESM1-5_r1i1p1f1_185001-210012_yearmonths.nc"
  dat <- read_ncdf(ncdf, var = "intpp", proxy = F, ignore_bounds = T) 
  ##data source
  #Cite: Buchanan, P. J. (2022). CMIP6 model vertically-integrated net primary production data (v1.0) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.7226186
  #temporal coverage: 1850-2100
  
  # Conversion from mol C/m2/s to global NPP (Pg C/y) 
  convert_to_globalNPP <- function(x){x * 136.57}
  
  # Conversion from global NPP (Pg C/y) to Chl from Maranon et al. 2014 (referred to by Campbell et al., 2021)
  convert_to_chl <- function(x) { 10^((log10(x) - 1.58)/1.29) }
  
  #extract data from 1870-1899
  ext <- st_apply(dat, 4, mean, na.rm=T)
  means <- ext$mean %>%  data.frame() %>% rename(NPP = ".") %>% mutate(year = 1850:2100) %>% mutate(NPP_global = convert_to_globalNPP(NPP)) %>% mutate(Chla = convert_to_chl(NPP_global)) %>% mutate(Chla_sqrt = sqrt(Chla))
  
  #plot means
  plot_annual <- ggplot(data = means, aes(x = year, y = chla)) + geom_point() + geom_smooth(method = 'lm')
  #1870-1899 [1] 2.007322e-07
  #2100 2.045731e-07
  
#02 to predict trophic group ratios using model for globalCPR
  #input final models
  Filter_mdl <- read_rds("Output/previousModels/mdl_export/Filter_mdl_zib.rds")
  Carni_mdl <- read_rds("Output/previousModels/mdl_export/Carni_mdl_zib_surveyLat.rds")
  Omni_mdl <- read_rds("Output/previousModels/mdl_export/Omni_mdl_zib_surveyLat.rds")
  
  #predictions
  pop_preds_FF<- predictions(Filter_mdl, 
                             newdata = datagrid(Chla_sqrt = means$chla,
                                                year = means$year),
                             re.form = NA, type = "response") 
  
  #plot
  ff_plot <- ggplot(data = pop_preds_FF) + theme_bw(base_size = 18) + 
    geom_ribbon(aes(x = year, y = estimate, 
                    ymin = conf.low, ymax = conf.high), alpha = 0.3, colour = "blue", linetype=2) + 
    geom_line(aes(x = year, y = estimate),colour="blue") + 
    labs(y = "Filter-feeders proportion", x = "Chl-a")  
  

#03 to predict using the 1870-1899 baseline value extracted from ACCESS-ESM1-5
  #determine the average 1870-1899
  baseline <- means %>% filter(year >= 1870 & year <= 1899) %>% summarise(mean = mean(chla))
  
  #low emission (SSP1-2.6)
  ssp_low <- -0.56
  
  #high emission (SSP5-8.5)
  ssp_high <- -2.99

  # Calculate future estimates of NPP
  new_NPP_SSP_high <- baseline + baseline * ssp_high/100
  new_NPP_SSP_low <- baseline + baseline * ssp_low/100
  
  #predictions
  ssp_preds_FF<- predictions(Filter_mdl, 
                             newdata = datagrid(Chla_sqrt = c(new_NPP_SSP_low, new_NPP_SSP_high)),
                             re.form = NA, type = "response") 
  
  ssp_preds_Omni <- predictions(Omni_mdl, 
                             newdata = datagrid(Chla_sqrt = c(new_NPP_SSP_low, new_NPP_SSP_high)),
                             re.form = NA, type = "response") 
  
  ssp_preds_Carni <- predictions(Carni_mdl, 
                             newdata = datagrid(Chla_sqrt = c(new_NPP_SSP_low, new_NPP_SSP_high)),
                             re.form = NA, type = "response") 
  
  #plot predictions from 1850 - 2100
  
  
##
  
  