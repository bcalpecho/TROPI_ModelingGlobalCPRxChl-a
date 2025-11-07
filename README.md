# Modeling Global CPR 

Repository for MSc thesis project on zooplankton trophic groups across the global productivity gradient.
The project used the combined data of CPR Surveys of AusCPR, North Atlantic CPR, North Pacific CPR, and Southern Ocean CPR. 

### Step-by-step preparation of data and modeling.

```/data_input``` contains the CPR abundance tables, trait table, and combined data of the CPR Surveys.
For the chlorophyll-a data, the netcdf files of OC-CCI can be downloaded through [OC-CCI](https://www.oceancolour.org/thredds/ncss/grid/CCI_ALL-v6.0-8DAY/dataset.html). 

* ```1_generate_traits.R``` assign trophic groups to provided zooplankton taxon list based on trait tables.
* ```2_extract_chla.R``` aggregates, extracts, and fill-up gaps of OC-CCI chlorophyll-a values.
* ```3_generate_completeDF.R``` finalizes the data frame composed of relative abundances, trait, and chlorophyll data.
* ```4_model_globalCPR.R``` generates the model predictions.
* ```5_predict_globalCPR.R``` substitutes the chl-a predictions under SSP scenarios to the model.
* ```0_wrangling_CPR.R``` additional functions for visualizing and preparing the zooplankton data.

We hope that this work provides insights into how zooplankton is likely to respond to ongoing and future accelerating climate-driven changes in phytoplankton abundance and cell size. 
If you have any questions or comments, please send a mail at bryan.alpecho@ulb.be. 
