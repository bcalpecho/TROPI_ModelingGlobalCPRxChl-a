# ---
# title: Wrangling CPR and trait table together
# author: Pata & Hunt (2023))
# date: 2025-
# output: csv file
# ---

##Adapted from
#'7-Database_subset_template'
# ---
# title: "Template for subsetting the database"
# author: "Patrick Pata"
# date: "06/30/2023"
# output: html_document
# Repository: https://github.com/Pelagic-Ecosystems/Zooplankton_trait_database
# ---
##

#### Load libraries and data ####
    
    packages <- c("tidyverse",
                  "openxlsx")
    setwd("C:/Users/power/OneDrive - Université Libre de Bruxelles/4 UQ/CurrentWork/functional traits/Pata & Hunt (2023) code")
    # Function to download the packages if necessary. Otherwise, these are loaded.
    package.check <- lapply(
      packages,
      FUN = function(x)
      {
        if (!require(x, character.only = TRUE))
        {
          install.packages(x, dependencies = TRUE,
                           repos = "http://cran.us.r-project.org")
          library(x, character.only = TRUE)
        }
      }
    )
    
    source("toolkit.R")
    
    s.format <- read.csv("data_input/trait_dataset_standard_format_20230628.csv")[-1,]
    
    trait.directory <- read.csv("data_input/trait_directory_20230628.csv") %>% 
      distinct(traitID, .keep_all = TRUE)
    
    # taxonomy table
    taxonomy <- read.csv("data_input/taxonomy_table_20230628.csv") 
    
    # stage table
    lifestagelist <- read.csv("data_input/lifestage_directory_20230628.csv") %>% 
      select(-c(majorgroup, notes))
    
    # # Load Level 1 dataset
    traits.lvl1 <- read.csv("data_input/Trait_dataset_level1/trait_dataset_level1-2023-06-28.csv")
    
    # Load Level 2 dataset
    traits.lvl2 <- read.csv("data_input/Trait_dataset_level2/trait_dataset_level2-2023-09-14.csv")
    
    # Change assocTemperature column class to allow row binds
    traits.lvl2$assocTemperature <- as.character(traits.lvl2$assocTemperature)
    traits.lvl2$verbatimTemperature <- as.character(traits.lvl2$verbatimTemperature)

#### List selected species and traits ####

  # List of traits to extract
  # trait.list <- c("bodyLengthMax","dryWeight","carbonWeight","respirationRate_15C",
  #                 "respirationRate_WSC_15C","nitrogenTotal","excretionRateN_15C",
  #                 "trophicGroup","feedingMode")
  # trait.list <- c("trophicGroup", "feedingMode","reproductionMode","verticalDistribution","dielVerticalMigration")
    trait.list <- c("trophicGroup")
  
  
  # List of species to extract
  # Need to know the species name that matches the species names in this database. Alternatively, can do this matching using AphiaIDs and the acceptedNameUsageID field.
  
  #unique aphia ID pre-25.04.2025 (792 total unique aphia)
    species.list <- read.csv("C:/Users/power/OneDrive - Université Libre de Bruxelles/4 UQ/CurrentWork/db/CPR on-process/cpr_all_aphia.csv", sep=",", header=T)
  
  # Match species list with taxonomy file to get taxonID
    # species.list <- species.list  %>%
    #   left_join(taxonomy, by = c("aphiaID")) %>% 
    #   filter(!is.na(taxonID)) %>% 
    #   distinct()
    
  # Match species list with taxonomy file to get taxonID AND keep aphiaID without taxonID
    species.list <- species.list  %>%
      left_join(taxonomy, by = c("aphiaID")) %>% 
      distinct(aphiaID, .keep_all = TRUE)
      # double check  
        unique <- unique(species.list$aphiaID) # 792 unique aphiaID

#continue to subset 
##following lists are made to quantify the classes of my taxa (for recording purposes)##
    taxon.list <- read.csv("C:/Users/power/Desktop/db/CPR on-process/cpr_all_taxonList.csv", sep=",", header=T) 
    
    #'NA' aphiaID - 12
    taxon.list2 <- taxon.list %>% 
      filter(is.na(aphiaID)) 
    write.csv(taxon.list2, "data_output/taxon_list_NA-aphiaID.csv", row.names = FALSE)
    
    #distinct non-'NA' aphiaID -> 792
    taxon.list3 <- taxon.list %>%   
      filter(!is.na(aphiaID)) %>% 
      distinct(aphiaID, .keep_all = T)
    
    #with aphiaID and non-distinct taxonID - 688 (inc. NA)
    taxon.list1.5 <- taxon.list %>% 
      left_join(taxonomy, by = c("aphiaID")) %>% #many-to-many relationship between taxonID and aphiaID
      filter(!is.na(taxonID)) %>% 
      distinct(taxonID, .keep_all = T)
    
    #with aphia ID and 'NA' taxonID - 408
    taxon.list1.5 <- taxon.list %>% 
      left_join(taxonomy, by = c("aphiaID")) %>% #many-to-many relationship between taxonID and aphiaID
      filter(is.na(taxonID))
    
    #distinct aphiaID and 'NA' taxonID -> 242
    taxon.list4 <- taxon.list %>%   
      filter(!is.na(aphiaID)) %>% 
      distinct(aphiaID, .keep_all = T) %>% 
      left_join(taxonomy, by = c("aphiaID")) %>%
      filter(is.na(taxonID)) 
    write.csv(taxon.list3, "data_output/taxon_list_NA-taxonID.csv", row.names = FALSE)
    
    #distinct aphiaID and non-distinct non-'NA' taxonID -> 687
    #distinct aphiaID and distinct non-'NA' taxonID -> 550
    taxon.list5 <- taxon.list %>%   
      filter(!is.na(aphiaID)) %>% 
      distinct(aphiaID, .keep_all = T) %>%
      left_join(taxonomy, by = c("aphiaID")) %>%
      filter(!is.na(taxonID)) %>% 
      distinct(taxonID, .keep_all = T) 
    write.csv(taxon.list5, "data_output/taxon_list_distinct-taxonID.csv", row.names = FALSE)
    
    #distinct aphiaID
    taxon.list6 <- taxon.list %>% 
      distinct(aphiaID, .keep_all = T) %>%
      left_join(taxonomy, by = c("aphiaID"))

##

# Subset the species and traits 
# This will return the trait values based on literature or derived using related traits.

      # species.list <- taxon.list5 #using list of distinct aphiaID and distinct non-'NA' taxonID
      # species.list <- read.csv("data_output/taxon_list_distinct-taxonID.csv", header = T, sep=",")
      
    # Create a subset of records based on the species and traits list.
    trait.subset <- traits.lvl2 %>% 
      filter(traitName %in% trait.list) %>% 
      filter(taxonID %in% species.list$taxonID) %>% 
      arrange(scientificName)
    
    # Evaluate how many species have trait information
    trait.subset.perc <- trait.subset %>% 
      distinct(traitName, taxonID) %>% 
      group_by(traitName) %>% 
      summarise(Nrecords = n(), .groups = "drop") %>% 
      mutate(Perc.sp = Nrecords / nrow(species.list) * 100)
    
    #TG.trait.subset
    TG.trait.subset <- trait.subset %>% 
      filter(traitName == "trophicGroup") %>% 
      select(taxonID, stageID, scientificName, order, family, genus, traitValue) %>% 
      distinct() %>% 
      mutate(traitValue = as.character(traitValue))
    write.csv(TG.trait.subset, "data_output/TG_trait_subset.csv", row.names = FALSE)

##

### Filling up the gaps ###

  #### Method 1:  Individually per species ####

  traits.lvl2.ed <- traits.lvl2
    
  # A. Trophic Group
    # For each categorical trait, find the binary versions
    cat.missing <- species.list %>% 
      filter(taxonID %notin% filter(trait.subset, traitName == "trophicGroup")$taxonID)
    # List the binary version of traits assigned to a categorical trait
    binary.list <- traits.lvl2.ed %>% 
      filter(grepl("TG.",traitName)) %>% 
      distinct(traitName)
    # Select level of generalization
    gen.level <- "order"
    # Loop through each of the binary traits.
    traits.generalized <- s.format
    for (trait in binary.list) {
      for (taxon in cat.missing$scientificName){
        gen <- getGroupLevelValue(taxon, trait, gen.level = gen.level, 
                                  trait.df = traits.lvl2, 
                                  taxonomy.df = taxonomy)
        if(is.data.frame(gen)) {
          gen <- gen %>%  
            mutate(traitValue = as.character(traitValue))
          traits.generalized <- traits.generalized %>% 
            bind_rows(gen)
        }
      }
    }

  # B. Feeding Mode
    # For each categorical trait, find the binary versions
    cat.missing <- species.list %>% 
      filter(taxonID %notin% filter(trait.subset, traitName == "feedingMode")$taxonID)
    # List the binary version of traits assigned to a categorical trait
    binary.list <- traits.lvl2 %>% 
      filter(grepl("FM.",traitName)) %>% 
      distinct(traitName)
    # Select level of generalization
    gen.level <- "order"
    # Loop through each of the binary traits.
    for (trait in binary.list) {
      for (taxon in cat.missing$scientificName){
        gen <- getGroupLevelValue(taxon, trait, gen.level = gen.level, 
                                  trait.df = traits.lvl2, 
                                  taxonomy.df = taxonomy)
        if(is.data.frame(gen)) {
          gen <- gen %>%  
            mutate(traitValue = as.character(traitValue))
          traits.generalized <- traits.generalized %>% 
            bind_rows(gen)
        }
      }
    }
    
    # Update the overall working trait dataset
    traits.lvl2.ed <- bind_rows(traits.lvl2, traits.generalized)


  #### Method 2: Generalize at broader taxonomic level ####
  cat.missing <- species.list %>% 
    filter(taxonID %notin% filter(trait.subset, traitName == "trophicGroup")$taxonID) %>% 
    mutate(traitName = "trophicGroup")
  
  # Generalize trophic group at family level. 
  family.bin.tg <- traits.lvl2 %>% 
    # Need to know the prefix of the binary trait.
    filter(grepl("TG.",traitName)) %>% 
    mutate(traitValue = as.numeric(traitValue)) %>% 
    # TODO: Choose the taxonomic level
    group_by(family, traitName) %>% 
    # Collect all instances of a categorical level per group
    summarise(traitValue = sum(traitValue),  
              .groups = "drop") %>% 
    group_by(family) %>% 
    mutate(ntaxa = sum(traitValue)) %>% 
    ungroup() %>% 
    
    # convert to proportions
    mutate(traitValue = traitValue / ntaxa)
    #write.csv(family.bin.tg, "data_output/family_bin_tg.csv", row.names = FALSE)
  
  missing.cat.tg.fam <- cat.missing %>% 
    filter(traitName %in% c("trophicGroup")) %>% 
    distinct(species, family) %>% 
    left_join(family.bin.tg, by = "family") %>% 
    filter(!is.na(traitValue)) %>% 
    select(-ntaxa) 

  #Decision point for traitValue if the proportion of family > 0.5 
  family.bin.tg2 <- family.bin.tg %>% 
    filter(traitValue > 0.5) %>% 
    select(family, traitName, traitValue) %>% 
    mutate(traitValue = ifelse(traitValue > 0.5, "yes", "no")) %>% 
    distinct()
    
  TraitValueByFamily_Generalized <- cat.missing %>% 
    filter(traitName %in% c("trophicGroup")) %>% 
    distinct(family) %>% 
    left_join(family.bin.tg2, by = "family")

#integrate the values of generalized traitValues to the trait subset
  #blank trait
  blank.trait.tg <- s.format 
  
  #add blank trait if the family is in the 'TraitValueByFamily_Generalized'
  library(dplyr)
  family.gen.list <- TraitValueByFamily_Generalized %>% 
    filter(traitValue == "yes") %>% 
    select(family, traitName)
  
  #### Evaluate and finalize the species by traits table ####
  # Regenerate the trait subset from the updated dataset
  trait.subset.missing.TG <- traits.lvl2.ed %>% 
    filter(taxonID %in% cat.missing$taxonID) %>% 
    filter(family %in% family.gen.list$family) 
  #none of them is in the list :( wrong because of missing 'family' attribute in 'family.gen.list'
  
  #Update trait list with generalized traits
  trait.subset.missing.TG.with.generalized <- trait.subset.missing.TG %>% 
    left_join(TraitValueByFamily_Generalized, by = c("family")) %>% #results to traitValue.y
    mutate(traitValue = ifelse(is.na(traitValue.x), traitValue.y, traitValue.x)) %>% 
    select(-c(traitValue.x, traitValue.y))
  
  # Evaluate progress of data gap estimation
  trait.subset %>% 
    distinct(traitName, taxonID) %>% 
    group_by(traitName) %>% 
    summarise(Nrecords = n(), .groups = "drop") %>% 
    mutate(Perc.sp = Nrecords / nrow(species.list) * 100)


##

library(dplyr)
#Integrate Pata & Hunt (2023) to species list

  import_TG.PataHunt <- function(species.list){
    # Read the trait dataset
    traits.lvl2.ed <- read.csv("data_output/TG_trait_subset.csv")
    
    # Select the relevant columns
    traits.lvl2.ed <- traits.lvl2.ed %>% 
      select(taxonID, scientificName, order, family, genus, traitValue) %>% 
      distinct()
    
    # Join trait with the species list
    species.list <- species.list %>% 
        left_join(traits.lvl2.ed, by = c("taxonID", "scientificName", "order", "family", "genus")) %>% 
        #rename column of traitValue to "Pata&Hunt.TG"
        mutate(traitValue = ifelse(is.na(traitValue), "not determined", traitValue))
    
    return(species.list)
  }

####Assign traits to the taxa list based on other references#

#order-level  
ff.order <- c("Pyrosomida","Salpida","Dolioida","Copelata","Aplousobranchia","Phlebobranchia","Stolidobranchia")
c.order <- c("Aphragmophora","Phragmophora","Carybdeida","Chirodropida","Anthoathecata","Leptothecata","Siphonophorae","Actinulida","Limnomedusae","Narcomedusae","Trachymedusae","Polypodiidea","Coronatae", "Rhizostomeae", "Semaeostomeae", "Stauromedusae", "Bivalvulida", "Multivalvulida", "Malacovalvulida", "Actiniaria", "Antipatharia", "Ceriantharia", "Corallimorpharia", "Scleractinia", "Zoantharia", "Malacacyonacea","Sclerocyonacea")
h.order <- c("Euphausiacea")

assign_FG.order <- function(zoop.list){
  zoop.list %>% 
    mutate(FG_order = case_when(order %in% ff.order ~ "Filter-feeders",
                                   order %in% c.order ~ "Carnivores",
                                   order %in% h.order ~ "Herbivores",
                                   .default = "not determined"
    )) %>% 
    relocate(FG_order, .after = traitValue)
}

#family-level
# c.family 
# h.family
# ff.family
# 
# assign_FG.family <- function(zoop.list){
#   zoop.list %>% 
#     mutate(FG_family = case_when(family %in% c.family ~ "Carnivores",
#                                  family %in% h.family ~ "Herbivores",
#                                  .default = "not determined")) %>% 
#     relocate(FG_family, .after = FG_order)
# }

#filter-or-not
assign_filter_or_not <- function(zoop.list){
  zoop.list %>% 
    mutate(filter_or_not = case_when(order %in% ff.order ~ "1",
                                     .default = "0")) %>% 
    relocate(filter_or_not, .after = FG_order)
}

species.list <- species.list %>% 
  import_TG.PataHunt() %>% 
  assign_FG.order() %>% 
  #assign_FG.family() %>% 
  assign_filter_or_not()

##Final Output
#match the cpr taxon list to be able to identify aphia ID without taxon ID by their cpr's given taxon
  taxonlist <- read.csv("C:/Users/power/OneDrive - Université Libre de Bruxelles/4 UQ/CurrentWork/db/CPR on-process/cpr_all_aphia-taxon.csv", sep=",", header=T)
    str(taxonlist)  
    str(species.list)
  species.list <- species.list %>% 
    left_join (taxonlist, by = c("aphiaID")) %>% 
    relocate(taxon, .after = acceptedNameUsage) %>% 
    rename(taxon_cpr_list = taxon)
  ?rename
  write.csv(species.list, "data_output/species_list_with_traits_2025.26.04.csv", row.names = FALSE)

  #which are duplicated
  species.list %>% 
    which(aphiaID %in% aphiaID[which(duplicated(aphiaID))])
  

  summary(duplicated(species.list$aphiaID))
#Thaliacea (3), Appendicularia (1), Ascidiacea (3)  
#chaetognaths (2), cubozoa (2), hydrozoa (3), hydrozoa (6), polypodiozoa (1), Scyphozoa (3), Staurozoa (1), Myxosporeae (2), Malacosporea (1), Hexacorallia (6), Octocorallia (3)
#Malacostraca (1)

#ctenophora, cnidaria, echinodermata, enteropneusta, hemichordata, mollusca, nematoda, nemertea, onychophora, phoronida, platyhelminthes, sipuncula
#order %in% c("Amphipoda","Isopoda","Ostracoda","Tanaidacea","Cumacea","Mysida","Euphausiacea","Decapoda") ~ "Carnivores",



#Next Step is to manually fill up the trait table

##
