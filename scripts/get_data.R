library(HMDHFDplus)
library(dplyr)
library(magrittr)
library(tibble)
library(purrr)

# Read environment variables containing HFD credentials
user <- Sys.getenv("HFD_user")
pass <- Sys.getenv("HFD_pass")

# Command argument for which lexis type to download.
cmd_arg <- commandArgs(trailingOnly = T)
lexis_type <- cmd_arg[1]

if(is.na(lexis_type)){
  stop("Please specify a Human Fertility Database Lexis code as a command line",
       " argument. Most commonly: \n",
       " 'VV' for Age-Cohort \n",
       " 'RR' for Age-Year")  
}

# Countries to download
country_codes <- c("GBRTENW", "FRATNP", "SWE","USA")


# for reproducibility - which date was the data are you using available!
dates <- map_df(country_codes,
                function(country) tibble(Country=country,
                                         HFD_date=getHFDdate(country)))

for (country in country_codes){
  
  # read in births and exposure data
  births_df <- readHFDweb(CNTRY = country,
                          item=paste0("births", lexis_type),
                          username = user, password=pass, fixup = T) %>%
    rename(Births=Total)
  

  expos_df <- readHFDweb(CNTRY = country,
                         item=paste0("expos", lexis_type),
                         username = user,
                         password=pass,fixup = T)
  
  # Call age reached during year (ARDY) "Age" for easier interpretation.
  if("ARDY" %in% colnames(expos_df)){
    expos_df %<>% rename(Age=ARDY)
    births_df %<>% rename(Age=ARDY)
  }
  
  path <- file.path("data",lexis_type,country)
  dir.create(path,recursive = T)
  saveRDS(expos_df, file.path(path, "expos_hfd.rds"))
  saveRDS(births_df, file.path(path, "births_hfd.rds"))
}

saveRDS(dates, file.path("data", paste0("HFD_date_", lexis_type,".rds")))

## get all data  ---------------------------------------------------------------
# for Myrskyla et al. and Schmertmann et al. models
all_countries <- HMDHFDplus::getHFDcountries()

# wrapper around HFDweb (adding country column)
get_hfd_item <- function(country,item, user,pass){
  asfr_df <- readHFDweb(CNTRY = country,
                        item=item,
                        username = user,password=pass,fixup = T) %>%
    mutate(Country=country)
  return(asfr_df)
}

asfrRR_all <- map_df(all_countries, get_hfd_item, 
                     "asfrRR", user, pass) # MM
asfrVV_all <- map_df(all_countries, get_hfd_item, 
                     "asfrVV", user, pass) # schmert
exposVV_all <- map_df(all_countries, get_hfd_item,
                      "exposVV", user, pass) # schmert

# chile has too few observations? 
asfrRR_all %<>% filter(Country !="CHL")
asfrVV_all %<>% filter(Country !="CHL")
exposVV_all %<>% filter(Country !="CHL")

saveRDS(asfrRR_all,"data/asfrRR_all.rds")
saveRDS(asfrVV_all,"data/asfrVV_all.rds")
saveRDS(exposVV_all,"data/exposVV_all.rds")

