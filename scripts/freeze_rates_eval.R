library(magrittr)
library(dplyr)
library(tidyr)
library(tibble)
library(rlang)
library(purrr)
source("R/data_functions.R")

country_codes <- c("GBRTENW", "FRATNP", "USA", "SWE")

T_holdback <-10
emp_rates_df <- get_empirical_rate_df()

emp_rates_df <- map_dfr(country_codes, function(cntry){
  er_df <- get_empirical_rate_df(country=cntry) %>% 
    mutate(Country=cntry)
  return(er_df)
})


emp_rates_df %<>% group_by(Country) %>%
  mutate(Jump_off=max(Year) - T_holdback)
emp_rates_freeze <- emp_rates_df %>% filter(Year==Jump_off)

emp_rates_freeze %<>% select(Age, Rate) %>% rename(Rate_freeze=Rate)


emp_rates_assess <- emp_rates_df %>% 
  filter(Year>Jump_off,
         Age>=14,
         Age<=50) %>% 
  left_join(emp_rates_freeze)

emp_rates_assess %<>% mutate(Bias_rate = Rate_freeze-Rate,
                             AE_rate = abs(Bias_rate),
                             SE_rate = Bias_rate**2,
                             PE_rate=AE_rate/(Rate_freeze+1e-06),
                             E_births=Rate_freeze * Exposure,
                             Bias_births= E_births-Births,
                             AE_births=abs(Bias_births),
                             PE_births=AE_births/Births,
                             SE_births=Bias_births**2)

freeze_assess_sum <- emp_rates_assess %>% 
  group_by(Country) %>% 
  summarise(WBias_rate=weighted.mean(Bias_rate,w=Births),
            WAE_rate=weighted.mean(AE_rate,w=Births),
            WRMSE_rate=sqrt(weighted.mean(SE_rate,w=Births)),
            WPE_rate=weighted.mean(PE_rate,w=Births),
            WE_births=weighted.mean(E_births,w=Births),
            WBias_births=weighted.mean(Bias_births,w=Births),
            WAE_births=weighted.mean(AE_births,w=Births),
            WPE_births=weighted.mean(PE_births,w=Births),
            WRMSE_births=sqrt(weighted.mean(SE_births,w=Births)),
            Bias_rate=mean(Bias_rate),
            AE_rate=mean(AE_rate),
            RMSE_rate=sqrt(mean(SE_rate)),
            PE_rate=mean(PE_rate),
            E_births=mean(E_births),
            Bias_births=mean(Bias_births),
            AE_births=mean(AE_births),
            PE_births=mean(PE_births),
            RMSE_births=sqrt(mean(SE_births)))
                               

freeze_assess_age <- emp_rates_assess %>% group_by(Age,Country) %>%
  summarise(WBias_rate=weighted.mean(Bias_rate,w=Births),
            WAE_rate=weighted.mean(AE_rate,w=Births),
            WRMSE_rate=sqrt(weighted.mean(SE_rate,w=Births)),
            WPE_rate=weighted.mean(PE_rate,w=Births),
            WE_births=weighted.mean(E_births,w=Births),
            WBias_births=weighted.mean(Bias_births,w=Births),
            WAE_births=weighted.mean(AE_births,w=Births),
            WPE_births=weighted.mean(PE_births,w=Births),
            WRMSE_births=sqrt(weighted.mean(SE_births,w=Births)),
            Bias_rate=mean(Bias_rate),
            AE_rate=mean(AE_rate),
            RMSE_rate=sqrt(mean(SE_rate)),
            PE_rate=mean(PE_rate),
            E_births=mean(E_births),
            Bias_births=mean(Bias_births),
            AE_births=mean(AE_births),
            PE_births=mean(PE_births),
            RMSE_births=sqrt(mean(SE_births)))



saveRDS(freeze_assess_sum, "data/freeze_assess_sum.Rds")
saveRDS(freeze_assess_age, "data/freeze_assess_age.Rds")
saveRDS(emp_rates_assess, "data/freeze_assess.Rds")
