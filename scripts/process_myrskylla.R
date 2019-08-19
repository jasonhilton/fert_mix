
library(purrr)
library(tidyr)
library(magrittr)
library(ggplot2)
library(tibble)
library(dplyr)

source("R/myrskylla_model_functions.R")

asfr_all <- readRDS("data/asfrRR_all.rds") %>% as_tibble()


## fit all models 
## filter ages!
fitted_df <- asfr_all %>% filter(Age %in% 14:49) %>% nest(-Country) %>% 
  mutate(model=map(data, fit_single_myrskylla_model, 5,10))

fitted_df %<>% mutate(sigma2=var(map(fitted_df$model, "epsilon") %>% unlist()))

fitted_df %<>%
  mutate(Predictions=pmap(fitted_df,
                          function(model,data,sigma2,...){
                            predict_myrskylla_model(model,data,sigma2, 10,5)
                          }))

fitted_df %<>% 
  mutate(Combined=map2(data, Predictions, 
                       function(data, Predictions){
                         left_join(data, 
                                   Predictions %>% 
                                    rename(ASFR_predictions=ASFR))
                         }))

assess_df <- fitted_df %>% select(Country, Combined) %>% unnest()

saveRDS(assess_df,"processed/MM_assess.rds")

# calculate RMSE

rmse_df <- assess_df %>% drop_na() %>% filter(abs(Quantile-0.5) < 0.01) %>%
  mutate(se= (ASFR_predictions-ASFR)**2) %>% 
  group_by(Country) %>% summarise(RMSE=sqrt(mean(se)))

# ggplot(rmse_df, aes(y=Country, x=RMSE)) + geom_point()  + 
#   scale_x_continuous(expand = expand_scale(mult=c(0,0.05)),limits = c(0,NA))
  

saveRDS(rmse_df,"processed/MM_assess.rds")


# coverage----------------------------------------------------------------------

coverage_df <- assess_df %>% drop_na() %>% select(-Interval) %>% 
  mutate(Quantile=paste0("Q", Quantile*100)) %>%
  spread(Quantile, ASFR_predictions) %>%
  mutate(I50=(ASFR>Q25 & ASFR<Q75),
         I80=(ASFR>Q10 & ASFR<Q90),
         I90=(ASFR>Q5 & ASFR<Q95),
         I95=(ASFR>Q2.5 & ASFR<Q97.5)) %>%
  group_by(Country) %>%
  summarise(I50=sum(I50)/n(),I80=sum(I80)/n(),I90=sum(I90)/n(), I95=sum(I95)/n())

coverage_df %<>% gather(Interval, Proportion, -Country)
# ggplot(coverage_df,aes(y=Country,x=Proportion,colour=Interval)) + geom_point()

saveRDS(coverage_df, "processed/MM_coverage.rds")
