library(HMDHFDplus)
library(MASS)
library(purrr)
library(magrittr)
library(dplyr)
library(tidyr)
library(Matrix)
library(bfcf)
library(ggplot2)

asfrVV_all <-readRDS("data/asfrVV_all.rds")
exposVV_all <- readRDS("data/exposVV_all.rds")

# penalty sizes
#ages <- 15:44
ages <- 14:49
hist_df <- get_historical_data(asfrVV_all,ages = ages)

n_ages <- length(ages)
penalties <- get_historical_penalties(hist_df)


data_all <-left_join(asfrVV_all, exposVV_all)
# do calibration

start_cohort <- 1950
country <- "GBRTENW"
forecast_joy <- 2006
forecast_horizon <- 2016

asfr_GBRTENW_df <- data_all %>% filter(Cohort>=start_cohort,
                                         ARDY %in% ages,
                                         Country==country)

n_forecast <- forecast_horizon - forecast_joy

last_observed <- forecast_joy - min(ages)

cohorts <- start_cohort:(last_observed + n_forecast)
n_cohorts <- length(cohorts)

last_complete <- find_last_complete_cohort(asfr_GBRTENW_df, ages=ages) - n_forecast
ind_last_complete <- which(cohorts %in% last_complete)

Kjs <- make_K_matrices(n_cohorts,ind_last_complete,ages, penalties)

cohort_inds <- which(cohorts %in% ((last_complete+1):(last_observed + n_forecast)))
calib <- calibrate_penalties(Kjs, ages, 
                             cohort_inds,n_iters = 30)


post_df <- get_posterior(calib$K,
                         data_cntry_df =  asfr_GBRTENW_df,
                         forecast_joy,
                         forecast_horizon,
                         ages=ages)

# ggplot(post_df %>% filter(Cohort %in% c(1980,1985,1990,1995,2000)),
#        aes(x=ARDY,y=ASFR_prediction)) +
#   geom_line() + facet_wrap(~Cohort) +
#   geom_point(aes(y=ASFR),colour="red") +
#   geom_vline(aes(xintercept=2006-Cohort))+
#   xlim(14,49)
# 
# plot_ages <- c(18,23,28,33,38,43)
# ggplot(post_df %>% filter(ARDY %in% plot_ages),
#        aes(x=Cohort,y=ASFR_prediction)) +
#   geom_line() + facet_wrap(~ARDY) +
#   geom_point(aes(y=ASFR),colour="red") +
#   geom_vline(aes(xintercept=2006-ARDY))



# ggplot(post_df %>% filter(ARDY %in% plot_ages),
#        aes(x=Cohort,y=ASFR_prediction)) +
#   geom_ribbon(aes(ymax=p95,ymin=p5), alpha=0.2) +
#   geom_ribbon(aes(ymax=p75,ymin=p25), alpha=0.2) +
#   geom_line() + facet_wrap(~ARDY) +
#   geom_point(aes(y=ASFR),colour="red") +
#   theme_bw()+
#   geom_vline(aes(xintercept=2006-ARDY))


# ggplot(post_df %>% filter(Cohort %in% c(1980,1985,1990,1995,2000)),
#        aes(x=ARDY,y=ASFR_prediction)) +
#   geom_ribbon(aes(ymax=p95,ymin=p5), alpha=0.25) +
#   geom_ribbon(aes(ymax=p75,ymin=p25), alpha=0.5) +
#   geom_line() + facet_wrap(~Cohort) +
#   geom_point(aes(y=ASFR),colour="red") +
#   geom_vline(aes(xintercept=2006-Cohort))+
#   theme_bw() +
#   xlim(14,45)


post_df %<>%
  mutate(SE=(ASFR_prediction-ASFR)**2,
         I90=((ASFR<p95 &
                 ASFR>p5)),
         I50=((ASFR<p75 &
                 ASFR>p25)))

# ggplot(post_df %>% group_by(Cohort) %>%
#          summarise(RMSE=sqrt(mean(SE,na.rm=T))),
#        aes(x=Cohort,y=RMSE))+ geom_line()
# 
# ggplot(post_df %>% filter(Year>2006) %>% group_by(ARDY) %>%
#          summarise(RMSE=sqrt(mean(SE,na.rm=T))),
#        aes(x=ARDY,y=RMSE))+ geom_line()


post_df %>% summarise(RMSE=sqrt(mean(SE,na.rm=T)),
                      I90=mean(as.numeric(I90),na.rm=T),
                      I50=mean(I50,na.rm=T))

post_df %>% group_by(V_ind) %>%
  summarise(RMSE=sqrt(mean(SE,na.rm=T)),
            I90=mean(as.numeric(I90),na.rm=T),
            I50=mean(I50,na.rm=T))

#ages <- 15:44
asfr_USA <- data_all %>% filter(Cohort>=start_cohort,
                    ARDY %in% ages,
                    Country=="USA")

asfr_FRATNP <- data_all %>% filter(Cohort>=start_cohort,
                                ARDY %in% ages,
                                Country=="FRATNP")

asfr_SWE <- data_all %>% filter(Cohort>=start_cohort,
                                   ARDY %in% ages,
                                   Country=="SWE")

asfr_USA %>% dim()
asfr_GBRTENW_df %>% dim()
asfr_FRATNP %>% dim()
asfr_SWE %>% dim()

select <- dplyr::select
asfr_USA %>% select(Cohort) %>% range()
asfr_SWE %>% select(Cohort) %>% range()
asfr_FRATNP %>% select(Cohort) %>% range()
asfr_GBRTENW_df %>% select(Cohort) %>% range()

# If we are using data with the same range
post_GBRTENW <- get_posterior(calib$K,
                              data_cntry_df =  asfr_GBRTENW_df,forecast_joy,
                              forecast_horizon,
                              ages=ages)
post_USA <- get_posterior(calib$K,
                         data_cntry_df =  asfr_USA,forecast_joy,
                         forecast_horizon, 
                         ages=ages)

post_SWE <- get_posterior(calib$K,
                          data_cntry_df =  asfr_SWE,forecast_joy,
                          forecast_horizon, ages=ages)

post_FRATNP <- get_posterior(calib$K,
                          data_cntry_df =  asfr_FRATNP,forecast_joy,
                          forecast_horizon,
                          ages=ages)

post_df <- rbind(post_GBRTENW %>% mutate(Country="England and Wales"),
                 post_SWE %>% mutate(Country="Sweden"),
                 post_FRATNP %>% mutate(Country="France"),
                 post_USA %>% mutate(Country="USA"))

post_df %<>%
  mutate(SE=(ASFR_prediction-ASFR)**2,
         I90=((ASFR<p95 &
                 ASFR>p5)),
         I50=((ASFR<p75 &
                 ASFR>p25)))


assess_df <- post_df %>% filter(!V_ind) %>% group_by(Country) %>%
  summarise(RMSE=sqrt(mean(SE,na.rm=T)),
            I90=mean(as.numeric(I90),na.rm=T),
            I50=mean(I50,na.rm=T))

saveRDS(assess_df,file = "processed/assess_schmert.rds")
# plot_ages <- c(18,23,28,33,38,43)
# ggplot(post_df %>% filter(Cohort %in% c(1980,1985,1990,1995,2000)),
#        aes(x=ARDY,y=ASFR_prediction)) +
#   geom_ribbon(aes(ymax=p95,ymin=p5), alpha=0.25) +
#   geom_ribbon(aes(ymax=p75,ymin=p25), alpha=0.5) +
#   geom_line() + facet_grid(Cohort~Country) +
#   geom_point(aes(y=ASFR),colour="red") +
#   geom_vline(aes(xintercept=2006-Cohort))+
#   theme_bw() +
#   xlim(14,49)
# 
# ggplot(post_df %>% filter(ARDY %in% plot_ages),
#        aes(x=Cohort,y=ASFR_prediction)) +
#   geom_ribbon(aes(ymax=p95,ymin=p5), alpha=0.2) +
#   geom_ribbon(aes(ymax=p75,ymin=p25), alpha=0.2) +
#   geom_line() + facet_grid(ARDY~Country, scales="free") +
#   geom_point(aes(y=ASFR),colour="red") +
#   theme_bw()+
#   geom_vline(aes(xintercept=2006-ARDY))
