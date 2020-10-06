#################################################################################################
#                                                                                               #
#  Parameters to load and filter the data                                                       #
#                                                                                               #
#################################################################################################

### Date
#
date_last      <- Sys.Date() - 2         # Data is loaded with one week (7 days) lag
day_zero       <- as.Date("2020-01-01")  # Date corresponding to day / time = 0
# Legacy
# date_cutoff    <- date_last  - 1         # Separates training and test sets (7 days to have a week of testing)
# date_pred      <- Sys.Date() - 2         # Date to predict 

### Data cutoffs
#
confirmed_cases_cutoff         <- 1000   # Filters out countries with less confirmed cases than
fatalities_cutoff              <- 25     # Filters out countries with less fatalities than
density_confirmed_cases_cutoff <- 0.001  # Filters out countries with less density of confirmed cases than
population_cutoff              <- 1e6    # Filters out countries with less population than

### Plotting 
#   Legacy
#   
# plot_flag      <- 0        # 0 : no output. This flag also rules the writting of summaries.
# plot_size      <- 500      # in px
# plot_font_size <- 22       # in px
# plot_quality   <- 75       # jpeg

### Projections and credible intervals
# 
projected_days <- 2 * (as.numeric(Sys.Date() - day_zero )) 
credible_int   <- c(0.025, 0.975)   # 95% credible interval 

### MCMC variables
#
num_chains     <- 3
num_adap       <- 1000

niter_update   <- 1e4      # Recommended 1e3-1e4
niter_sim_init <- 5e3      # Number of MCMC steps for the first iteration
niter_sim_max  <- 1e7      # Exits with (geq). Recommended: 1e7
thinning_init  <- 1        # 100
epsilon_0      <- 2        # Maximum change in the position of the inflection point. Exits with (leq)



#################################################################################################
#                                                                                               #
#  Databases from `covid19datahub/R`                                                            #
#                                                                                               #
#################################################################################################

### Full dataset `full_data`
#   Dates from  `day_zero` to `date_last`

full_data <- covid19(verbose = FALSE, raw = FALSE, level = 1, start = day_zero, end = date_last)

### Selected and fortified dataset `sf_dat`
#   Select columns of interest, rename them, and calculate the number of days elapsed from 
#   day_zero = Date[1].
#   Replaces NAs from the lag by zero which is their consistent value.
#

sf_dat <-full_data %>% 
  select(id                                    , 
         Country = administrative_area_level_1 , 
         Date    = date                        , 
         Confirmed_Cases = confirmed           , 
         Fatalities      = deaths              , 
         Recovered       = recovered           , 
         Population      = population) %>% 
  group_by(Country) %>%
  mutate(Active_Cases            = Confirmed_Cases - Fatalities - Recovered ,
         Density_Active_Cases    = Active_Cases    / Population             ,
         Density_Confirmed_Cases = Confirmed_Cases / Population             , 
         Density_Fatalities      = Fatalities      / Population             ,
         Density_Recovered       = Recovered       / Population             ,
         Daily_Confirmed_Cases   = Confirmed_Cases - lag(Confirmed_Cases)   ,
         Daily_Fatalities        = Fatalities   - lag(Fatalities)           ,
         Daily_Recovered         = Recovered    - lag(Recovered)            , 
         Daily_Active_Cases      = Active_Cases - lag(Active_Cases)         ,
         Time                    = as.numeric(Date - Date[1])) %>% 
  replace_na(list(Daily_Confirmed_Cases = 0, Daily_Fatalities = 0, Daily_Recovered = 0, Daily_Active_Cases = 0)) %>%
  ungroup() 


### Filtered and fortified dataset `latest_global_figures`
#   Summarizes the different cases worldwide.
#
latest_global_figures <- sf_dat %>%
  filter(Date == date_last) %>%
  summarize( total_confirmed_cases = sum(Confirmed_Cases) ,
             total_fatalities      = sum(Fatalities)      ,
             total_recovered       = sum(Recovered)       ,
             total_active_cases    = sum(Confirmed_Cases - Fatalities - Recovered)  ,
             total_population      = sum(Population) ) %>%
  mutate( global_density_confirmed_cases = total_confirmed_cases / total_population ,
          global_density_fatalities      = total_fatalities      / total_population ,
          global_density_recovered       = total_recovered       / total_population ,
          global_density_active_cases    = total_active_cases    / total_population ,
          global_fatality_rate           = total_fatalities / total_confirmed_cases ,
          global_recovery_rate           = total_recovered  / total_confirmed_cases ) 


### Filtered dataset `dat`
#   Includes countries consistent with criterion `Pop+Dens` 
#   See "utilities" file and "data cutoffs" (above)
#

dat <- sf_dat %>% 
  filter(Country %in% include_countries(sf_dat, "Pop+Dens"))



