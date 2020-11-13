#################################################################################################
#                                                                                               #
#  Libraries and remote installs                                                                #
#                                                                                               #
#################################################################################################

remotes::install_github("covid19datahub/R")
packages_needed <- c("globals" ,  "COVID19"  , "dplyr"   , "tidyr"    , "stringr", "rjags"    , 
                     "quantmod",  "ggplot2"  , "ggthemes", "extrafont", "scales" , "lubridate", 
                     "future"  ,  "iterators", "parallel", "rngtools" , "doRNG"  , "foreach"  ,
                     "doFuture" )

packages_to_be_installed <- packages_needed[!(packages_needed %in% installed.packages()[,"Package"])]
if(length(packages_to_be_installed)) install.packages(packages_to_be_installed) 
invisible(lapply(packages_needed, require, character.only = TRUE))


#################################################################################################
#                                                                                               #
#  Future and DoFuture                                                                          #
#                                                                                               #
#################################################################################################

registerDoFuture()
plan(multisession, workers = detectCores() - 1)
options(future.globals.maxSize = 20000 * 1024^2)



#################################################################################################
#                                                                                               #
#  Functions to filter the covid19 database according a given criterion                         #
#                                                                                               #
#################################################################################################

### Function `exclude_countries()`
#
#    Selects countries (to exclude) using a criterion based on a minimum number of
#   "Fatalities", "Confirmed_Cases", or "Pop+Dens" (Density of Population and 
#    Population)
#    Returns a vector of countries names (char)

exclude_countries <- function(df,type){
  if (type == "Fatalities"){
    countries_out <- setdiff( unique(df$Country) , unique(subset(df, Fatalities > fatalities_cutoff)$Country) )
  }
  else if (type == "Confirmed_Cases"){
    countries_out <- setdiff( unique(df$Country) , unique(subset(df, Confirmed_Cases > confirmed_cases_cutoff)$Country) )
  }
  else if (type == "Pop+Dens"){
    countries_out <- setdiff( unique(df$Country) , unique(subset(df, Confirmed_Cases > confirmed_cases_cutoff & Density_Confirmed_Cases > density_confirmed_cases_cutoff & Population > population_cutoff)$Country) )    
  }
  else{
    stop(type,"is not a valid type: Must be Fatalities, Confirmed_Cases, or Pop+Dens")
  }
  invisible(countries_out)
}


### Function `include_countries()`
#
#    Selects countries (to include) using a criterion based on a minimum number of
#   "Fatalities", "Confirmed_Cases", or "Pop+Dens" (Density of Population and 
#    Population)
#    Returns a vector of countries names (char)

include_countries <- function(df, type){
  if (type == "Fatalities"){
    countries_in <-as.character(sort(unique(subset(df, Fatalities > fatalities_cutoff)$Country)))
  }
  else if (type == "Confirmed_Cases"){
    countries_in <-as.character(sort(unique(subset(df, Confirmed_Cases > confirmed_cases_cutoff)$Country)))
  }
  else if (type == "Pop+Dens"){
    countries_in <-as.character(sort(unique(subset(df, Confirmed_Cases > confirmed_cases_cutoff & Density_Confirmed_Cases > density_confirmed_cases_cutoff & Population > population_cutoff)$Country)))
  }
  else{
    stop(type,"is not a valid type: Must be Fatalities, Confirmed_Cases, or Pop+Dens")
  }
  invisible(countries_in)
}


### Function `filter_countries_by_wave_1()`
#
#   Returns (invisible) named list.
#
#   Separates a list (vector) of countries `country_list` present in the dataframe 
#   `df` into two groups: first_wave and second_wave, which are returned in a named list.
#
#   `type` is direct cases, not the daily_cases 
#
#   Uses `findPeaks()` and `findValleys()` from package `quantmod`.

filter_countries_by_wave_1 <- function(df, country_list, type){
  first_wave_countries  <- c()
  second_wave_countries <- c()
  
  for (i in 1:length(country_list)){
    df %>% 
      filter(Country %in% country_list[i]) %>% 
      select(!! rlang::sym(c(type))) %>% 
      rollmean(k = 41, align = "center") %>% 
      smooth.spline( spar = 0.7 ) -> tmp.spline
    
    findPeaks(diff(tmp.spline$y))   -> tmp.peaks
    findValleys(diff(tmp.spline$y)) -> tmp.valleys
    
    criterion <- length(c(tmp.peaks, tmp.valleys ) )
    
    if (criterion >= 3){
      second_wave_countries <- c(second_wave_countries, country_list[i])
    }
    else{
      first_wave_countries <- c(first_wave_countries, country_list[i])
    }
  }
  return_list <- list(first_wave  = first_wave_countries, 
                      second_wave = second_wave_countries)
  invisible(return_list)
}


### Function `filter_countries_by_wave_2()`
#
#   Returns (invisible) named list.
#
#   Separates a list (vector) of countries `country_list` present in the dataframe 
#   `df` into three groups: first_wave, second_wave, and first2second_wave, which
#   are returned in a named list.
#
#   `type` is transformed into daily_cases (see below) as the function uses the 
#   second derivative to classify.
#
#   Uses `findPeaks()` and `findValleys()` from package `quantmod`.

filter_countries_by_wave_2 <- function(df, country_list, type){
  first_wave_countries        <- c()
  second_wave_countries       <- c()
  first2second_wave_countries <- c()
  
  if (str_detect(type,"Daily")==TRUE){
    print("Not valid type. Stopping")
    stop()
  }
  
  type <- paste0("Daily_",type)
  
  for (i in 1:length(country_list)){
    dat %>%
      filter(Country %in% country_list[i]) %>%
      select(!! rlang::sym(c(type))) %>%
      rollmean(k = 41, align = "center") %>%
      smooth.spline( spar = 0.85) -> tmp.spline
    
    findPeaks(abs(diff(tmp.spline$y)))   -> tmp.peaks
    findValleys(abs(diff(tmp.spline$y))) -> tmp.valleys
    
    criterion <- length( c( tmp.peaks[tmp.peaks > 50], tmp.valleys[tmp.valleys > 50] ) )
    
    if (criterion > 3){
      second_wave_countries <- c(second_wave_countries, country_list[i])
    }
    else if (criterion < 3){
      first_wave_countries <- c(first_wave_countries, country_list[i])
    }
    else {
      first2second_wave_countries <- c(first2second_wave_countries, country_list[i])
    }
  }
  
  return_list <- list(first_wave        = first_wave_countries, 
                      second_wave       = second_wave_countries, 
                      first2second_wave = first2second_wave_countries)
  invisible(return_list)  
}


### Function `filter_countries_by_wave_3()`
#
#   Returns (invisible) an ordered vector of 1 (first wave), 2 (second wave) and/or 12 
#   (first to second wave) from list (vector) of countries `country_list`, present in the 
#   dataframe `df`.
#
#   `type` is turned into daily_cases (see below) as the function uses the second 
#   derivative to classify.
#
#   Uses `findPeaks()` and `findValleys()` from package `quantmod`.

filter_countries_by_wave_3 <- function(df, country_list, type){
  waves_countries <- c()
  
  if (str_detect(type,"Daily")==TRUE){
    print("Not valid type. Stopping")
    stop()
  }
  
  type <- paste0("Daily_",type)
  
  for (i in 1:length(country_list)){
    dat %>%
      filter(Country %in% country_list[i]) %>%
      select(!! rlang::sym(c(type))) %>%
      rollmean(k = 41, align = "center") %>%
      smooth.spline( spar = 0.85) -> tmp.spline
    
    findPeaks(abs(diff(tmp.spline$y)))   -> tmp.peaks
    findValleys(abs(diff(tmp.spline$y))) -> tmp.valleys
    
    criterion <- length( c( tmp.peaks[tmp.peaks > 50], tmp.valleys[tmp.valleys > 50] ) )
    
    if (criterion > 3){
      waves_countries <- c(waves_countries, 2) # Second wave
    }
    else if (criterion < 3){
      waves_countries <- c(waves_countries, 1) # First wave
    }
    else {
      waves_countries <- c(waves_countries, 12) # First to second wave
    }
  }
  
  invisible(waves_countries)  
}


### Function `rank_countries_by_wave()`
#
#   Returns (invisible) a dataframe with the countries in `x` classified by 
#   their wave status for the confirmed cases, fatalities and recovered.
#   
#   Uses `filter_countries_by_wave_3()`

rank_countries_by_wave <- function(x, df){
  
  x <- sort(x)
  confirmed_cases_by_wave <- filter_countries_by_wave_3(df, x, "Confirmed_Cases")
  fatalities_by_wave      <- filter_countries_by_wave_3(df, x, "Fatalities")
  recovered_by_wave       <- filter_countries_by_wave_3(df, x, "Recovered")
  
  bind_cols(Country              = x                       , 
            Confirmed_Cases_Wave = confirmed_cases_by_wave , 
            Fatalities_Wave      = fatalities_by_wave      ,
            Recovered_Wave       = recovered_by_wave
  ) -> return_df
  
  invisible(return_df)
}


### Function `get_latest_figures()`
#
#   Returns the last (per `date_last`) cumulative figures from dataframe `df`

get_latest_figures <- function(df){
  df %>%
    group_by(Country) %>%
    filter(Date == date_last) %>%
    ungroup()
}


### Function `get_wave()`
#
#   Returns the wave for `country` and  `type` from `df`.
#   `df` output by `rank_countries_by_wave()`
#   first2second wave -> first wave  : `fswave = 1` (default)
#   first2second wave -> second wave : `fswave = 2`

get_wave <- function(df, country, type, fswave = 1){
  
  tmp_sum <- countries_by_wave %>%
    filter(Country == country) %>%
    select(-Country) %>%
    rowSums(.) 
  
  if (tmp_sum == 14) {fswave <- 1}
  else if(tmp_sum == 15){fswave <- 1}
  else if(tmp_sum == 16){fswave <- 2}
  else if(tmp_sum == 25){fswave <- 1}
  else if(tmp_sum == 26){fswave <- 2}
  else{fswave <- 1}
  
  df %>%
    replace(., .==12, fswave) %>% # first2second wave -> first wave
    filter(Country == country) %>% 
    pull(!! rlang::sym(c(paste0(type,"_Wave")))) 
}



#################################################################################################
#                                                                                               #
# Misc functions. Usually little helpers                                                        #
#                                                                                               #
#################################################################################################

### Function `get_prec()`
#
#   Calculates the precision of a vector x

get_prec <- function(x){ 1/var(x) }


### Function `get_inflection_points()`
#
#   Extracts the position of the peaks of a vector (~ts)
#   * Needs package `quantmod`
#   * Returns a vector 

get_inflection_points <- function(x){findPeaks(diff(x))}



#################################################################################################
#                                                                                               #
# JAGS functions                                                                                #
#                                                                                               #
#################################################################################################

### Function build_model_string()
#
#   Creates the model string for JAGS 
#   Input: normal priors' parameters, wave = 1,2 defaults to 1.
#   Returns : a string containing the model.

build_model_string <- function(type                , 
                               mean_logasymptote_1 , 
                               mean_location_1     , 
                               mean_growth_rate_1  , 
                               prec_logasymptote_1 = 1e-4 , 
                               prec_location_1     = 1e-4 , 
                               prec_growth_rate_1  = 1e-4 , 
                               mean_logasymptote_2 = mean_logasymptote_1 , 
                               mean_location_2     = mean_location_1     , 
                               mean_growth_rate_2  = mean_growth_rate_1  , 
                               prec_logasymptote_2 = 1e-4 , 
                               prec_location_2     = 1e-4 , 
                               prec_growth_rate_2  = 1e-4 , 
                               wave                = 1
                               ){
  
  if (wave == 1){
    model_string = paste0(" model {

    for (i in 1:length(",type,")) {
        ",type,"[i] ~ dnorm(mu[i], prec)
        mu[i] = asymptote_1 * exp( -location_1 * exp( -growth_rate_1 * Time[i]) )
    }
    
        asymptote_1     = exp(logasymptote_1)
        logasymptote_1  ~ dnorm(",mean_logasymptote_1,",", prec_logasymptote_1," )     
        location_1      ~ dnorm(",mean_location_1,    ",", prec_location_1,    " )     
        growth_rate_1   ~ dnorm(",mean_growth_rate_1, ",", prec_growth_rate_1, " )     
        prec            ~ dgamma(1.0, 1.0)
} "
    )
  }
  else if (wave == 2){
    model_string = paste0(" model {

    for (i in 1:length(",type,")) {
        ",type,"[i] ~ dnorm(mu[i], prec)
        mu[i] = asymptote_1 * exp( -location_1 * exp( -growth_rate_1 * Time[i]) ) + asymptote_2 * exp( -location_2 * exp( -growth_rate_2 * Time[i]) )
    }
    
        asymptote_1     = exp(logasymptote_1)
        logasymptote_1  ~ dnorm(",mean_logasymptote_1,",", prec_logasymptote_1," )     
        location_1      ~ dnorm(",mean_location_1,    ",", prec_location_1,    " )     
        growth_rate_1   ~ dnorm(",mean_growth_rate_1, ",", prec_growth_rate_1, " ) T(0,1) 
        
        asymptote_2     = exp(logasymptote_2)
        logasymptote_2  ~ dnorm(",mean_logasymptote_2,",", prec_logasymptote_2," )     
        location_2      ~ dnorm(",mean_location_2,    ",", prec_location_2,    " )     
        growth_rate_2   ~ dnorm(",mean_growth_rate_2, ",", prec_growth_rate_2, " ) T(0,1)       
        
        prec            ~ dgamma(1.0, 1.0)
} "
    ) 
  }
  return(model_string)
}


### Function `run_mcmc_simulation()`
#
#   Runs the MCMC simulation w/JAGS (legacy).
#   * Needs a model_string, keeps  priors unchanged during the `while` loop.
#   * Updates thinning parameter depending on the number of MCMC steps.
#   * Needs wave = 1,2 defaulting to 1.
#   Returns a list with model_sim (MCMC) , its average model_csim, final iter.
#   Updated by `run_mcmc_simulation_2`.

run_mcmc_simulation <- function( index_country    , 
                                 model_string     , 
                                 model_parameters , 
                                 number_chains         = num_chains     , 
                                 number_adap           = num_adap       , 
                                 number_thin_init      = thinning_init  , 
                                 number_steps_update   = niter_update   , 
                                 number_steps_run_init = niter_sim_init , 
                                 number_steps_run_max  = niter_sim_max  , 
                                 epsilon_init          = projected_days , 
                                 epsilon_max           = epsilon_0      ,
                                 wave                  = 1
                                 ){
  
  print(">>>>>>>>>>>>>>>>>>>>>>>>>>>")
  print("")
  print(paste("Country",by_country[index_country]))
  print("")
  print(">>>>>>>>>>>>>>>>>>>>>>>>>>>")
  
  set.seed(1970)
  
  # Subset the country
  jags_datafile <- as.list( subset(dat , Country == by_country[index_country]) )
  
  # Change in the location of the inflection point
  inf_pt_init <- epsilon_init  
  epsilon     <- epsilon_init
  niter_sim   <- number_steps_run_init
  thinning    <- number_thin_init
  
  # Begin while loop to converge the MCMC steps
  while (epsilon >= epsilon_max && niter_sim <= number_steps_run_max ){
    
    # Make sure that for a single wave case we start with three chains
    if (wave == 1){
      number_chains <- 3
    }
    
    model_jags <- jags.model(
                             textConnection(model_string) ,
                             data     = jags_datafile     ,
                             n.chains = number_chains     ,
                             n.adap   = number_adap
                             )
    
    update(
           model_jags , 
           number_steps_update
           ) 
    
    model_sim <- coda.samples(
                              model          = model_jags       ,
                              variable.names = model_parameters ,
                              n.iter         = niter_sim        ,
                              thin           = thinning         ,
                              quiet          = TRUE
                              )
    
    model_csim <- do.call(rbind, model_sim) 
    
    # Calculate the inflection point and distance
    fit_range <- 0:projected_days
    inf_pt    <- get_inflection_points(gompertz_eval(model_csim, fit_range, wave = wave))
    inf_pt    <- mean(inf_pt) 
    epsilon   <- abs( inf_pt - inf_pt_init )
    
    # Monitor 
    print("HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH")
    print(paste("number of MCMC steps",   niter_sim))
    print(paste("Distance to convergence",  epsilon))
    print(paste("Thinning every", thinning,"sample"))
    print(paste("Number of chains", number_chains  ))
    print("HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH")
    
    # Update the position of the inflection point and the number of iterations
    inf_pt_init <- inf_pt
    niter_sim   <- 10*niter_sim
    
    # Update the thinning as function of niter_sim
    if (niter_sim > 1e3 && niter_sim < 1e5){
      if (wave == 1){
        thinning      <- 1 * thinning
        number_chains <- 3
      }
      else if(wave == 2){
        thinning      <- 5 * thinning
        number_chains <- 1
      }
    }
    
    if (niter_sim > 1e4 && niter_sim < 1e6){
      if (wave == 1){
        thinning      <- 1 * thinning
        number_chains <- 3
      }
      else if(wave == 2){
        thinning      <- 5 * thinning
        number_chains <- 1
      }
    }
    
    if (niter_sim > 1e5 && niter_sim < 1e7){
      if (wave == 1){
        thinning      <- 1 * thinning
        number_chains <- 3
      }
      else if(wave == 2){
        thinning      <- 5 * thinning
        number_chains <- 1
      }
    }
    
    if (niter_sim > 1e6 && niter_sim < 1e8){
      if (wave == 1){
        thinning      <- 1 * thinning
        number_chains <- 3
      }
      else if(wave == 2){
        thinning      <- 5 * thinning
        number_chains <- 1
      }
    }
    
  } # end while loop
  return_list<-list(model_sim, model_csim, niter_sim)
  return(return_list)
}


### Function `run_mcmc_simulation_2()`
#
#   Runs the MCMC simulation w/JAGS.
#   * Needs a list with the parameters to create the model_string
#   * Updates 
#     ** priors
#     ** thinning parameter depending on the number of MCMC steps.
#   Returns a list with model_sim (MCMC) , its average model_csim, final iter.

run_mcmc_simulation_2 <- function(index_country                          , 
                                  list_model_string                      , 
                                  number_chains         = num_chains     , 
                                  number_adap           = num_adap       , 
                                  number_thin_init      = thinning_init  , 
                                  number_steps_update   = niter_update   , 
                                  number_steps_run_init = niter_sim_init , 
                                  number_steps_run_max  = niter_sim_max  , 
                                  epsilon_init          = projected_days , 
                                  epsilon_max           = epsilon_0
                                  ){
  
  set.seed(1970)
  
  # Subset the country
  jags_datafile <- as.list( subset(dat , Country == by_country[index_country]) )
  
  # Change in the location of the inflection point
  inf_pt_init <- epsilon_init  
  epsilon     <- epsilon_init
  niter_sim   <- number_steps_run_init
  thinning    <- number_thin_init
  
  # Parameters to build the model string 
  type                <- list_model_string[["type"]]
  mean_logasymptote_1 <- list_model_string[["mean_logasymptote_1"]] 
  mean_location_1     <- list_model_string[["mean_location_1"]] 
  mean_growth_rate_1  <- list_model_string[["mean_growth_rate_1"]] 
  prec_logasymptote_1 <- list_model_string[["prec_logasymptote_1"]] 
  prec_location_1     <- list_model_string[["prec_location_1"]] 
  prec_growth_rate_1  <- list_model_string[["prec_growth_rate_1"]]
  mean_logasymptote_2 <- list_model_string[["mean_logasymptote_2"]]
  mean_location_2     <- list_model_string[["mean_location_2"]] 
  mean_growth_rate_2  <- list_model_string[["mean_growth_rate_2"]] 
  prec_logasymptote_2 <- list_model_string[["prec_logasymptote_2"]] 
  prec_location_2     <- list_model_string[["prec_location_2"]] 
  prec_growth_rate_2  <- list_model_string[["prec_growth_rate_2"]]
  wave                <- list_model_string[["wave"]]  
  
  print(">>>>>>>>>>>>>>>>>>>>>>>>>>>")
  print(">>                         ")
  print(paste(">> Country", by_country[index_country]))
  print(paste(">>> Case"  , type))
  print(paste(">>> Wave"  , wave ))
  print(">>                         ")
  print(">>>>>>>>>>>>>>>>>>>>>>>>>>>")
  
  
  # Begin while loop to converge the MCMC steps
  while (epsilon >= epsilon_max && niter_sim <= number_steps_run_max ){
    
    model_string <- build_model_string(type                = type                , 
                                       
                                       mean_logasymptote_1 = mean_logasymptote_1 ,
                                       mean_location_1     = mean_location_1     ,
                                       mean_growth_rate_1  = mean_growth_rate_1  ,
                                       prec_logasymptote_1 = prec_logasymptote_1 ,
                                       prec_location_1     = prec_location_1     ,
                                       prec_growth_rate_1  = prec_growth_rate_1  ,
                                       
                                       mean_logasymptote_2 = mean_logasymptote_2 ,
                                       mean_location_2     = mean_location_2     ,
                                       mean_growth_rate_2  = mean_growth_rate_2  ,
                                       prec_logasymptote_2 = prec_logasymptote_2 ,
                                       prec_location_2     = prec_location_2     ,
                                       prec_growth_rate_2  = prec_growth_rate_2  ,
                                       
                                       wave                = wave
                                       )
    
    model_jags <- jags.model(
                             textConnection(model_string) ,
                             data     = jags_datafile     ,
                             n.chains = number_chains     ,
                             n.adap   = number_adap       ,
                             inits    = list(.RNG.name="base::Wichmann-Hill", .RNG.seed=1970) ,
                             )
    
    update(
           model_jags, 
           number_steps_update
           ) 
    
    if (wave == 2){
      model_parameters <- c("logasymptote_1" ,
                            "location_1"     ,
                            "growth_rate_1"  , 
                            "logasymptote_2" ,
                            "location_2"     ,
                            "growth_rate_2"  , 
                            "prec"
                             )
    }
    
    if (wave == 1){
      model_parameters <- c("logasymptote_1" ,
                            "location_1"     ,
                            "growth_rate_1"  , 
                            "prec"
                             ) 
    }
    
    model_sim <- coda.samples(
                              model          = model_jags       ,
                              variable.names = model_parameters ,
                              n.iter         = niter_sim        ,
                              thin           = thinning         ,
                              quiet          = TRUE
                              )
    

    model_csim <- do.call(rbind, model_sim) 
    
    # Calculate the inflection point and distance
    fit_range <- 0:projected_days
    inf_pt    <- get_inflection_points(gompertz_eval(model_csim, fit_range, wave = wave))
    inf_pt    <- mean(inf_pt) 
    epsilon   <- abs( inf_pt - inf_pt_init )
    
    # Monitor 
    print("HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH")
    print(paste("number of MCMC steps",  niter_sim ))
    print(paste("Distance to convergence", epsilon ))
    print(paste("Thinning every", thinning,"sample"))
    print(paste("Number of chains", number_chains  ))
    print("HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH")
    
    # Update the position of the inflection point and the number of iterations
    inf_pt_init <- inf_pt
    niter_sim   <- 5 * niter_sim
    
    # Update the thinning as function of niter_sim
    if (niter_sim > 1e3 && niter_sim < 1e5){
      thinning      <- 1 * thinning
      number_chains <- 3
    }
    
    if (niter_sim > 1e4 && niter_sim < 1e6){
      thinning      <- 1 * thinning
      number_chains <- 3
    }
    
    if (niter_sim > 1e5 && niter_sim < 1e7){
      thinning      <- 1 * thinning
      number_chains <- 3
    }
    
    if (niter_sim > 1e6 && niter_sim < 1e8){
      thinning      <- 1 * thinning
      number_chains <- 3
    }
    
    # Update priors
    updated_priors <- as.data.frame(model_csim) %>% 
      select(-prec) %>% 
      summarize_all(list(mean = mean, prec = get_prec)) 
    
    mean_logasymptote_1 <- updated_priors$logasymptote_1_mean 
    mean_location_1     <- updated_priors$location_1_mean 
    mean_growth_rate_1  <- updated_priors$growth_rate_1_mean 
    prec_logasymptote_1 <- 0.5 * updated_priors$logasymptote_1_prec 
    prec_location_1     <- 0.5 * updated_priors$location_1_prec 
    prec_growth_rate_1  <- 0.5 * updated_priors$growth_rate_1_prec 
    
    mean_logasymptote_2 <- updated_priors$logasymptote_2_mean 
    mean_location_2     <- updated_priors$location_2_mean 
    mean_growth_rate_2  <- updated_priors$growth_rate_2_mean 
    prec_logasymptote_2 <- 0.5 * updated_priors$logasymptote_2_prec 
    prec_location_2     <- 0.5 * updated_priors$location_2_prec 
    prec_growth_rate_2  <- 0.5 * updated_priors$growth_rate_2_prec 
    
  } # end while loop
  
  return_list<-list( niter = niter_sim  , 
                     csim  = model_csim , 
                     wave  = wave )
  
  return( return_list )
  
}



#################################################################################################
#                                                                                               #
# Functions to get characteristics of the Gompertz curve                                        #
#                                                                                               #
#################################################################################################

### Function `gompertz_eval()` 
#
#   Evaluates the Gompertz function from the average posterior values collected in `model_csim`, for a natural time range 
#   in days collected in vector `date_eval` starting from `date_zero` which  sets the variable `Time=0`. 
#   Returns a vector to be interpreted as ts (but not a ts object).

gompertz_eval <- function(model_csim, date_eval, wave = 1){
  if (class(date_eval) == "Date"){
    date_eval <- as.numeric(as.Date(date_eval)-as.Date(day_zero), "days") 
  }
  
  geval <- exp(colMeans(model_csim)[["logasymptote_1"]]) * 
      exp( -colMeans(model_csim)[["location_1"]] * exp(-as.numeric(colMeans(model_csim)[["growth_rate_1"]], units="days") * date_eval) )

  
  if (wave == 2){
    geval <- geval + exp(colMeans(model_csim)[["logasymptote_2"]]) * 
      exp( -colMeans(model_csim)[["location_2"]] * exp(-as.numeric(colMeans(model_csim)[["growth_rate_2"]] ,units="days") * date_eval) )
  }
  
  return(geval)
}


### Function `comp_deriv_gompertz()` 
#
#   Evaluates the components to the derivative of the Gompertz function from the average posterior values 
#   collected in `model_csim`, for numeric values collected in vector `x` (corresponding to the variable `Time`) 
#   Returns a list with two vectors one for each component to be interpreted as ts (but not a ts object).
#   Only for wave 2

comp_deriv_gompertz <- function(model_csim, x){
  if (class(x) == "Date"){
    x <- as.numeric(as.Date(x)-as.Date(day_zero), "days") 
  }
  
  comp_1 <- exp(colMeans(model_csim)[["logasymptote_1"]]) * colMeans(model_csim)[["location_1"]] * colMeans(model_csim)[["growth_rate_1"]] /
    exp( colMeans(model_csim)[["growth_rate_1"]] * x + colMeans(model_csim)[["location_1"]] * exp(-colMeans(model_csim)[["growth_rate_1"]] * x ) )
  
  comp_2 <-  exp(colMeans(model_csim)[["logasymptote_2"]]) * colMeans(model_csim)[["location_2"]] * colMeans(model_csim)[["growth_rate_2"]] /
    exp( colMeans(model_csim)[["growth_rate_2"]] * x + colMeans(model_csim)[["location_2"]] * exp(-colMeans(model_csim)[["growth_rate_2"]] * x ) )
  
  return_list <- list( comp_1 = comp_1 , 
                       comp_2 = comp_2  )
  
  return(return_list)
}



#################################################################################################
#                                                                                               #
# Functions to get credible intervals from MCMC posteriors                                      #
#                                                                                               #
#################################################################################################

### Function `model_credible_interval()`
#
#   Returns the credible interval for a Gompertz function from the posterior values of its 
#   parameters `model_csim` for a single time value `x` in days units from `date_zero`.

model_credible_interval <- function(model_csim, x, wave){
  if (class(x) == "Date"){
    x <- as.numeric(as.Date(x)-as.Date(day_zero), "days") 
  } 
  
  mci <- quantile(
    exp(  model_csim[,"logasymptote_1"]) * 
      exp( -model_csim[,"location_1"] * exp(-model_csim[,"growth_rate_1"] * x) ) ,
    probs=credible_int, names = FALSE  )
  
  if (wave == 2){
    mci <- mci + quantile(
      exp(  model_csim[,"logasymptote_2"]) * 
        exp( -model_csim[,"location_2"] * exp(-model_csim[,"growth_rate_2"] * x) ) ,
      probs=credible_int, names = FALSE  )
  }
  
  return(mci)
}


### Function `deriv_model_credible_interval()`
#
#   Returns the credible interval for the first derivative of a Gompertz function from the posterior 
#   values of its parameters `model_csim` for a single time value `x` in days units from `date_zero`.

deriv_model_credible_interval <- function(model_csim, x, wave){
  if (class(x) == "Date"){
    x <- as.numeric(as.Date(x)-as.Date(day_zero), "days") 
  } 
  dmci <- quantile(
    exp( model_csim[,"logasymptote_1"] )  * model_csim[,"location_1"] * model_csim[,"growth_rate_1"] / 
      exp( model_csim[,"growth_rate_1"] * x + model_csim[,"location_1"] * exp(-model_csim[,"growth_rate_1"] * x ) ) ,
    probs = credible_int, names = FALSE )
  
  if (wave == 2){
    dmci <- dmci + quantile(
      exp( model_csim[,"logasymptote_2"] )  * model_csim[,"location_2"] * model_csim[,"growth_rate_2"] / 
        exp( model_csim[,"growth_rate_2"] * x + model_csim[,"location_2"] * exp(-model_csim[,"growth_rate_2"] * x ) ) ,
      probs = credible_int, names = FALSE )
  }
  
  return(dmci)
}


### Function `get_credible_intervals()`
#
#   Returns a matrix `dim(2, length(x))` with the credible interval for  Gompertz function 
#   and its derivative from the posterior values of its parameters `model_csim` for a range 
#   `x` in days units from `date_zero`.
#
#   Depends on functions `model_credible_interval`  and `deriv_model_credible_interval`
#

get_credible_intervals <- function(model_csim, x, wave){
  # ci       <- sapply(x, model_credible_interval,       model_csim = model_csim, wave = wave )
  # ci_deriv <- sapply(x, deriv_model_credible_interval, model_csim = model_csim, wave = wave )
  
  ci       <- vapply(x, model_credible_interval,       model_csim = model_csim, wave = wave, FUN.VALUE = double(2) )
  ci_deriv <- vapply(x, deriv_model_credible_interval, model_csim = model_csim, wave = wave, FUN.VALUE = double(2) )
  
  return_list <- list( model_lower_ci   = ci[1,]       , 
                       model_upper_ci   = ci[2,]       , 
                       d_model_lower_ci = ci_deriv[1,] , 
                       d_model_upper_ci = ci_deriv[2,]  )
  
  return(return_list)
}



#################################################################################################
#                                                                                               #
# Functions to visualize the COVID19 related data                                               #
#                                                                                               #
#################################################################################################

### Function `top_10_ranking()`
#
#   For a dataframe `df` and `type` outputs a ggplot with the top 10 and bottom 10
#   ranking of countries present in `df`.
#
#   Uses `stringr` for character replacement and pattern detection.

top_10_ranking <- function(df, type){
  
  ltype <- str_replace_all(str_to_lower(type), "_", " ")
  
  if (str_detect(ltype, "density") == TRUE){
    ltype %>%
      str_replace_all(.,  "density", " ") %>%
      str_trim() -> ltype 
    
    my_factor <- 1e0
    my_text   <- paste(ltype, "(per capita)")
    }
  else if(str_detect(ltype, "daily") == TRUE){
    my_factor <- 1e3
    my_text   <- paste(ltype, "(in thousands)")
    }
  else {
    my_factor <- 1e6
    my_text   <- paste(ltype, "(in millions)")
    }
  
  my_factor <- as.character(my_factor)
  
  rank <- data.frame("Rank" = c(rep("Low", 10), rep("High",10)))
  
  
  df %>% 
    group_by(Country) %>% 
    arrange(!! rlang::sym(c(type))) %>%
    ungroup() %>%
    slice(c(1:10, (n()-9):n() )) %>%
    select(!!! rlang::syms(c("Country", type))) %>%
    bind_cols(.,rank) -> df_tmp
  
  df_tmp %>%
    ggplot(aes_string(x = paste0(type,"/", my_factor), y = "Country", color = "Rank" )) + 
    geom_point(size = 4) +
    geom_segment(aes(xend = 0, yend = Country), size = 2) +
    scale_color_manual(values=c("#E41A1C", "#377EB8")) + # Set1 
    scale_x_continuous("", position = "top") + # limits = c(0, 7), expand = c(0,0.1)
    labs(title   = paste("Highest and lowest number of", my_text), 
         caption = paste("Source: A. Diaz Ortiz, Data from COVID-19 Data Hub as of",as.character(date_last))) +
    theme_tufte(base_family = "sans") +
    theme( 
      axis.line.y  = element_blank() ,
      axis.ticks.y = element_blank() ,
      axis.title   = element_blank() ,
      axis.text    = element_text(color = "black"),
      legend.position = "none"
    ) +
    scale_y_discrete(limits = df_tmp$Country) 
}


### Function `plot_cases_2()`
#
#   Plots the data, the gompertz fit, together with credible intervals using ggplot
#    
#   Input variables are :
#     * type = char with the variable name, i.e., "Daily_Confirmed_Cases".
#     * df = dataframe, i.e., dat
#     * x = range to fit. numeric
#     * model_csim = mcmc posteriors
#     * wave = 1 o 2 depending on country

plot_cases_2 <- function(country_index, type, df, x, model_csim, wave){
  
  min_date   <- min(as.Date(x, origin = day_zero)) 
  max_date   <- max(as.Date(x, origin = day_zero)) 
  date_x_lim <- as.character(c(min_date, max_date + 150))
  
  country_pop <- df %>%
    filter(Country == by_country[country_index]) %>%
    select(Population) %>%
    slice_tail(n = 1) 
  
  ylabel    <- type %>% 
    str_replace_all(.,  "_", " ") %>%
    str_trim() %>%
    str_to_sentence(., locale = "en") 
  
  # Select the y labels and data helpers
  
  if (type == "Active_Cases"){ 
    
    gfit  <- model_csim
    if (gfit[length(gfit)] < 0){
      label_max_gfit <-  "lim. 0+"
    }
    else{
    label_max_gfit <- as.character( prettyNum(round(gfit[length(gfit)], 0), big.mark = " ") )
    }
    
    inflexion_pts <- get_inflection_points(gfit)
    inflexion_pts <- as.Date(round(inflexion_pts, 0), origin = day_zero)
    
    data_help <- bind_cols( Time = x                             , 
                            Date = as.Date(x, origin = day_zero) , 
                            type = gfit                      
                            )
    names(data_help) <- c("Time", "Date",  type)
  }
  else{  
    gfit    <- gompertz_eval(model_csim, x, wave = wave)
    gfit_c1 <- gompertz_eval(model_csim, x, wave = 1)
    gfit_c2 <- gfit - gfit_c1
    
    label_max_gfit    <- as.character( prettyNum(round(gfit[length(gfit)],       0), big.mark = " ") )
    label_max_gfit_c1 <- as.character( prettyNum(round(gfit_c1[length(gfit_c1)], 0), big.mark = " ") )
    label_max_gfit_c2 <- as.character( prettyNum(round(gfit_c2[length(gfit_c2)], 0), big.mark = " ") )
    
    inflexion_pts <- get_inflection_points(gfit)
    inflexion_pts <- as.Date(round(inflexion_pts, 0), origin = day_zero)
    
    ci_ribbon <- get_credible_intervals(model_csim, x, wave = wave)
    
    data_help <- bind_cols( Time        = x                             , 
                            Date        = as.Date(x, origin = day_zero) , 
                            Lower_ci    = ci_ribbon[["model_lower_ci"]] , 
                            Upper_ci    = ci_ribbon[["model_upper_ci"]] ,
                            type        = gfit                          ,
                            Component_1 = gfit_c1                       ,
                            Component_2 = gfit_c2
                            )
    
    names(data_help) <- c("Time", "Date", "Lower_ci", "Upper_ci", type, "Component_1", "Component_2")
    }
  
  if (wave == 2){
    
    df %>%
      filter(Country == by_country[country_index]) %>%
      ggplot(aes_string(x = "Date", y = type)
      ) +
      geom_point( shape = 22          , 
                  size  = 3           , 
                  color = "#440154FF" , 
                  fill  = "#440154FF" ,  
                  alpha = 0.35
      ) +
      geom_line( data     = data_help , 
                 linetype = "solid"   , 
                 color    = "#E41A1C" , 
                 size     = 0.75
      ) +
      geom_ribbon( data  = data_help                                 , 
                   aes(x = Date, ymin =  Lower_ci, ymax =  Upper_ci) , 
                   fill  = "#E41A1C"                                     , 
                   alpha = 0.2
      ) +
      geom_line( data     = data_help           , 
                 aes(x = Date, y = Component_1) ,
                 linetype = "solid"             , 
                 color    = "#377EB8"           , 
                 size     = 0.75                ,
                 alpha    = 0.70
      ) +  
      geom_line( data     = data_help           , 
                 aes(x = Date, y = Component_2) ,
                 linetype = "solid"             , 
                 color    = "#377EB8"           , 
                 size     = 0.75                ,
                 alpha    = 0.70
      ) +
      labs( title   = by_country[country_index] ,  
            y       = ylabel                    , 
            caption = paste("Source: A. Diaz Ortiz, Data from COVID-19 Data Hub as of",as.character(date_last))
      ) +
      theme_tufte( base_family = "sans" , 
                   base_size   = 13
      )  +
      annotate( "text"                           ,
                label = label_max_gfit_c1        ,
                x     = max_date + 60            ,
                y     = gfit_c1[length(gfit_c1)] ,
                vjust = "top"                    ,
                size  = 3                        ,
                color = "grey40"                 ,
                angle = 0
      )  +
      annotate( "text"                           ,
                label = label_max_gfit_c2        ,
                x     = max_date + 60            ,
                y     = gfit_c2[length(gfit_c2)] ,
                vjust = "top"                    ,
                size  = 3                        ,
                color = "grey40"                 ,
                angle = 0
      ) +
      annotate( "text"                     ,
                label = label_max_gfit     ,
                x     = max_date + 60      ,
                y     = gfit[length(gfit)] ,
                vjust = "bottom"           ,
                size  = 3                  ,
                color = "grey40"           ,
                angle = 0
      ) + 
      scale_y_continuous(labels = label_number()) 
  }
  else if (wave == 1){
    
    if( type == "Active_Cases" ){
      df %>%
        filter(Country == by_country[country_index]) %>%
        ggplot(aes_string(x = "Date", y = type)) +
        geom_point( shape = 22          , 
                    size  = 3           , 
                    color = "#440154FF" , 
                    fill  = "#440154FF" ,  
                    alpha = 0.35
        ) +
        geom_line( data     = data_help , 
                   linetype = "solid"   , 
                   color    = "#E41A1C" , 
                   size     = 0.75
        ) +
        labs( title   = by_country[country_index] ,  
              y       = ylabel                    , 
              caption = paste("Source: A. Diaz Ortiz, Data from COVID-19 Data Hub as of",as.character(date_last))
        ) +
        theme_tufte( base_family = "sans" , 
                     base_size   = 13
        )  +
        annotate( "text"                     ,
                  label = label_max_gfit     ,
                  x     = max_date + 60      ,
                  y     = gfit[length(gfit)] ,
                  vjust = "bottom"           ,
                  size  = 3                  ,
                  color = "grey40"           ,
                  angle = 0
        ) +
        scale_y_continuous(labels = label_number())
    }
    else{
    df %>%
      filter(Country == by_country[country_index]) %>%
      ggplot(aes_string(x = "Date", y = type)) +
      geom_point( shape = 22          , 
                  size  = 3           , 
                  color = "#440154FF" , 
                  fill  = "#440154FF" ,  
                  alpha = 0.35
      ) +
      geom_line( data     = data_help , 
                 linetype = "solid"   , 
                 color    = "#E41A1C" , 
                 size     = 0.75
      ) +
      geom_ribbon( data  = data_help                                 , 
                   aes(x = Date, ymin =  Lower_ci, ymax =  Upper_ci) , 
                   fill  = "#E41A1C"                                     , 
                   alpha = 0.2
      ) +
      labs( title   = by_country[country_index] ,  
            y       = ylabel                    , 
            caption = paste("Source: A. Diaz Ortiz, Data from COVID-19 Data Hub as of",as.character(date_last))
      ) +
      theme_tufte( base_family = "sans" , 
                   base_size   = 13
      )  +
      annotate( "text"                     ,
                label = label_max_gfit     ,
                x     = max_date + 60      ,
                y     = gfit[length(gfit)] ,
                vjust = "bottom"           ,
                size  = 3                  ,
                color = "grey40"           ,
                angle = 0
      ) + 
      scale_y_continuous(labels = label_number())
    }
  }
}


### Function `plot_daily_cases_2()`
#
#   Plots the data, the gompertz fit, together with credible intervals using ggplot
#    
#   Input variables are :
#     * type       = char with the variable name, i.e., "Daily_Confirmed_Cases".
#     * df         = dataframe, i.e., dat
#     * x          = range to fit. numeric
#     * model_csim = mcmc posteriors
#     * wave       = 1 o 2 depending on country

plot_daily_cases_2 <- function(country_index, type, df, x, model_csim, wave){
  
  min_date   <- min(as.Date(x, origin = day_zero)) 
  max_date   <- max(as.Date(x, origin = day_zero)) 
  date_x_lim <- as.character(c(min_date, max_date + 150))
  
  dcases <- df %>%
    filter(Country == by_country[country_index]) %>%
    slice_max(n = 1, order_by = !! rlang::sym(c(type))) %>%
    select(!! rlang::sym(c(type))) 
  
  max_dcases <- dcases[[type]]
  
  ylabel    <- type %>% 
    str_replace_all(.,  "_", " ") %>%
    str_trim() %>%
    str_to_sentence(., locale = "en") 
  
  # Select the y labels and data helpers
  if (type == "Daily_Active_Cases"){ 
    
    gfit          <- model_csim
    dgfit         <- diff(gfit)
    max_dgfit     <- max(diff(gfit))
    
    inflexion_pts <- get_inflection_points(gfit)
    inflexion_pts <- as.Date(round(inflexion_pts, 0), origin = day_zero)
    if (length(inflexion_pts) > 1){wave <- 2}
    
    y_pos     <- max(max_dgfit, max_dcases)
    
    data_help <- bind_cols( Time     = x[-1]                               , 
                            Date     = as.Date(x[-1], origin = day_zero)   , 
                            type = dgfit  
                            )
    
    names(data_help) <- c("Time", "Date", type)
  }
  else {
    gfit          <- gompertz_eval(model_csim, x, wave = wave)
    dgfit         <- diff(gfit)
    max_dgfit     <- max(diff(gfit))
    
    inflexion_pts <- get_inflection_points(gfit)
    inflexion_pts <- as.Date(round(inflexion_pts, 0), origin = day_zero)
    
    ci_ribbon     <- get_credible_intervals(model_csim, x, wave = wave)
    max_ci_ribbon <- max(ci_ribbon[["d_model_upper_ci"]][-1])
    
    y_pos      <- max(max_ci_ribbon, max_dgfit, max_dcases)

    data_help <- bind_cols( Time     = x[-1]                               , 
                            Date     = as.Date(x[-1], origin = day_zero)   , 
                            Lower_ci = ci_ribbon[["d_model_lower_ci"]][-1] , 
                            Upper_ci = ci_ribbon[["d_model_upper_ci"]][-1] ,
                            type = dgfit  
                            )         
    
    names(data_help) <- c("Time", "Date", "Lower_ci", "Upper_ci", type)
  
  }

  if (wave == 2){
    if ( type == "Daily_Active_Cases" ){
      df %>%
        filter(Country == by_country[country_index]) %>%
        ggplot(aes_string(x = "Date", y = type)) +
        geom_point( shape = 22          , 
                    size  = 3           , 
                    color = "#440154FF" , 
                    fill  = "#440154FF" ,  
                    alpha = 0.35
        ) +
        geom_line( data     = data_help , 
                   linetype = "solid"   , 
                   color    = "#E41A1C" , 
                   size     = 0.75
        ) +
        labs( title   = by_country[country_index] ,  
              y       = ylabel                    , 
              caption = paste("Source: A. Diaz Ortiz, Data from COVID-19 Data Hub as of",as.character(date_last))
        ) +
        theme_tufte( base_family = "sans" , 
                     base_size   = 13
        )  +
        geom_segment(aes(x    = inflexion_pts[1] , 
                         xend = inflexion_pts[1] ,
                         y    = 1.15 * y_pos     , 
                         yend = 1.05 * y_pos
        ),
        color    = "grey40" , 
        linetype = "dotted"
        ) +
        geom_segment(aes(x    = inflexion_pts[length(inflexion_pts)] , 
                         xend = inflexion_pts[length(inflexion_pts)] ,
                         y    = 1.15 * y_pos     , 
                         yend = 1.05 * y_pos
        ),
        color    = "grey40" , 
        linetype = "dotted"
        ) +
        annotate( "text"                                 , 
                  label = as.character(inflexion_pts[1]) ,
                  x     = inflexion_pts[1]               , 
                  y     = 1.2 * y_pos                    ,
                  vjust = "bottom"                       , 
                  size  = 3                              , 
                  color = "grey40"                       , 
                  angle = 0
        ) +
        annotate( "text"                                 , 
                  label = as.character(inflexion_pts[length(inflexion_pts)]) ,
                  x     = inflexion_pts[length(inflexion_pts)]               , 
                  y     = 1.2 * y_pos                    ,
                  vjust = "bottom"                       ,  
                  size  = 3                              , 
                  color = "grey40"                       , 
                  angle = 0
        ) +
        scale_y_continuous(labels = label_number())
    }
    else{
    df %>%
      filter(Country == by_country[country_index]) %>%
      ggplot(aes_string(x = "Date", y = type)) +
      geom_point( shape = 22          , 
                  size  = 3           , 
                  color = "#440154FF" , 
                  fill  = "#440154FF" ,  
                  alpha = 0.35
      ) +
      geom_line( data     = data_help , 
                 linetype = "solid"   , 
                 color    = "#E41A1C" , 
                 size     = 0.75
      ) +
      geom_ribbon( data  = data_help                                 , 
                   aes(x = Date, ymin =  Lower_ci, ymax =  Upper_ci) , 
                   fill  = "#E41A1C"                                     , 
                   alpha = 0.2
      ) +
      labs( title   = by_country[country_index] ,  
            y       = ylabel                    , 
            caption = paste("Source: A. Diaz Ortiz, Data from COVID-19 Data Hub as of",as.character(date_last))
      ) +
      theme_tufte( base_family = "sans" , 
                   base_size   = 13
      )  +
      geom_segment(aes(x    = inflexion_pts[1] , 
                       xend = inflexion_pts[1] ,
                       y    = 1.15 * y_pos     , 
                       yend = 1.05 * y_pos
                       ),
                   color    = "grey40" , 
                   linetype = "dotted"
      ) +
      geom_segment(aes(x    = inflexion_pts[length(inflexion_pts)] , 
                       xend = inflexion_pts[length(inflexion_pts)] ,
                       y    = 1.15 * y_pos     , 
                       yend = 1.05 * y_pos
                       ),
                     color    = "grey40" , 
                     linetype = "dotted"
      ) +
      annotate( "text"                                 , 
                label = as.character(inflexion_pts[1]) ,
                x     = inflexion_pts[1]               , 
                y     = 1.2 * y_pos                    ,
                vjust = "bottom"                       , 
                size  = 3                              , 
                color = "grey40"                       , 
                angle  = 0
      ) +
      annotate( "text"                                 , 
                label = as.character(inflexion_pts[length(inflexion_pts)]) ,
                x     = inflexion_pts[length(inflexion_pts)]               , 
                y     = 1.2 * y_pos                    ,
                vjust = "bottom"                       ,  
                size = 3                               , 
                color = "grey40"                       , 
                angle = 0
      ) +
      scale_y_continuous(labels = label_number())
    }
  }
  else if (wave == 1){
    if ( type == "Daily_Active_Cases" ){
      df %>%
        filter(Country == by_country[country_index]) %>%
        ggplot(aes_string(x = "Date", y = type)) +
        geom_point( shape = 22          , 
                    size  = 3           , 
                    color = "#440154FF" , 
                    fill  = "#440154FF" ,  
                    alpha = 0.35
        ) +
        geom_line( data     = data_help , 
                   linetype = "solid"   , 
                   color    = "#E41A1C" , 
                   size     = 0.75
        ) +
        labs( title   = by_country[country_index] ,  
              y       = ylabel                    , 
              caption = paste("Source: A. Diaz Ortiz, Data from COVID-19 Data Hub as of",as.character(date_last))
        ) +
        theme_tufte( base_family = "sans" , 
                     base_size   = 13
        )  +
        geom_segment(aes(x    = inflexion_pts[1] , 
                         xend = inflexion_pts[1] ,
                         y    = 1.15 * y_pos     , 
                         yend = 1.05 * y_pos
        ),
        color    = "grey40" , 
        linetype = "dotted"
        ) +
        annotate( "text"                                 , 
                  label = as.character(inflexion_pts[1]) ,
                  x     = inflexion_pts[1]               , 
                  y     = 1.2 * y_pos                     ,
                  vjust = "middle"                       , 
                  size  = 3                              , 
                  color = "grey40"                       , 
                  angle = 0
        ) +
        scale_y_continuous(labels = label_number())
    }
    else{
    df %>%
      filter(Country == by_country[country_index]) %>%
      ggplot(aes_string(x = "Date", y = type)) +
      geom_point( shape = 22          , 
                  size  = 3           , 
                  color = "#440154FF" , 
                  fill  = "#440154FF" ,  
                  alpha = 0.35
      ) +
      geom_line( data     = data_help , 
                 linetype = "solid"   , 
                 color    = "#E41A1C" , 
                 size     = 0.75
      ) +
      geom_ribbon( data  = data_help                                 , 
                   aes(x = Date, ymin =  Lower_ci, ymax =  Upper_ci) , 
                   fill  = "#E41A1C"                                     , 
                   alpha = 0.2
      ) +
      labs( title   = by_country[country_index] ,  
            y       = ylabel                    , 
            caption = paste("Source: A. Diaz Ortiz, Data from COVID-19 Data Hub as of",as.character(date_last))
      ) +
      theme_tufte( base_family = "sans" , 
                   base_size   = 13
      )  +
      geom_segment(aes(x    = inflexion_pts[1] , 
                       xend = inflexion_pts[1] ,
                       y    = 1.15 * y_pos     , 
                       yend = 1.05 * y_pos
                       ),
                   color    = "grey40" , 
                   linetype = "dotted"
      ) +
      annotate( "text"                                 , 
                label = as.character(inflexion_pts[1]) ,
                x     = inflexion_pts[1]               , 
                y     = 1.2 * y_pos                     ,
                vjust = "middle"                       , 
                size  = 3                              , 
                color = "grey40"                       , 
                angle = 0
      ) +
      scale_y_continuous(labels = label_number())
    } 
  } 
} 


### Function `visualize_cases_mcmc()`
#
#   Outputs a list of plots from the MCMC data 
#    
#   Input variables are :
#     * country_index    = index of by_country vector
#     * type             = char with the variable name, i.e., "Daily_Confirmed_Cases"
#     * df               = dataframe, i.e., dat
#     * x                = range to fit. numeric
#     * list_output_mcmc = list with output of `run_mcmc_simulations_2()`
#     * wave             = 1 o 2 depending on country

visualize_cases_mcmc <- function(country_index, type, df, x, list_output_mcmc, wave = 1){
  
  dtype      <- paste0("Daily_",type)
  model_csim <- list_output_mcmc[["csim"]]
  wave       <- list_output_mcmc[["wave"]]
  
  min_date   <- min(as.Date(x, origin = day_zero)) 
  max_date   <- max(as.Date(x, origin = day_zero)) 
  date_x_lim <- as.character(c(min_date, max_date + 150))
  
  case_y_lim <- df %>%
    filter(Country == by_country[country_index]) %>%
    select(!! rlang::sym(c(type))) %>%
    range() 
  
  dcase_y_lim <- df %>%
    filter(Country == by_country[country_index]) %>%
    select(!! rlang::sym(c(dtype))) %>%
    unname() %>%
    unlist() %>%
    quantile(probs = c(0, 0.99))
  
  case_x_lim <- df %>%
    select(Time) %>%
    range() %>%
    lubridate::as_date(origin = day_zero)
  
  pc  <- plot_cases_2(country_index, type, df, x, model_csim, wave)          + 
            scale_x_date(labels = date_format("%Y-%m"))                      +
            coord_cartesian(xlim = lubridate::as_date(date_x_lim))
  pdc <- plot_daily_cases_2(country_index, dtype, df, x, model_csim, wave)   + 
            scale_x_date(labels = date_format("%Y-%m"))                      +
            coord_cartesian(xlim = lubridate::as_date(date_x_lim))
  
  pcz  <- plot_cases_2(country_index, type, df, x, model_csim, wave)         + 
            scale_x_date(labels = date_format("%Y-%m"))                      +
            coord_cartesian(xlim = case_x_lim, ylim =  case_y_lim )
  
  pdcz <- plot_daily_cases_2(country_index, dtype, df, x, model_csim, wave)  + 
            scale_x_date(labels = date_format("%Y-%m"))                      +
            coord_cartesian(xlim = case_x_lim, ylim = dcase_y_lim )
  
  return_list<-list(Case = pc, Daily_Case = pdc, Case_Zoom = pcz, Daily_Case_Zoom = pdcz )
  return(return_list)
}


### Function `screen_top_10()`
#
#   Returns screen plots for all possible cases in `mycases` plus population
#   produced by function `top_10_ranking()` and dataframe `df`

screen_top_10 <- function(df){
  mycases <- c("Confirmed_Cases", "Fatalities", "Recovered", "Active_Cases")
  
  for (i in mycases){
    print(top_10_ranking(df = df , type = i))
    print(top_10_ranking(df = df , type = paste0("Daily_",i)))
    print(top_10_ranking(df = df , type = paste0("Density_",i)))
  }
  
  print(top_10_ranking(df = latest_figures , type = "Population"))
}


### Function `store_top_10()`
#   
#   Writes to disk `top_10_ranking()` plots for all cases in `mycases` + `Population`.
#   * save_path = path to store on disk
#   * res = resolution in dpi (see ggsave) defaults to "print".

store_top_10 <- function( df, save_path, res= "print" ){
  mycases <- c("Confirmed_Cases", "Fatalities", "Recovered", "Active_Cases")
  
  for (type in mycases){
    
    cpath_tmp <- file.path( save_path , paste0( "top10_", type, "_", date_last , ".tiff") )
    ggsave( file   = cpath_tmp , 
            plot   = top_10_ranking(df = df , type = type) ,
            width  = 7.29 , height = 4.51 , unit   = "in" , dpi = res 
    )
    
    cpath_tmp <- file.path( save_path , paste0( "top10_", "Daily_", type, "_", date_last , ".tiff") )
    ggsave( file  = cpath_tmp ,
            plot  = top_10_ranking(df = df , type = paste0("Daily_", type)) ,
            width = 7.29 , height = 4.51 , unit ="in", dpi = res 
    )
    
    cpath_tmp  <- file.path( save_path , paste0( "top10_", "Density_", type, "_", date_last , ".tiff") )
    ggsave( file  = cpath_tmp , 
            plot  = top_10_ranking(df = df , type = paste0("Density_",type)) ,
            width = 7.29, height = 4.51, unit   ="in", dpi = res )
  }
  
  type <- "Population"
  cpath_tmp <- file.path( save_path , paste0( "top10_", type, "_", date_last , ".tiff") )
  ggsave( file   = cpath_tmp , 
          plot   = top_10_ranking(df = df , type = type) ,
          width  = 7.29 , height = 4.51 , unit   = "in" , dpi = res 
  )
}


### Function `screen_plots()`
#
#   Returns plots for all unique combinations of `by_cases` and `by_country`
#   stored in list `list_vis` (output from `visualize_cases_mcmc()`)

screen_plots <- function(list_vis){
  
  invisible(mapply(function(type, country) print(list_vis[[type]][[country]][["Case"]] ), 
                   expand.grid( by_cases, by_country, stringsAsFactors = FALSE )$Var1, 
                   expand.grid( by_cases, by_country, stringsAsFactors = FALSE )$Var2))
  
  invisible(mapply(function(type, country) print(list_vis[[type]][[country]][["Daily_Case"]] ) , 
                   expand.grid( by_cases, by_country, stringsAsFactors = FALSE )$Var1, 
                   expand.grid( by_cases, by_country, stringsAsFactors = FALSE )$Var2))
  
  invisible(mapply(function(type, country) print(list_vis[[type]][[country]][["Case_Zoom"]] ) , 
                   expand.grid( by_cases, by_country, stringsAsFactors = FALSE )$Var1, 
                   expand.grid( by_cases, by_country, stringsAsFactors = FALSE )$Var2))
  
  invisible(mapply(function(type, country) print(list_vis[[type]][[country]][["Daily_Case_Zoom"]] ) , 
                   expand.grid( by_cases, by_country, stringsAsFactors = FALSE )$Var1, 
                   expand.grid( by_cases, by_country, stringsAsFactors = FALSE )$Var2))
}


### Function `store_plots()`
#
#   Writes to disk plots for all unique combinations of `by_cases` and `by_country`
#   stored in list `list_vis` (output from `visualize_cases_mcmc()`)
#   * save_path = path to store on disk
#   * res = resolution in dpi (see ggsave) defaults to "print".

store_plots <- function(list_vis, save_path, res = "print"){
  for (country in by_country){
    for (type in by_cases){
      
      dpath_tmp  <- file.path( save_path , 
                               paste0( country , "_Daily_" , type , "_", date_last , ".tiff") 
      )
      
      ggsave( list_vis[[type]][[country]][["Daily_Case"]] , 
              file   = dpath_tmp , 
              width  = 7.29      , 
              height = 4.51      , 
              unit   = "in"      , 
              dpi    = res
      )
      
      cpath_tmp  <- file.path( save_path , 
                               paste0( country , "_" , type , "_" , date_last , ".tiff") 
      )
      
      ggsave( list_vis[[type]][[country]][["Case"]] , 
              file   = cpath_tmp, 
              width  = 7.29     , 
              height = 4.51     , 
              unit   ="in"      , 
              dpi    = res
      )
      
      czpath_tmp  <- file.path( save_path , 
                                paste0( country , "_" , type , "_Zoom_" , date_last , ".tiff") 
      )
      
      ggsave( list_vis[[type]][[country]][["Case_Zoom"]] , 
              file   = czpath_tmp, 
              width  = 7.29     , 
              height = 4.51     , 
              unit   ="in"      , 
              dpi    = res
      )
      
      dzpath_tmp  <- file.path( save_path , 
                                paste0( country , "_Daily_" , type , "_Zoom_", date_last , ".tiff") 
      )
      
      ggsave( list_vis[[type]][[country]][["Daily_Case_Zoom"]] , 
              file   = dzpath_tmp , 
              width  = 7.29      , 
              height = 4.51      , 
              unit   = "in"      , 
              dpi    = res
      )  
      
    }
  }
}


