#
### Libraries
#
remotes::install_github("covid19datahub/R")

library("COVID19")
library("dplyr")
library("rjags")
library("MLmetrics")
library("doRNG")
library("foreach")
library("doFuture")


#
### Parameters
## Data
date_last                      <- Sys.Date() - 7    # Data is loaded with one week (7 days) lag
date_cutoff                    <- date_last  - 1    # Separates training and test sets (7 days to have a week of testing)
date_pred                      <- Sys.Date() - 2    # Date to predict 

## Cutoffs
confirmed_cases_cutoff         <- 1000              # Filters out countries with less confirmed cases than
fatalities_cutoff              <- 25                # Filters out countries with less fatalities than
density_confirmed_cases_cutoff <- 0.001             # Filters out countries with less density of confirmed cases than
population_cutoff              <- 1e6               # Filters out countries with less population than

## Plotting 
plot_flag      <- 0        # 0 : no output. This flag also rules the writting of summaries.
plot_size      <- 500      # in px
plot_font_size <- 22       # in px
plot_quality   <- 75       # jpeg

## Projections and credible intervals
projected_days <- 2*(as.numeric(Sys.Date() - as.Date("2020-01-01"))) # Entire range of days (data + prediction)
credible_int   <- c(0.025, 0.975)   # 95% credible interval 

## MCMC variables
num_chains     <- 3
num_adap       <- 1000

niter_update   <- 1e3      # Recommended 1e3
niter_sim_init <- 1e3      # Number of MCMC steps for the first iteration
niter_sim_max  <- 1e7      # Exits with (geq). Recommended: 1e7
thinning_init  <- 1        # 100
epsilon_0      <- 2        # Maximum change in the position of the inflection point. Exits with (leq)


#
### Function to cut off countries using a criterion based on a minimum number of
#   "Fatalities", "Confirmed_Cases", or "Pop+Dens" (Density of Population and 
#    Population)
#
exclude_countries <- function(df,type){
  if (type == "Fatalities"){
    countries_out <- setdiff( unique(df$Country) , unique(subset(df, Fatalities > fatalities_cutoff)$Country) )
  }
  else if (type == "Confirmed_Cases"){
    countries_out <- setdiff( unique(df$Country) , unique(subset(df, Confirmed_Cases > confirmed_cases_cutoff)$Country) )
  }
  else if (type == "Pop+Dens"){
    countries_out <- setdiff( unique(df$Country) , unique(subset(df, Density_Confirmed_Cases > density_confirmed_cases_cutoff & Population > population_cutoff)$Country) )    
  }
  else{
    stop(type,"is not a valid type: Must be Fatalities, Confirmed_Cases, or Pop+Dens")
  }
  invisible(countries_out)
}


include_countries <- function(df, type){
  if (type == "Fatalities"){
    countries_in <-as.character(sort(unique(subset(df, Fatalities > fatalities_cutoff)$Country)))
  }
  else if (type == "Confirmed_Cases"){
    countries_in <-as.character(sort(unique(subset(df, Confirmed_Cases > confirmed_cases_cutoff)$Country)))
  }
  else if (type == "Pop+Dens"){
    countries_in <-as.character(sort(unique(subset(df, Density_Confirmed_Cases > density_confirmed_cases_cutoff & Population > population_cutoff)$Country)))
  }
  else{
    stop(type,"is not a valid type: Must be Fatalities, Confirmed_Cases, or Pop+Dens")
  }
  invisible(countries_in)
}


#
### Function to build the model string for the MCMC simulation w/JAGS
#
build_model_string <- function(type, mean_logasymptote, mean_location, mean_growth_rate, prec_logasymptote=1e-4, prec_location=1e-4, prec_growth_rate=1e-4){
  model_string = paste0(" model {

    for (i in 1:length(",type,")) {
        ",type,"[i] ~ dnorm(mu[i], prec)
        mu[i] = asymptote*exp(-location*exp(-growth_rate*Time[i]))
    }
    
        asymptote = exp(logasymptote)
  
        logasymptote  ~ dnorm(",mean_logasymptote,",", prec_logasymptote," )     
        location      ~ dnorm(",mean_location,    ",", prec_location,    " )     
        growth_rate   ~ dnorm(",mean_growth_rate, ",", prec_growth_rate, " )     
        
        prec          ~ dgamma(1.0, 1.0)
} "
  )
  return(model_string)
}

#
### Function to run the MCMC simulation w/JAGS
#
run_mcmc_simulation <- function(index_country, 
                                model_string, 
                                model_parameters, 
                                number_chains = num_chains, 
                                number_adap = num_adap, 
                                number_thin_init = thinning_init, 
                                number_steps_update = niter_update, 
                                number_steps_run_init = niter_sim_init, 
                                number_steps_run_max = niter_sim_max, 
                                epsilon_init = projected_days, 
                                epsilon_max = epsilon_0){
  
  print(">>>>>>>>>>>>>>>>>>>>>>>>>>>")
  print("")
  print(paste("Country",by_country[index_country]))
  print("")
  print(">>>>>>>>>>>>>>>>>>>>>>>>>>>")
  
  
  set.seed(1970)
  
  # Subset the country
  jags_datafile <- as.list( subset(train_by_country , Country == by_country[index_country]) )
  
  # Change in the location of the inflection point
  inf_pt_init <- epsilon_init  
  epsilon     <- epsilon_init
  niter_sim   <- number_steps_run_init
  thinning    <- number_thin_init
  
  # Begin while loop to converge the MCMC steps
  while (epsilon >= epsilon_max && niter_sim <= number_steps_run_max ){
    
    model_jags <- jags.model(
      textConnection(model_string) ,
      data     = jags_datafile ,
      n.chains = number_chains ,
      n.adap   = number_adap
    )
    
    update(
      model_jags, 
      number_steps_update
    ) 
    
    model_sim <- coda.samples(
      model          = model_jags ,
      variable.names = model_parameters ,
      n.iter         = niter_sim ,
      thin           = thinning ,
      quiet          = TRUE
    )
    
    model_csim <- do.call(rbind, model_sim) 
    
    # Calculate the inflection point and distance
    inf_pt  <-inflection_point(model_csim)
    epsilon <- abs( inf_pt - inf_pt_init )
    
    # Monitor 
    print("HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH")
    print(paste("number of MCMC steps",niter_sim))
    print(paste("Distance to convergence",epsilon))
    print(paste("Thinning every",thinning,"sample"))
    print("HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH")
    
    # Update the position of the inflection point and the number of iterations
    inf_pt_init <- inf_pt
    niter_sim   <- 10*niter_sim
    
    # Update the thinning as function of niter_sim
    if (niter_sim > 1e3 && niter_sim < 1e5){
      thinning <- 10*thinning
    }
    
    if (niter_sim > 1e4 && niter_sim < 1e6){
      thinning <- 10*thinning
    }
    
    return_list<-list(model_sim, model_csim, niter_sim)
  } # end while loop
  return(return_list)
}




#
### Functions to get characteristics of the Gompertz curve and its derivative
#
growth_rate <- function(csim_data_gompertz){
  colMeans(csim_data_gompertz)[[1]] 
}

location <- function(csim_data_gompertz){
  colMeans(csim_data_gompertz)[[2]] 
}

asymptote <- function(csim_data_gompertz){
  exp( colMeans(csim_data_gompertz)[[3]] )
}

inflection_point <- function(csim_data_gompertz){
  log( location(csim_data_gompertz) ) / growth_rate(csim_data_gompertz)
}

peak_daily_rate <- function(csim_data_gompertz){
  as.character( as.Date(day_zero) + round(inflection_point(csim_data_gompertz), 0) ) 
}

days_first_fatality <- function(csim_data_gompertz, date_1st_fatality){
  round(inflection_point(csim_data_gompertz), 0) - date_1st_fatality 
}

days_first_confirmed_case <- function(csim_data_gompertz, date_1st_confirmed_case){
  round(inflection_point(csim_data_gompertz), 0) - date_1st_confirmed_case 
}

gompertz_eval <- function(csim_data_gompertz, date_eval){
  if (class(date_eval) == "Date"){
    date_eval <- as.numeric(as.Date(date_eval)-as.Date(day_zero), "days") 
  }
  exp(colMeans(csim_data_gompertz)[[3]]) * exp( -colMeans(csim_data_gompertz)[[2]] * exp(-as.numeric(colMeans(csim_data_gompertz)[[1]], units="days") * date_eval) )
}


#
### Functions to get credible intervals from MCMC posteriors
#
model_credible_interval <- function(csim_data_gompertz, date_pred){
  if (class(date_pred) == "Date"){
    date_pred <- as.numeric(as.Date(date_pred)-as.Date(day_zero), "days") 
  } 
  quantile(exp(csim_data_gompertz[,3]) * exp(-csim_data_gompertz[,2] * exp(-csim_data_gompertz[,1] * date_pred )), 
           probs=credible_int, names = FALSE  )
}

deriv_model_credible_interval <- function(csim_data_gompertz, date_pred){
  if (class(date_pred) == "Date"){
    date_pred <- as.numeric(as.Date(date_pred)-as.Date(day_zero), "days") 
  } 
  quantile(exp( csim_data_gompertz[,3] ) * csim_data_gompertz[,1] * csim_data_gompertz[,2] / (exp(csim_data_gompertz[,1] * date_pred + csim_data_gompertz[,2] * exp(-csim_data_gompertz[,1] * date_pred ) )),
           probs=credible_int, names = FALSE )
}



#
### Functions to create and fill the summary dataframe
#
create_dummy_summary_df <- function(){
  data.frame(
    Country          = rep("", length(by_country)) ,
    Growth_Rate      = rep(NA, length(by_country)) ,
    Location         = rep(NA, length(by_country)) , 
    Log_Asymptote    = rep(NA, length(by_country)) ,
    FATorConfCases   = rep(NA, length(by_country)) ,
    Peak_daily_rate  = rep("", length(by_country)) ,
    Days_1st_Fat     = rep(NA, length(by_country)) ,
    Days_1st_CCase   = rep(NA, length(by_country)) ,
    Predicted_Day    = rep(NA, length(by_country)) ,
    Lower_CI         = rep(NA, length(by_country)) ,
    Upper_CI         = rep(NA, length(by_country)) ,
    MSE              = rep(NA, length(by_country)) ,           
    MSPE             = rep(NA, length(by_country)) ,
    Precision        = rep(NA, length(by_country)) ,
    Niter            = rep(NA, length(by_country)) ,
    stringsAsFactors = FALSE
  )       
  
}

summary_by_country <- function(type, csim_data_gompertz, country, date_pred, data_train, data_test, converged_niter_sim){
  # Fatalities or Confirmed_Cases
  if (type == "Fatalities"){
    y_obs_fit  <- subset(data_train, Country==by_country[i])$Fatalities
    y_obs_pred <- subset(data_test,  Country==by_country[i])$Fatalities 
  }
  else if (type == "Confirmed_Cases"){
    y_obs_fit  <- subset(data_train, Country==by_country[i])$Confirmed_Cases
    y_obs_pred <- subset(data_test,  Country==by_country[i])$Confirmed_Cases
  }
  else{
    stop(type,"is not a valid type: Must be Fatalities or Confirmed_Cases")
  }
   
  # Fitting and prediction x and y from data
  x_obs_fit  <- as.numeric(subset(data_train, Country==by_country[i])$Time, "days")
  x_obs_pred <- as.numeric(subset(data_test,  Country==by_country[i])$Time, "days")
  

  # Fitting and prediction from model
  y_fit  <- gompertz_eval(csim_data_gompertz, x_obs_fit)
  y_pred <- gompertz_eval(csim_data_gompertz, x_obs_pred)
  
  # Date of first Fatality / Confirmed Case for country[i]
  date_first_fatality        <- subset(data_train, Fatalities > 0      & Country==by_country[i] )$Time[1]
  date_first_confirmed_case  <- subset(data_train, Confirmed_Cases > 0 & Country==by_country[i] )$Time[1] 
  
  list(
    country ,                                       # Country
    growth_rate(csim_data_gompertz) ,               # Growth_Rate
    location(csim_data_gompertz) ,                  # Location
    log(asymptote(csim_data_gompertz)) ,            # Log_Asymptote
    round(asymptote(csim_data_gompertz), 0) ,       # Asymptote, rounded to nearest integer
    peak_daily_rate(csim_data_gompertz) ,           # Peak_daily_rate
    days_first_fatality(csim_data_gompertz, date_first_fatality) ,             # Days_1st_Fat
    days_first_confirmed_case(csim_data_gompertz, date_first_confirmed_case) , # Days_1st_CCase
    
    # Round to nearest integer to avoid fractional days
    round(gompertz_eval(csim_data_gompertz, date_pred), 0) ,                 # Pred at date_pred
    round(model_credible_interval(csim_data_gompertz, date_pred)[1], 0) ,    # Lower CI
    round(model_credible_interval(csim_data_gompertz, date_pred)[2], 0) ,    # Upper CI
    
    MSE(y_fit /asymptote(csim_data_gompertz) ,  y_obs_fit/asymptote(csim_data_gompertz) ) ,  # MSE
    MSE(y_pred/asymptote(csim_data_gompertz) , y_obs_pred/asymptote(csim_data_gompertz) ) ,  # MSPE
    
    colMeans(csim_data_gompertz)[[4]] , # Precision
    
    # Final number of iterations for a converged inflection point locus
    # The factor of 1/10 is necessary to get the last value executed in the while loop
    converged_niter_sim/10                   
  )
}


plot_gompertz_validation <- function(type, data_train, data_test, formula, csim_data_gompertz, projected_days, plot_flag = 0, label_x_axis = paste("Days elapsed since",as.character(day_zero))){
  
  # Date of first Fatality / Confirmed Case for country[i]
  date_first_fatality        <- subset(data_train, Fatalities > 0      & Country==by_country[i] )$Time[1]
  date_first_confirmed_case  <- subset(data_train, Confirmed_Cases > 0 & Country==by_country[i] )$Time[1] 
  
  # Max value of the response variable
  if (type == "Fatalities"){
    max_y <- max(subset(data_train, Country==by_country[i] )$Fatalities)
  }
  else if (type == "Confirmed_Cases"){
    max_y <- max(subset(data_train, Country==by_country[i] )$Confirmed_Cases)
  }
  else {
    max_y <- 0
  }
  
  # Credible intervals polygon
  xaxis_ci       <- as.numeric(1:projected_days) 
  lower_ci       <- sapply(xaxis_ci, model_credible_interval,       csim_data_gompertz = csim_data_gompertz)[1,]
  upper_ci       <- sapply(xaxis_ci, model_credible_interval,       csim_data_gompertz = csim_data_gompertz)[2,]
  lower_ci_deriv <- sapply(xaxis_ci, deriv_model_credible_interval, csim_data_gompertz = csim_data_gompertz)[1,]
  upper_ci_deriv <- sapply(xaxis_ci, deriv_model_credible_interval, csim_data_gompertz = csim_data_gompertz)[2,]
  
  
  # To output plot
  if (plot_flag == 1){
    # Fatalities or Confirmed_Cases
    if (type == "Fatalities"){
      plot_file = paste0(results_dir,"Plots/",by_country[i],"_","fatalities_v_", date_cutoff,".jpg")
      jpeg(plot_file, width = plot_size, height = plot_size, units = "px", pointsize = plot_font_size, quality = plot_quality)
    }
    else if (type == "Confirmed_Cases"){
      plot_file = paste0(results_dir,"Plots/",by_country[i],"_","confirmed_cases_v_", date_cutoff,".jpg")
      jpeg(plot_file, width = plot_size, height = plot_size, units = "px", pointsize = plot_font_size, quality = plot_quality)
      
    }
    else{
      stop(type,"is not a valid type: Must be Fatalities or Confirmed_Cases")
    }
  }
  
  plot(
    formula,
    data     = subset(data_train, Country==by_country[i]),
    type     = "b",
    pch      = 22 ,
    bg       = "dark green" ,
    col      = "dark green" ,
    
    xlab     = label_x_axis,
    # ylab     = "Fatalities" ,
    
    xlim     = c( 0, projected_days ) ,
    ylim     = c( 0, 1.05*max(upper_ci, max_y) )
  )
  
  points(
    formula ,
    data=subset(data_test, Country==by_country[i]) ,
    type = "b" ,
    pch  = 22  ,
    bg   = "dark red" ,
    col  = "dark red"
  )
  
  curve( 
    gompertz_eval(csim_data_gompertz, x) ,
    lwd = 2, 
    lty = 1, 
    col = 'red', 
    add = TRUE 
  )
  
  # Locus of the inflection point
  abline(
    v = inflection_point(csim_data_gompertz) ,
    lwd = 4, 
    lty = 3, 
    col = adjustcolor("grey", alpha.f = 0.9)
  )
  
  # Locus of the first fatality
  abline(
    v = date_first_fatality ,
    lwd = 4, 
    lty = 3, 
    col = adjustcolor("purple", alpha.f = 0.5)
  )
  
  # Locus of the first confirmed case
  abline(
    v = date_first_confirmed_case ,
    lwd = 4, 
    lty = 3, 
    col = adjustcolor("blue", alpha.f = 0.5)
  )
  
  polygon(
    x   = c(xaxis_ci,rev(xaxis_ci)), 
    y   = c(lower_ci,rev(upper_ci)), 
    col = adjustcolor("red", alpha.f = 0.20), 
    border = NA
  )
  
  legend(
    "bottomright",
    legend = c( by_country[i], round(asymptote(csim_data_gompertz),0) ) ,
    bty    = "n",
    cex    = 1.0
  )
  
  if (plot_flag == 1){
    dev.off()
  }
  
}

plot_gompertz <- function(type, data_train, formula, csim_data_gompertz, projected_days, plot_flag = 0, label_x_axis = paste("Days elapsed since",as.character(day_zero))){
 
  # Date of first Fatality / Confirmed Case for country[i]
  date_first_fatality        <- subset(data_train, Fatalities > 0      & Country==by_country[i] )$Time[1]
  date_first_confirmed_case  <- subset(data_train, Confirmed_Cases > 0 & Country==by_country[i] )$Time[1] 
  
  # Max value of the response variable
  if (type == "Fatalities"){
    max_y <- max(subset(data_train, Country==by_country[i] )$Fatalities)
  }
  else if (type == "Confirmed_Cases"){
    max_y <- max(subset(data_train, Country==by_country[i] )$Confirmed_Cases)
  }
  else {
    max_y <- 0
  }
  
  # Credible intervals polygon
  xaxis_ci       <- as.numeric(1:projected_days) 
  lower_ci       <- sapply(xaxis_ci, model_credible_interval, csim_data_gompertz = csim_data_gompertz)[1,]
  upper_ci       <- sapply(xaxis_ci, model_credible_interval, csim_data_gompertz = csim_data_gompertz)[2,]
  
  # To output plot
  if (plot_flag == 1){
    # Fatalities or Confirmed_Cases
    if (type == "Fatalities"){
      plot_file = paste0(results_dir,"Plots/",by_country[i],"_","fatalities_", date_cutoff,".jpg")
      jpeg(plot_file, width = plot_size, height = plot_size, units = "px", pointsize = plot_font_size, quality = plot_quality)
    }
    else if (type == "Confirmed_Cases"){
      plot_file = paste0(results_dir,"Plots/",by_country[i],"_","confirmed_cases_", date_cutoff,".jpg")
      jpeg(plot_file, width = plot_size, height = plot_size, units = "px", pointsize = plot_font_size, quality = plot_quality)
      
    }
    else{
      stop(type,"is not a valid type: Must be Fatalities or Confirmed_Cases")
    }
  }
  
  plot(
    formula,
    data     = subset(data_train, Country==by_country[i]),
    type     = "b",
    pch      = 22 ,
    bg       = "dark green" ,
    col      = "dark green" ,
    
    xlab     = label_x_axis,
    
    xlim     = c( 0, projected_days ) ,
    ylim     = c( 0, 1.05*max(upper_ci, max_y) )
  )
  
  curve( 
    gompertz_eval(csim_data_gompertz, x) ,
    lwd = 2, 
    lty = 1, 
    col = 'red', 
    add = TRUE 
  )
  
  # Locus of the inflection point
  abline(
    v = inflection_point(csim_data_gompertz) ,
    lwd = 4, 
    lty = 3, 
    col = adjustcolor("grey", alpha.f = 0.9)
  )
  
  # Locus of the first fatality
  abline(
    v = date_first_fatality ,
    lwd = 4, 
    lty = 3, 
    col = adjustcolor("purple", alpha.f = 0.5)
  )
  
  # Locus of the first confirmed case
  abline(
    v = date_first_confirmed_case ,
    lwd = 4, 
    lty = 3, 
    col = adjustcolor("blue", alpha.f = 0.5)
  )
  
  polygon(
    x   = c(xaxis_ci,rev(xaxis_ci)), 
    y   = c(lower_ci,rev(upper_ci)), 
    col = adjustcolor("red", alpha.f = 0.20), 
    border = NA
  )
  
  legend(
    "bottomright",
    legend = c( by_country[i], round(asymptote(csim_data_gompertz),0) ) ,
    bty    = "n",
    cex    = 1.0
  )
  
  if (plot_flag == 1){
    dev.off()
  }
}


plot_gompertz_daily_cases_validation <- function(type, data_train, data_test, formula, csim_data_gompertz, projected_days, plot_flag = 0, label_x_axis = paste("Days elapsed since",as.character(day_zero))){

  # Date of first Fatality / Confirmed Case for country[i]
  date_first_fatality        <- subset(data_train, Fatalities > 0      & Country==by_country[i] )$Time[1]
  date_first_confirmed_case  <- subset(data_train, Confirmed_Cases > 0 & Country==by_country[i] )$Time[1] 
  
  # Credible intervals polygon
  xaxis_ci       <- as.numeric(1:projected_days) 
  lower_ci_deriv <- sapply(xaxis_ci, deriv_model_credible_interval, csim_data_gompertz = csim_data_gompertz)[1,]
  upper_ci_deriv <- sapply(xaxis_ci, deriv_model_credible_interval, csim_data_gompertz = csim_data_gompertz)[2,]
  
  # To output plot
  if (plot_flag == 1){
    # Fatalities or Confirmed_Cases
    if (type == "Fatalities"){
      plot_file = paste0(results_dir,"Plots/",by_country[i],"_","daily_fatalities_v_", date_cutoff,".jpg")
      jpeg(plot_file, width = plot_size, height = plot_size, units = "px", pointsize = plot_font_size, quality = plot_quality)
    }
    else if (type == "Confirmed_Cases"){
      plot_file = paste0(results_dir,"Plots/",by_country[i],"_","daily_confirmed_cases_v_", date_cutoff,".jpg")
      jpeg(plot_file, width = plot_size, height = plot_size, units = "px", pointsize = plot_font_size, quality = plot_quality)
      
    }
    else{
      stop(type,"is not a valid type: Must be Fatalities or Confirmed_Cases")
    }
  }
  
  plot(
    formula,
    data     = subset(data_train, Country==by_country[i]),
    type     = "b",
    pch      = 22 ,
    bg       = "dark green" ,
    col      = "dark green" ,
    xlab     = label_x_axis ,
    xlim     = c( 0, projected_days )
  )
  
  points(
    formula ,
    data=subset(data_test, Country==by_country[i]) ,
    type = "b" ,
    pch  = 22  ,
    bg   = "dark red" ,
    col  = "dark red"
  )
  
  curve(
    gompertz_eval(csim_data_gompertz, x) * growth_rate(csim_data_gompertz) * location(csim_data_gompertz) * exp(-growth_rate(csim_data_gompertz) * x ) ,
    lwd  = 2 , 
    lty  = 1 , 
    col  = 'red',
    add = TRUE
  )
  
  polygon(  
    x   = c(xaxis_ci,rev(xaxis_ci)) ,
    y   = c(lower_ci_deriv,rev(upper_ci_deriv)) , 
    col = adjustcolor("red", alpha.f = 0.2), 
    border = NA
  )
  
  # Locus of the inflection point
  abline(
    v = inflection_point(csim_data_gompertz) ,
    lwd = 4, 
    lty = 3, 
    col = adjustcolor("grey", alpha.f = 0.9)
  )
  
  # Locus of the first fatality
  abline(
    v = date_first_fatality ,
    lwd = 4, 
    lty = 3, 
    col = adjustcolor("purple", alpha.f = 0.5)
  )
  
  # Locus of the first confirmed case
  abline(
    v = date_first_confirmed_case ,
    lwd = 4, 
    lty = 3, 
    col = adjustcolor("blue", alpha.f = 0.5)
  )
  
  legend(
    "topright",
    legend = c( by_country[i], round( asymptote(csim_data_gompertz) * growth_rate(csim_data_gompertz) * exp(-1), 0) ) ,
    bty    = "n",
    cex    = 1.0
  )
  
  if (plot_flag == 1){
    dev.off()
  } 
}

plot_gompertz_daily_cases <- function(type, data_train, formula, csim_data_gompertz, projected_days, plot_flag = 0, label_x_axis = paste("Days elapsed since",as.character(day_zero))){
 
  # Date of first Fatality / Confirmed Case for country[i]
  date_first_fatality        <- subset(data_train, Fatalities > 0      & Country==by_country[i] )$Time[1]
  date_first_confirmed_case  <- subset(data_train, Confirmed_Cases > 0 & Country==by_country[i] )$Time[1] 
  
  ### The daily cases
  xaxis_ci       <- as.numeric(1:projected_days) 
  lower_ci_deriv <- sapply(xaxis_ci, deriv_model_credible_interval, csim_data_gompertz = csim_data_gompertz)[1,]
  upper_ci_deriv <- sapply(xaxis_ci, deriv_model_credible_interval, csim_data_gompertz = csim_data_gompertz)[2,]
  
  # To output plot
  if (plot_flag == 1){
    # Fatalities or Confirmed_Cases
    if (type == "Fatalities"){
      plot_file = paste0(results_dir,"Plots/",by_country[i],"_","daily_fatalities_", date_cutoff,".jpg")
      jpeg(plot_file, width = plot_size, height = plot_size, units = "px", pointsize = plot_font_size, quality = plot_quality)
    }
    else if (type == "Confirmed_Cases"){
      plot_file = paste0(results_dir,"Plots/",by_country[i],"_","daily_confirmed_cases_", date_cutoff,".jpg")
      jpeg(plot_file, width = plot_size, height = plot_size, units = "px", pointsize = plot_font_size, quality = plot_quality)
      
    }
    else{
      stop(type,"is not a valid type: Must be Fatalities or Confirmed_Cases")
    }
  }
  
  
  plot(
    formula,
    data     = subset(data_train, Country==by_country[i]),
    type     = "b",
    pch      = 22 ,
    bg       = "dark green" ,
    col      = "dark green" ,
    xlab     = label_x_axis ,
    xlim     = c( 0, projected_days )
  )
  
  curve(
    gompertz_eval(csim_data_gompertz, x) * growth_rate(csim_data_gompertz) * location(csim_data_gompertz) * exp(-growth_rate(csim_data_gompertz) * x ) ,
    lwd  = 2 , 
    lty  = 1 , 
    col  = 'red',
    add = TRUE
  )
  
  polygon(  
    x   = c(xaxis_ci,rev(xaxis_ci)) ,
    y   = c(lower_ci_deriv,rev(upper_ci_deriv)) , 
    col = adjustcolor("red", alpha.f = 0.2), 
    border = NA
  )
  
  # Locus of the inflection point
  abline(
    v = inflection_point(csim_data_gompertz) ,
    lwd = 4, 
    lty = 3, 
    col = adjustcolor("grey", alpha.f = 0.9)
  )
  
  # Locus of the first fatality
  abline(
    v = date_first_fatality ,
    lwd = 4, 
    lty = 3, 
    col = adjustcolor("purple", alpha.f = 0.5)
  )
  
  # Locus of the first confirmed case
  abline(
    v = date_first_confirmed_case ,
    lwd = 4, 
    lty = 3, 
    col = adjustcolor("blue", alpha.f = 0.5)
  )
  
  
  legend(
    "topright",
    legend = c( by_country[i], round( asymptote(csim_data_gompertz) * growth_rate(csim_data_gompertz) * exp(-1)  ,0 ) ) ,
    bty    = "n",
    cex    = 1.0
  )
  
  if (plot_flag == 1){
    dev.off()
  } 
} 



