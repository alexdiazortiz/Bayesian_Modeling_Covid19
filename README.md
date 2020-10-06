# Bayesian Modeling of Covid19 Spread
A Bayesian approach of the modeling COVID19 one and two-wave spread based on the Gompertz equation via Markov Chain Monte Carlo (MCMC) simulations.

It needs:
* JAGS (https://sourceforge.net/projects/mcmc-jags/)
* RStudio ( https://rstudio.com)
* R (https://www.r-project.org)

It uses:
* the R Interface to COVID-19 Data Hub  of Guidotti and Ardia to retrive current COVID19 data worldwide. 

   - See : https://cran.r-project.org/web/packages/COVID19/index.html and "COVID-19 Data Hub," Journal of Open Source Software, 5(51), 2376. doi: 10.21105/joss.02376.
   
* RJAGS package of Martyn Plummer to perform MCMC to obtain posterior estimates 
   - See : https://cran.r-project.org/web/packages/rjags/index.html
   
* Packages for parallel computations (see full list in the "utilities" file):
   - foreach
   - doFuture
   - doRNG
     
It performs:
* A nonlinear parametric fit of the COVID19 spread (i.e., fatalities, confirmed cases, recovered) by country.

It provides:
* Both daily and cumulative spreading (fatalities, confirmed cases, recovered and active cases) plots. 
