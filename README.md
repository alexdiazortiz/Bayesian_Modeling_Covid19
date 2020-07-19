# Bayesian Modeling of Covid19 Spread
A Bayesian approach of the modeling COVID19 spread based on the Gompertz equation. 

It uses
* the R Interface to COVID-19 Data Hub  of Guidotti and Ardia to retrive current COVID19 data worldwide. 

   - See : https://cran.r-project.org/web/packages/COVID19/index.html and "COVID-19 Data Hub," Journal of Open Source Software, 5(51), 2376. doi: 10.21105/joss.02376.
* RJAGS package of Martyn Plummer to perform MCMC to obtain posterior estimates 
   - See : https://cran.r-project.org/web/packages/rjags/index.html
   
   
It performs:
* A nonlinear parametric fit of the COVID19 spread (i.e., fatalities and confirmed cases) by country.

It provides:
* A dataframe summarizing the model parameters by country for further analysis.
* Spread (fatalities and confirmed cases) plots.
