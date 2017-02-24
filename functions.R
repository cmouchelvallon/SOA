library(tidyverse)

meanMW <- function(MWseed, mseed, MWsoa, msoa) {
  return((MWseed*mseed + MWsoa*msoa)/(mseed+msoa))
}

Kp <- function(Temp, MW, pvap, accoeff) {
  return((8.2e-5*Temp)/(1e6*MW*accoeff*pvap))
}

to_minimise_M <- function(M, Kp, Tmass ){
  
}
prod_df <- tibble::tibble(
  iter = 0,
  Name = c("Prod1","Prod2","Prod3"),
  yield = c(1.0, 0.5, 0.1),
  pvap = c(1e-10, 1e-13, 1e-7), 
  molw = c(100, 150, 200),
  accoeff = c(1, 1, 1)
)

MWseed <- 250
mseed <- 10
delta_misoprene <- 1
Temp <- 298

Figuess <- 0.1
prod_df <- prod_df %>%
  mutate(
    totalmass = yield*delta_misoprene,
    Fi = Figuess)

convergence <- FALSE
# while(!convergence){
  SOA <- prod_df %>%
    mutate(aermass = totalmass*Fi) %>%
    
    summarise(MSOA = sum(aermass),
              MWsoa = weighted.mean(molw, aermass*Fi)
    )
  
  MW <- meanMW(MWseed, mseed, SOA$MWsoa, SOA$MSOA)
  prod_df <- prod_df %>%
    mutate(Kp = Kp(Temp, MW, pvap, accoeff))
  
  
  M_before <- SOA$MSOA + mseed
  optim(M_before)
# }
