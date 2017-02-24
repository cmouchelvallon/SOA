library(tidyverse)

meanMW <- function(MWseed, mseed, MWsoa, msoa) {
  return((MWseed*mseed + MWsoa*msoa)/(mseed+msoa))
}

Kp <- function(Temp, MW, pvap, accoeff) {
  return((8.2e-5*Temp)/(1e6*MW*accoeff*pvap))
}

to_minimise_M <- function(M, Kp, Tmass ){
  return(abs(1-sum(Kp*Tmass/(1+Kp*M))))
  
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

reltol <- 0.01

Figuess <- 0.1
old_prod_df <- prod_df %>%
  mutate(
    totalmass = yield*delta_misoprene,
    Fi = Figuess)

final_df <- old_prod_df

convergence <- FALSE
while(!convergence){
  SOA <- old_prod_df %>%
    mutate(aermass = totalmass*Fi) %>%
    
    summarise(MSOA = sum(aermass),
              MWsoa = weighted.mean(molw, aermass*Fi)
    )
  
  MW <- meanMW(MWseed, mseed, SOA$MWsoa, SOA$MSOA)
  old_prod_df <- old_prod_df %>%
    mutate(Kp = Kp(Temp, MW, pvap, accoeff))
  
  
  M_before <- SOA$MSOA
  
  M <- (optimize(f = to_minimise_M, interval = c(0,mseed*2), Kp = old_prod_df$Kp, Tmass = old_prod_df$totalmass))$minimum
  
  Mtot <- M + mseed
  print(abs(M-M_before)/M_before)
  if(abs(M-M_before)/M_before <= reltol) {
    convergence <- TRUE
  }
  new_prod_df <- old_prod_df %>%
    mutate(
      iter = iter+1,
      Fi = Kp*Mtot*totalmass/(1+Kp*Mtot))
  
  final_df <- bind_rows(final_df, new_prod_df)
  
  old_prod_df <- new_prod_df
}
