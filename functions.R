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


find_equilibrium <- function(prod_df, delta_misoprene, MWseed, mseed, Temp) {
  
  reltol <- 0.001
  
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
    new_prod_df <- old_prod_df %>%
      mutate(
        iter = iter+1,
        Fi = Kp*Mtot*totalmass/(1+Kp*Mtot))
    
    final_df <- bind_rows(final_df, new_prod_df)
    
    if (all(abs(new_prod_df$Fi-old_prod_df$Fi)/old_prod_df$Fi < reltol)) {
      convergence <- TRUE
    }
    
    old_prod_df <- new_prod_df
  }
  
  return(final_df)

}


prod_df <- tibble::tibble(
  iter = 0,
  Name = c("Prod1"),
  yield = c(1.0),
  pvap = c(1e-10), 
  molw = c(250),
  accoeff = c(1)
)

prod_df <- tibble::tibble(
  iter = 0,
  Name = c("Prod1","Prod2","Prod3"),
  yield = c(1.0, 0.5, 0.1),
  pvap = c(1e-10, 1e-13, 1e-7), 
  molw = c(100, 150, 200),
  accoeff = c(1, 1, 1)
)

prod_df <- tibble::tibble(
  iter = 0,
  Name = c("Prod1","Prod2","Prod3","Prod4", "Prod5","Prod6", "Prod7","Prod8"),
  yield = c(1.0, 0.5, 0.1, 0.2, 0.4, 0.7,0.8, 0.3),
  pvap = c(1e-10, 1e-13, 1e-7, 1e-5, 1e-12, 1e-11, 1e-9, 1e-7), 
  molw = c(100, 150, 200, 250, 120, 140, 190, 180),
  accoeff = 1
)


MWseed <- 250
mseed <- 10
delta_misoprene <- 1
Temp <- 298
mseed_sensitivity <- NULL
for (mseed in c(1:9 %o% 10^(-2:2))){
  
  mseed_sensitivity <- bind_rows(mseed_sensitivity, 
                                 find_equilibrium(prod_df, delta_misoprene, MWseed, mseed, Temp) %>%
                                   mutate(Mseed = mseed)
  )
}

SOAyield <- mseed_sensitivity %>%
  group_by(Mseed) %>%
  filter(iter == max(iter)) %>%
  mutate(SOAyield = totalmass*Fi)

ggplot(SOAyield) +
  geom_area(aes(x =Mseed, y = SOAyield, fill = Name))+
  scale_x_log10(name = "Mseed [ug m-3]")
