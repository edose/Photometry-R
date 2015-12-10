##### Model.R   Get model (via mixed-model regression) for one Astronight's Comp and other standard stars.
##### TODO : write function model(), whose ,
#####          and whose output is a model object from lmer(). 

model <- function (df_master, filter="V", maxUncert=0.02, fit_transform=FALSE, fit_extinction=TRUE) {
  # Input is the Astronight's master data frame (or CSV file of it).
  # Returns lmer model object (of class merMod) of fit to comparison stars only.
  require(dplyr)
  require(magrittr)
  df_model <- df_master %>%
    filter(Filter==filter) %>%
    filter(MagUncertainty<=maxUncert) %>%
    filter(Saturated==FALSE) %>%
    filter(StarType=="Comp")
  
  formula_string <- "InstMag ~ offset(CatMag) + JD_mid + "
  offset<- rep(0,nrow(df_model))
  if (fit_transform) {
    formula_string <- paste(formula_string, ""      sep="")
  } else {
    transform <- list(V=-0.0259,R=+0.0319,I=+0.0064)[filter] %>% unlist() # default (not fit) value.
  }
  if (fit_extinction) {
    
  } else {
    extinction <- list(V=0.2,R=0.1,I=0.08)[filter] %>% unlist() # default (not fit) value.
  }
  
  
  
  return (thisModel)
}