##### Model.R   Get model (via mixed-model regression) for one Astronight's Comp and other standard stars.

modelAll <- function(df_master, maxMagUncertainty=0.02, 
                     fit_transform=FALSE, fit_extinction=TRUE, fit_starID=FALSE) {
  filters <- df_master$Filter %>% unique()
  masterModelList <- list()
  print(filters)
  for (filter in filters) {
    filterList <- modelOneFilter(df_master, filter, maxMagUncertainty, 
                                 fit_transform, fit_extinction, fit_starID)
    masterModelList[[filter]] <- filterList
  }
  return (masterModelList)
}

modelOneFilter <- function (df_master, filter="V", maxMagUncertainty=0.02, 
                            fit_transform=FALSE, fit_extinction=TRUE, fit_starID=FALSE) {
  # Input is the Astronight's master data frame.
  # Returns lmer model object (of class merMod) of fit to comparison stars only.
  require(dplyr)
  require(magrittr)
  df_model <- df_master %>%
    filter(StarType=="Comp") %>%
    filter(Filter==filter) %>%
    filter(UseInModel==TRUE) %>%
    filter(MagUncertainty<=maxMagUncertainty) %>%
    filter(Saturated==FALSE) %>%
    filter(!is.na(CI)) %>%
    filter(CI<=1.6)

  formula_string <- "InstMag ~ offset(CatMag) + (1|JD_mid)"
  thisOffset <- rep(0,nrow(df_model))
  if (fit_transform) {
    transform <- NA
    formula_string <- paste0(formula_string, " + CI")
  } else {
    transform <- list(V=-0.0259,R=+0.0319,I=+0.0064)[filter] %>% unlist() # default (not fit) value.
    thisOffset <- thisOffset + transform * df_model$CI
  }
  if (fit_extinction) {
    extinction <- NA
    formula_string <- paste0(formula_string, " + Airmass")
  } else {
    extinction <- list(V=0.2,R=0.1,I=0.08)[filter] %>% unlist() # default (not fit) value.
    thisOffset <- thisOffset + extinction * df_model$Airmass
  }
  if (fit_starID) {
    formula_string <- paste0(formula_string, " + (1|ModelStarID)")
  }
  
  thisFormula <- as.formula(formula_string)
  df_model <- df_model %>%
    mutate(JD_mid=as.factor(JD_mid)) %>%
    mutate(ModelStarID=as.factor(ModelStarID))
  
  # Run model.
  require(lme4)  
  thisModel <- lmer (thisFormula, data=df_model, REML=FALSE, offset=thisOffset)
  
  # Initial stats & diagnostics.

  return (list(model=thisModel, filter=filter, transform=transform, extinction=extinction))
}

################################################################################################
##### Below are support-only functions, not called by user. ####################################

