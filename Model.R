##### Model.R   Get model (via mixed-model regression) for one Astronight's Comp and other standard stars.
#####  In testing 20151228.
##### Typical usage: m <- modelAll(df_master)

modelAll <- function(df_master, maxMagUncertainty=0.02, maxColorIndex=2,
                     vignette=TRUE, vignette4=FALSE,
                     fit_transform=FALSE, fit_extinction=TRUE, fit_starID=FALSE) {
  require(dplyr)
  filters <- df_master$Filter %>% unique()
  masterModelList <- list()
  for (filter in filters) {
    filterList <- modelOneFilter(df_master, filter, maxMagUncertainty, maxColorIndex,
                                 vignette, vignette4,
                                 fit_transform, fit_extinction, fit_starID)
    masterModelList[[filter]] <- filterList
  }
  return (masterModelList)
}

################################################################################################
##### Below are test or support-only functions, rarely or not typically called by user. ########

modelOneFilter <- function (df_master, filter="V", maxMagUncertainty=0.02, maxColorIndex=2,
                            fit_vignette=TRUE, fit_vignette4=FALSE,
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
    filter(!is.na(Airmass)) %>%
    filter(!is.na(CatMag)) %>%
    filter(!is.na(CI)) %>%
    filter(CI<=maxColorIndex)
  
  formula_string <- "InstMag ~ offset(CatMag) + (1|JD_mid)"
  thisOffset <- rep(0,nrow(df_model))
  if (fit_transform) {
    transform <- NA
    formula_string <- paste0(formula_string, " + CI")
  } else {
    transform <- list(V=-0.0259,R=+0.0319,I=+0.0064)[filter] %>% unlist() # user-given (not fitted) value.
    thisOffset <- thisOffset + transform * df_model$CI
  }
  if (fit_extinction) {
    extinction <- NA
    formula_string <- paste0(formula_string, " + Airmass")
  } else {
    extinction <- list(V=0.2,R=0.1,I=0.08)[filter] %>% unlist() # user-given (not fitted) value.
    thisOffset <- thisOffset + extinction * df_model$Airmass
  }
  if (fit_vignette) {
    formula_string <- paste0(formula_string, " + Vignette")
  }
  if (fit_vignette4) {
    formula_string <- paste0(formula_string, " + Vignette4")
  } 
  if (fit_starID) {
    formula_string <- paste0(formula_string, " + (1|ModelStarID)")
  }
  
  # Run model.
  thisFormula <- as.formula(formula_string)
  df_model <- df_model %>%
    mutate(JD_mid=as.factor(JD_mid)) %>%
    mutate(ModelStarID=as.factor(ModelStarID))  # need factors for lmer() random effects.
  require(lme4)  
  thisModel <- lmer (thisFormula, data=df_model, REML=FALSE, offset=thisOffset)
  df_model <- df_model %>%
    mutate(JD_mid=as.character(JD_mid)) %>%
    mutate(ModelStarID=as.character(ModelStarID))  # undo the factoring done just before running model.
  
  # Construct output data frame "obs" (one row per observation included in model).
  obs <- df_model %>% 
    mutate(Residual=residuals(thisModel)) %>%
    mutate(Fitted=fitted(thisModel))

  # Construct output data frame "image" (one row per image included in model).
  # image <- ranef(thisModel)$"JD_mid"   # extract from lmer::ranef() list. 
  image <- ranef(thisModel)$"JD_mid"
  image$JD_mid <- row.names(image)        # dplyr drops row names silently (and too greedily).
  row.names(image) <- NULL
  names(image)[names(image)=="(Intercept)"] <- "Cirrus"  # rename column (dplyr fails on this).
  image <- image %>% select(JD_mid, everything())
  # Also need to capture (from df_model) per-image: the image name, airmass, #obs used, object, exposure.

  
  # Construct output data frame "star" (one row per comp star included in model).
  # ... need an included observation count per star, too.
  
  
  # Extract scalars & parameter values.
  sigma <- sigma(thisModel)
  ZeroPoint  <- fixef$"(Intercept)"
  extinction <- ifelse(fit_extinction, fixef$Airmass,   extinction)
  transform  <- ifelse(fit_transform,  fixef$CI,        transform)
  vignette   <- ifelse(fit_vignette,   fixef$Vignette,  0)
  vignette4  <- ifelse(fit_vignette4,  fixef$Vignette4, 0)
  
  # Initial stats & diagnostics.
  
  
  
  
  return (list(model=thisModel, obs=obs, image=image, star=star, 
               filter=filter, transform=transform, extinction=extinction))
  
  # return (list(model=thisModel, obs=df_model, filter=filter, transform=transform, extinction=extinction))
}

mock_df_master <- function (df_master, stdevMag=0.01) {
  require(dplyr)
  ZeroPoint  <- list(I=-19,R=-20,V=-21)             %>% unlist() # arbitrary but realistic.
  Extinction <- list(I=0.08,R=0.1,V=0.2)            %>% unlist() # same as modelAll() defaults.
  Transform  <- list(V=-0.0259,R=+0.0319,I=+0.0064) %>% unlist() # same as modelAll() defaults.
  vignetteCorner <- list(I=0.09,R=0.1,V=+0.06)      %>% unlist() # InstMag increase at image corner.
  df_mock <- df_master %>%
    mutate(InstMagSaved=InstMag) %>%
    mutate(InstMag=CatMag + 
             ZeroPoint[Filter] + 
             Extinction[Filter]*Airmass + 
             Transform[Filter]*CI + 
             vignetteCorner[Filter]*Vignette + 
             rnorm(InstMag,0,stdevMag))
  return(df_mock)
}
