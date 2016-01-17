##### Model.R   Get model (via mixed-model regression) for one Astronight's Comp and other standard stars.
#####  In testing 20151228.
##### Typical usage: m <- modelAll(df_master)

modelAll <- function(AN_top_folder="J:/Astro/Images/C14", AN_rel_folder, 
                     maxMagUncertainty=0.02, maxColorIndex=2,
                     vignette=TRUE, vignette4=FALSE,
                     fit_transform=FALSE, fit_extinction=TRUE, fit_starID=FALSE) {
  require(dplyr, quietly=TRUE)
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

omitObs <- function (AN_top_folder="J:/Astro/Images/C14", AN_rel_folder) {
  ##### Reads AN folder's df_master and omit.txt, returns df_filtered with requested observations omitted.
  ##### Typically called by Model.R::modelOneFilter() and Predict.R::predictOneFilter().
  ##### Tested OK 20160117.
  ##### Typical usage: omitObs(AN_rel_folder="20151216")
  require(stringi, quietly=TRUE)
  require(dplyr, quietly=TRUE)
  AN_folder <- make_safe_path(AN_top_folder, AN_rel_folder)
  photometry_folder <- make_safe_path(AN_folder, "Photometry")
  
  df_master <- get_df_master(AN_top_folder, AN_rel_folder) # begin with df_master, all rows.
  df_filtered <- df_master
  cat("omitObs() begins with", nrow(df_filtered), "rows.\n")

  # Get omit file for this AN folder, then parse directive lines.
  omit_file <- make_safe_path(photometry_folder, "omit.txt")
  if (!file.exists(omit_file)) {
    cat("Omit file", omit_file, "does not exist.")
    return (NA)
  } else {
    lines <- readLines(omit_file, warn=FALSE)   # read last line even without EOL character(s).
    for (iLine in 1:length(lines)) {
      lines[iLine] <- lines[iLine] %>% 
        strsplit(";",fixed=TRUE) %>% unlist() %>% first() %>% trimws()  # remove comments
    }
    directiveLines <- lines[stri_detect_regex(lines,'^#')] # detect and collect directive text lines.
  }
  
  # Process directive lines to curate df_omitted (i.e., to omit observations etc as requested).
  for (thisLine in directiveLines) {
    directive <- thisLine %>% strsplit("[ \t]",fixed=FALSE) %>% unlist() %>% 
      first() %>% trimws() %>% toupper()
    value     <- thisLine %>% substring(nchar(directive)+1) %>% trimws()  # all but the directive
    if (is.na(directive)) {directive <- ""}
    if (directive=="#OBS") {
      parms <- value %>% strsplit(",",fixed=TRUE) %>% unlist() %>% trimws()
      if (length(parms)>=2) {
        FITSname <- paste0(parms[1], ".fts")
        starID <- parms[2]
        df_filtered <- df_filtered %>% filter(!(FITSfile==FITSname & StarID==starID))
      }
    }
    if (directive=="#STAR") {
      parms <- value %>% strsplit(",",fixed=TRUE) %>% unlist() %>% trimws()
      if (length(parms) %in% 2:3) {
        object <- parms[1]
        starID <- parms[2]
        filter <- ifelse(length(parms)==3, parms[3], NA)
        if (is.na(filter)) {
          df_filtered <- df_filtered %>% filter(!(Object==object & StarID==starID))
        } else {
          df_filtered <- df_filtered %>% filter(!(Object==object & StarID==starID & Filter==filter))
        }
      }
    }
    if (directive=="#IMAGE") {
      FITSname <- paste0(value, ".fts")
      df_filtered <- df_filtered %>% filter(!FITSfile==FITSname)
    }
  }
  cat("omitObs() ends with", nrow(df_filtered), "rows.\n")
  if (TRUE) {  # set to FALSE (or remove) when testing has been completed.
    thisDiff <- setdiff(df_master, df_filtered)
    cat(nrow(thisDiff), "rows were removed:\n")
    print(thisDiff %>% select(Serial, FITSfile, StarID))
  }
  
  return (df_filtered)
}


################################################################################################
##### Below are test or support-only functions, rarely or not typically called by user. ########

get_df_master <- function(AN_top_folder="J:/Astro/Images/C14", AN_rel_folder) {
  source("C:/Dev/Photometry/$Utility.R")
  AN_folder   <- make_safe_path(AN_top_folder, AN_rel_folder)
  photometry_folder <- make_safe_path(AN_folder, "Photometry")
  df_master_path <- make_safe_path(photometry_folder, "df_master.Rdata")
  load(df_master_path, verbose=TRUE)
  return (df_master)
}

modelOneFilter <- function (AN_top_folder="J:/Astro/Images/C14", AN_rel_folder, 
                            filter, maxMagUncertainty=0.02, maxColorIndex=2, saturatedADU=54000,
                            fit_vignette=TRUE, fit_vignette4=FALSE,
                            fit_transform=FALSE, fit_extinction=TRUE, fit_starID=FALSE) {
  # Input is the Astronight's master data frame.
  # Returns lmer model object (of class merMod) of fit to comparison stars only.
  require(dplyr, quietly=TRUE)
  require(magrittr, quietly=TRUE)

  df_model <- omitObs(AN_top_folder, AN_rel_folder) %>% # omitObs() returns df w/user-requested obs removed.
    filter(StarType=="Comp") %>%
    filter(Filter==filter) %>%
    filter(UseInModel==TRUE) %>%
    filter(MagUncertainty<=maxMagUncertainty) %>%
    filter(!is.na(CI)) %>%
    filter(CI<=maxColorIndex) %>%
    filter(MaxADU<=saturatedADU) %>%
    filter(!is.na(Airmass)) %>%
    filter(!is.na(CatMag))
  
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
  require(dplyr, quietly=TRUE)
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
