##### Model.R   Get model (via mixed-model regression) for one Astronight's Comp and other standard stars.
#####  Modified June 2016 to offer an option to set a max allowed error (from AAVSO chart),
#####    in which case the models use only standards and high-quality comps.
#####  
##### Typical usages:
#####   listV <- modelOneFilter(AN_rel_folder="20151216", filter="V", maxCatMagError=0.01)
#####   [or to get star effect:] modelOneFilter(AN_rel_folder="20151216", filter="V", fit_starID=TRUE)$star 
#####        %>% arrange(desc(abs(CatalogEffect)))
#####   For each model, curate its input points by editing the omit.txt.
#####   make_masterModelList(AN_rel_folder="20151216", modelLists=list(listV, listR, listI))

modelOneFilter <- function (AN_top_folder="J:/Astro/Images/C14", AN_rel_folder=NULL, 
                            filter=NULL, maxInstMagSigma=0.03, maxColorIndex=2.5, saturatedADU=54000,
                            maxCatMagError=NULL,  # NULL means "don't use to filter observations"
                            fit_skyBias=TRUE, fit_vignette=TRUE, fit_XY=FALSE,
                            fit_transform=FALSE, fit_extinction=TRUE, fit_starID=FALSE) {
  # Inputs are: (1) the Astronight's master data frame (as stored in /Photometry), and
  #             (2) the omit.txt file of observations to omit (also stored in /Photometry).
  # Returns model list, including: lmer model, df_obs, df_image, scalar results.
  # Usage if run manually: listV <- modelOneFilter(AN_rel_folder="20151216", filter="V")
  require(dplyr, quietly=TRUE)
  require(magrittr, quietly=TRUE)
  if (is.null(AN_rel_folder)) {stop(">>>>> You must provide a AN_rel_folder parm, ",
                                  "e.g., AN_rel_folder='20151216'.")}
  if (is.null(filter)) {stop(">>>>> You must provide a filter parm, e.g., filter='V'.")}

  df_model <- make_df_model(AN_top_folder, AN_rel_folder, 
                            filter, maxInstMagSigma, maxColorIndex, saturatedADU,
                            maxCatMagError)
  formula_string <- "InstMag ~ offset(CatMag) + (1|JD_mid)"
  thisOffset <- rep(0,nrow(df_model))
  if (fit_transform) {
    # Case of fitting transform values from the current data.
    transform <- NA
    formula_string <- paste0(formula_string, " + CI")
  } else {
    # Case of user-given (not fitted) transform value.
    # These improved values are from standards run of AN 20151218.
    transform <- list(V=-0.0199,R=+0.0398,I=+0.0432)[filter] %>% unlist()
    thisOffset <- thisOffset + transform * df_model$CI
    cat(">>>>> Transform (Color Index) not fit: value fixed at ", transform, "\n", sep="")
    
  }
  if (fit_extinction) { 
    extinction <- NA
    formula_string <- paste0(formula_string, " + Airmass")
  } else {
    extinction <- list(V=0.167,R=0.128,I=0.103)[filter] %>% unlist() # user-given (fit values for DSW site)
    thisOffset <- thisOffset + extinction * df_model$Airmass
    cat(">>>>> Extinction (Airmass) not fit: value fixed at ", extinction, "\n", sep="")
  }
  if (fit_skyBias) {
    if (all(df_model$SkyBias==0)) {
      cat(">>>>> ALL SkyBias values in df_model's ", nrow(df_model), 
          " observations are 0, so removing SkyBias term from model.\n", sep="")
    } else {
      formula_string <- paste0(formula_string, " + SkyBias")
    }
  }
  if (fit_vignette) {
    formula_string <- paste0(formula_string, " + Vignette")
  }
  if (fit_XY) {
    formula_string <- paste0(formula_string, " + X1024 + Y1024")
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
  image <- ranef(thisModel)$"JD_mid"
  image$JD_mid <- row.names(image)        # dplyr drops row names silently (and too greedily).
  row.names(image) <- NULL
  names(image)[names(image)=="(Intercept)"] <- "CirrusEffect"  # rename column (dplyr fails on this).
  image <- image %>%
    left_join(df_model[,c("FITSfile","JD_mid")] %>% unique(), by="JD_mid") %>%
    select(JD_mid, everything())

  # Construct output data frame "star" (one row per comp star included in model, zeroes otherwise).
  # ... need an included observation count per star, too.
  df_starCount <- data.frame(ModelStarID=df_model$ModelStarID, stringsAsFactors=FALSE) %>%
    group_by(ModelStarID) %>%
    summarize(nObs=n()) %>%
    as.data.frame()  
  if (fit_starID) {
    star <- ranef(thisModel)$"ModelStarID"
    star$ModelStarID <- row.names(star)  # dplyr drops row names silently (and too greedily).
    row.names(star) <- NULL
    names(star)[names(star)=="(Intercept)"] <- "CatalogEffect"  # rename column (dplyr fails on this).
    star <- star %>% 
      left_join(df_starCount, by="ModelStarID") %>%
      select(ModelStarID, CatalogEffect, everything())
  } else {
    star <- df_starCount %>% 
      mutate(CatalogEffect=0) %>%
      select(ModelStarID, CatalogEffect, everything())
  }
  
  # Extract scalar results.
  extinction <- ifelse(fit_extinction, fixef(thisModel)["Airmass"],   extinction)
  transform  <- ifelse(fit_transform,  fixef(thisModel)["CI"],        transform)
  vignette   <- ifelse(fit_vignette,   fixef(thisModel)["Vignette"],  0)

  modelList <- list(model=thisModel, obs=obs, AN=AN_rel_folder, image=image, star=star, filter=filter, 
                    transform=transform, fit_transform=fit_transform,
                    extinction=extinction, fit_extinction=fit_extinction,
                    vignette=vignette)
  cat("modelOneFilter('", filter, "') completed on ", nrow(df_model), " observations.\n", sep="")
  # source('C:/Dev/Photometry/Plots.R')
  # modelPlots(modelList)
  print(summary(thisModel))
  cat(nrow(df_model), "observations -->", 
      "sigma =", format(1000*sigma(thisModel), nsmall=1, digits=3), "mMag.")
  return (modelList)
}

make_masterModelList <- function(AN_top_folder="J:/Astro/Images/C14", AN_rel_folder=NULL, modelLists=NULL) {
  ##### The last stage in modeling, run once all comp models are OK. 
  #####    Saves to AN's Photometry folder, in one list of lists (the "masterModelList"). Nothing returned.
  ##### Needs testing (20160119).
  ##### Typical usage: saveAllModels(AN_rel_folder="20151216", modelLists=list(listV,listI))
  
  if (is.null(AN_rel_folder)) {stop(">>>>> You must provide a AN_rel_folder parm, ",
                                  "e.g., AN_rel_folder='20151216'.")}
  if (is.null(modelLists)) {stop(">>>>> You must provide a list of modelLists, ",
                                       "e.g., modelLists=list(listV,listI).")}
  if (class(modelLists[[1]])!="list") {stop(">>>>> You must enter modelLists as =list() and NOT =c().")}
  
  masterModelList <- list()
  for (modelList in modelLists) {
    filter <- modelList$filter
    masterModelList[[filter]] <- modelList
  }
  AN_folder   <- make_safe_path(AN_top_folder, AN_rel_folder)
  photometry_folder <- make_safe_path(AN_folder, "Photometry")
  masterModelList_path <- make_safe_path(photometry_folder, "masterModelList.Rdata")
  save(masterModelList, file=masterModelList_path, precheck=FALSE)
  
  # Make a template-only pre-predict.txt file if one doesn't already exist.
  pre_predictPath <- make_safe_path(photometry_folder, "pre-predict.txt")
  if (!file.exists(pre_predictPath)) {
    lines <- c(
      paste0(";----- This is pre-predict.txt for AN folder ", AN_rel_folder),
      paste0(";----- Select comp stars (by FOV, filter, & StarID) from input to predictAll()."),
      paste0(";----- Example directive lines:\n;"),
      paste0(";#COMPS  Obj, V, 132, 133 144    ; to KEEP from FOV 'Obj': ",
             "comp stars '132' '133' and '144' in filter 'V'"),
      paste0(";\n;----- Add your directive lines:\n;\n\n")
    )
    writeLines(lines, con=pre_predictPath)
  }
  cat("saveAllModels() has saved this AN's masterModelList to:\n   ", masterModelList_path, 
      "   and has written a 'pre-predict.txt' template file.",
      "\n   This AN is now ready for predictAll().\n")
}


################################################################################################
##### Below are test or support-only functions, rarely or not typically called by user. ########

make_df_model <- function (AN_top_folder="J:/Astro/Images/C14", AN_rel_folder=NULL, 
                           filter=NULL, maxInstMagSigma=0.03, maxColorIndex=2.5, saturatedADU=54000,
                           maxCatMagError=NULL) {
  
  df_omitted <- omitObs(AN_top_folder, AN_rel_folder) # returns df w/user-requested obs removed.
  if (is.na(df_omitted[[1]][[1]]))  {stop(">>>>> omit.txt file does not exist.")}
  df_model <- df_omitted %>%
    filter(Filter==filter) %>%
    filter(StarType=="Comp") %>%
    filter(!is.na(CatMag)) %>%
    filter(!is.na(CI)) %>%
    filter(!is.na(Airmass)) %>%
    filter(InstMagSigma<=maxInstMagSigma) %>%
    filter(CI<=maxColorIndex) %>%
    filter(MaxADU_Ur<=saturatedADU)
  eligibleCatMagErrors <- df_model$CatMagError[!is.na(df_model$CatMagError)]
  if (is.null(maxCatMagError)) {
    hist(eligibleCatMagErrors, breaks=16, 
         main=paste("Filter =", filter, "   //   median", median(eligibleCatMagErrors)))
  } else {
    df_model <- df_model %>% filter(CatMagError <= maxCatMagError)
  }
  return (df_model)
}

omitObs <- function (AN_top_folder="J:/Astro/Images/C14", AN_rel_folder) {
  ##### Reads AN folder's df_master and omit.txt, returns df_filtered with requested observations omitted.
  ##### Typically called by Model.R::modelOneFilter() and Predict.R::predictAll().
  ##### Tested OK 20160117.
  ##### Typical usage: omitObs(AN_rel_folder="20151216")
  require(stringi, quietly=TRUE)
  require(dplyr, quietly=TRUE)
  source('C:/Dev/Photometry/$Utility.R')
  
  AN_folder <- make_safe_path(AN_top_folder, AN_rel_folder)
  photometry_folder <- make_safe_path(AN_folder, "Photometry")
  
  df_master <- get_df_master(AN_top_folder, AN_rel_folder) # begin with df_master, all rows.
  df_filtered <- df_master
  cat("omitObs() begins with", nrow(df_filtered), "rows, ")
  
  # Get omit file for this AN folder, then parse directive lines.
  omit_file <- make_safe_path(photometry_folder, "omit.txt")
  if (!file.exists(omit_file)) {
    cat("\n>>>>> Omit file", omit_file, "does not exist.\n\n")
    return (NA)
  } else {
    lines <- readLines(omit_file, warn=FALSE)   # read last line even without EOL character(s).
    for (iLine in 1:length(lines)) {
      lines[iLine] <- lines[iLine] %>% 
        strsplit(";",fixed=TRUE) %>% unlist() %>% first() %>% trimws()  # remove comments
    }
    # Detect and collect directive text lines.
    directiveLines <- lines[stri_detect_regex(lines,'^#')]
    directiveLines <- directiveLines[!is.na(directiveLines)]
  }
  
  # Process directive lines to curate df_omitted (i.e., to omit observations etc as requested).
  for (thisLine in directiveLines) {
    directive <- thisLine %>% strsplit("[ \t]",fixed=FALSE) %>% unlist() %>% 
      first() %>% trimws() %>% toupper()
    value     <- thisLine %>% substring(nchar(directive)+1) %>% trimws()  # all but the directive
    if (is.na(directive)) {directive <- ""}
    if (directive=="#OBS") {
      parms <- value %>% strsplit(",",fixed=TRUE) %>% unlist() %>% trimws()
      if (length(parms)==2) {
        FITSname <- paste0(parms[1], ".fts")
        starID <- parms[2]
        df_filtered <- df_filtered %>% filter(!(FITSfile==FITSname & StarID==starID))
      } else {
        cat(paste(">>>>> Can't parse line: ", thisLine))
      }
    }
    if (directive=="#STAR") {
      parms <- value %>% strsplit(",",fixed=TRUE) %>% unlist() %>% trimws()
      if (length(parms) %in% 2:3) {
        object <- parms[1]
        starID <- parms[2]
        filter <- ifelse(length(parms)==3, parms[3], NA)
        if (is.na(filter)) {
          df_filtered <- df_filtered %>% filter(!(FOV==object & StarID==starID))
        } else {
          df_filtered <- df_filtered %>% filter(!(FOV==object & StarID==starID & Filter==filter))
        }
      } else {
        cat(paste(">>>>> Can't parse line: ", thisLine))
      }
    }
    if (directive=="#SERIAL") {
      serials_to_omit <- stri_extract_all_words(value) %>% 
        unlist() %>% strsplit(",") %>% unlist() %>% as.numeric()
      df_filtered <- df_filtered %>% filter(!Serial %in% serials_to_omit)
    }
    if (directive=="#IMAGE") {
      FITSname <- paste0(value, ".fts")
      df_filtered <- df_filtered %>% filter(!FITSfile==FITSname)
    }
    if (directive=="#JD") {
      parms <- value %>% strsplit(",",fixed=TRUE) %>% unlist() %>% trimws()
      if (length(parms)==2) {
        minJD <- as.numeric(parms[1])
        maxJD <- as.numeric(parms[2])
        if (minJD>0 | maxJD < 2) {
          floorJD <- floor(min(as.numeric(df_filtered$JD_mid)))
          df_filtered <- df_filtered %>% 
            filter((JD_mid <= minJD+floorJD) | (JD_mid >= maxJD+floorJD))
        } else {
          cat(paste(">>>>> minJD and/or maxJD values suspect: ", thisLine))
        }
      } else {
        cat(paste(">>>>> Can't parse line: ", thisLine))
      }
    }
  }
  cat("ends with", nrow(df_filtered), "rows.\n")
  if (TRUE) {  # set to FALSE (or remove) when testing has been completed.
    thisDiff <- setdiff(df_master, df_filtered)
    cat("omitObs() removed", nrow(thisDiff), "observations.\n")
    # if (nrow(thisDiff)>=1) { print(thisDiff %>% select(Serial, FITSfile, StarID)) }
  }
  return (df_filtered)
}

get_df_master <- function(AN_top_folder="J:/Astro/Images/C14", AN_rel_folder) {
  source("C:/Dev/Photometry/$Utility.R")
  AN_folder   <- make_safe_path(AN_top_folder, AN_rel_folder)
  photometry_folder <- make_safe_path(AN_folder, "Photometry")
  df_master_path <- make_safe_path(photometry_folder, "df_master.Rdata")
  load(df_master_path)
  return (df_master)
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
