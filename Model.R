##### Model.R   Get model (via mixed-model regression) for one Astronight's Comp and other standard stars.
#####  Tests OK 20160124.
##### Typical usages: listV <- modelOneFilter(AN_rel_folder="20151216", filter="V")
#####                 omitSerial(AN_rel_folder="20151216", serial=c(123,32))
#####                 make_masterModelList(AN_rel_folder="20151216", modelLists=c(listV, listR, listI))

modelOneFilter <- function (AN_top_folder="J:/Astro/Images/C14", AN_rel_folder=NULL, 
                            filter=NULL, maxMagUncertainty=0.03, maxColorIndex=2.5, saturatedADU=54000,
                            fit_vignette=TRUE, fit_vignette4=FALSE, fit_XY=FALSE,
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

  df_model <- omitObs(AN_top_folder, AN_rel_folder) %>% # returns df w/user-requested obs removed.
    filter(StarType=="Comp") %>%
    filter(Filter==filter) %>%
    filter(UseInModel==TRUE) %>%
    filter(MagUncertainty<=maxMagUncertainty) %>%
    filter(!is.na(CI)) %>%
    filter(CI<=maxColorIndex) %>%
    filter(MaxADU_Ur<=saturatedADU) %>%
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
  vignette4  <- ifelse(fit_vignette4,  fixef(thisModel)["Vignette4"], 0)
  
  modelList <- list(model=thisModel, obs=obs, AN=AN_rel_folder, image=image, star=star, filter=filter, 
                    transform=transform, fit_transform=fit_transform,
                    extinction=extinction, fit_extinction=fit_extinction,
                    vignette=vignette, vignette4=vignette4)
  cat("modelOneFilter('", filter, "') completed on ", nrow(df_model), " observations.\n", sep="")
  # source('C:/Dev/Photometry/Plots.R')
  # modelPlots(modelList)
  print(summary(thisModel))
  cat(nrow(df_model), "observations -->", 
      "sigma =", format(1000*sigma(thisModel), nsmall=1, digits=3), "mMag.")
  return (modelList)
}

omitSerial <- function (AN_top_folder="J:/Astro/Images/C14", AN_rel_folder=NULL, serial=NULL) {
  ##### User utility to quickly add lines to omit.txt from Serial numbers (df_master) only.
  ##### Tested OK 20160117.
  ##### Typical usage: omitSerial(AN_rel_folder="20151216", serial=c(1220,923,88,727))
  
  if (is.null(AN_rel_folder)) {stop(">>>>> You must provide a AN_rel_folder parm, ",
                                  "e.g., AN_rel_folder='20151216'.")}
  if (is.null(serial)) {
    stop(">>>>> No serial values given.")
  }
  require(stringi, quietly=TRUE)
  require(dplyr, quietly=TRUE)
  AN_folder <- make_safe_path(AN_top_folder, AN_rel_folder)
  photometry_folder <- make_safe_path(AN_folder, "Photometry")
  df_master <- get_df_master(AN_top_folder, AN_rel_folder)
  if(any(!(serial %in% df_master$Serial))) {
    stop(">>>>> at least one serial number given is not in df_master for ", AN_rel_folder)
  }
  lines <- paste0(";\n;Next ", length(serial), " lines added via omitSerial() at ", Sys.time(), "...\n")
  for (thisSerial in serial) {
    thisFITS <- df_master$FITSfile[df_master$Serial==thisSerial] %>%
      strsplit(".f", fixed=TRUE) %>% unlist() %>% first() %>% trimws()
    thisStarID <- df_master$StarID[df_master$Serial==thisSerial]
    lines <- c(lines, paste0("#OBS   ", thisFITS, ", ", thisStarID, "  ;  Serial=", thisSerial, "\n"))
  }
  omitPath <- make_safe_path(photometry_folder, "omit.txt")
  cat(lines, file=omitPath, sep="", append=TRUE)
  cat(lines) # output to Console
}

make_masterModelList <- function(AN_top_folder="J:/Astro/Images/C14", AN_rel_folder=NULL, modelLists=NULL) {
  ##### The last stage in modeling, run once all comp models are OK. 
  #####    Saves to AN's Photometry folder, in one list of lists (the "masterModelList"). Nothing returned.
  ##### Needs testing (20160119).
  ##### Typical usage: saveAllModels(AN_rel_folder="20151216", modelLists=list(listV,listI))
  
  if (is.null(AN_rel_folder)) {stop(">>>>> You must provide a AN_rel_folder parm, ",
                                  "e.g., AN_rel_folder='20151216'.")}
  if (is.null(modelLists)) {stop(">>>>> You must provide a vector of modelLists, ",
                                       "e.g., vectorOfModelLists=list(listV,listI).")}
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
  cat("saveAllModels() has saved this AN's masterModelList to:\n   ", masterModelList_path, 
      "\n   Now ready for predictAll().\n")
}


################################################################################################
##### Below are test or support-only functions, rarely or not typically called by user. ########

omitObs <- function (AN_top_folder="J:/Astro/Images/C14", AN_rel_folder) {
  ##### Reads AN folder's df_master and omit.txt, returns df_filtered with requested observations omitted.
  ##### Typically called by Model.R::modelOneFilter() and Predict.R::predictOneFilter().
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
  load(df_master_path, verbose=TRUE)
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
