##### PredictTargets.R   Predict magnitudes of Target (unknown) stars in Astronight's master data frame.
#####    Uses this Astronight's lmer() model object etc from Model.R::modelAll().
##### V 0.12 tested OK 20160213.
##### Typical workflow:
#####    Ensure masterModelList is ready to go from ListV, etc via Model.R::make_masterModelList().
#####    df_predictions <- predictAll(AN_rel_folder="20151216")
#####    AAVSO(AN_rel_folder="20151216", software_version="0.0.0")
#####    Examine AAVSO report, edit report_map.txt for #SERIAL & #COMBINE directives, rerun AAVSO().
#####    Submit/upload AAVSOreport-nnnnnnnn.txt to AAVSO; check for proper upload.
#####    Set all /Photometry files to read only (in Windows).

predictAll <- function (AN_top_folder="J:/Astro/Images/C14", AN_rel_folder=NULL, 
                        saturatedADU=54000, maxInstMagSigma=0.05, CI_filters=c("V","I")) {
  ##### Performs ALL steps to get transformed magnitudes for Check and Target observations in df_master.
  ##### Returns big data frame df_predict (which is transformed).
  ##### Requires that AN folder contain R files df_master.Rdata (from Input.R::make_df_master()) and 
  #####    masterModelList.Rdata (from Model.R::make_masterModelList()).
  ##### Typical usage: df_predictions <- predictAll(AN_rel_folder="20151216")
  if (is.null(AN_rel_folder)) {stop(">>>>> You must provide a AN_rel_folder, ",
                                  "e.g., AN_rel_folder='20151216'.")}
  source("C:/Dev/Photometry/$Utility.R")
  AN_folder   <- make_safe_path(AN_top_folder, AN_rel_folder)
  photometry_folder <- make_safe_path(AN_folder, "Photometry")
  path_df_master <- make_safe_path(photometry_folder, "df_master.Rdata")
  load(path_df_master)
  path_masterModelList <- make_safe_path(photometry_folder, "masterModelList.Rdata")
  load(path_masterModelList)
  
  # Select Target and Check rows for which CatMag to be predicted (all filters).
  # Make df of all needed inputs (newdata) for all filters & all (untransformed) lme4::predict() calls.
  require(dplyr, quietly=TRUE)
  df_predict_input <- df_master %>%
    filter(StarType %in% c("Check","Target")) %>%
    filter(MaxADU_Ur<=saturatedADU) %>%
    filter(InstMagSigma<=maxInstMagSigma) %>%
    mutate(CI=ifelse(is.na(CI),0,CI)) %>%
    select(Serial, ModelStarID, StarID, Chart, Xcentroid, Ycentroid, InstMag, InstMagSigma, StarType,
           JD_mid, Filter, Airmass, CI, Vignette) %>%
    mutate(CatMag=0) # arbitrarily chosen, to complete model.

  # Get untransformed predicted Inst Mags (via running predict(), collecting all results).
  df_predictions <- data.frame()
  require(lme4, quietly = TRUE)
  for (thisModelList in masterModelList) {
    allowed_JD_mid <- unique(as.character(model.frame(thisModelList$model)$JD_mid))
    df_thisFilter <- df_predict_input %>% 
      filter(Filter==thisModelList$filter) %>%
      filter(JD_mid %in% allowed_JD_mid) # don't predict for JDs outside model.
    predictOutput <- predict(thisModelList$model, newdata=df_thisFilter, re.form= ~(1|JD_mid))
    # Add extinction*airmass term if not included in model (as formula term "Airmass").
    if (!thisModelList$fit_extinction) {
      predictOutput <- predictOutput + extinction * df_thisFilter$Airmass
    }
    df_predictions_thisFilter <- df_thisFilter %>%
      mutate(PredictedInstMag = predictOutput)
    df_predictions <- rbind(df_predictions, df_predictions_thisFilter)
  }

  # Impute UntransformedMags from PredictedInstMag and InstMag (and CatMag=0).
  df_predictions <- df_predictions %>%
    mutate(UntransformedMag = InstMag - PredictedInstMag)
  
  # Impute Color Index values for target observations by interpolation.
  transforms <- c(masterModelList[[CI_filters[1]]]$transform, masterModelList[[CI_filters[2]]]$transform)
  df_predictions <- imputeTargetCIs(df_predictions, CI_filters, transforms) # Fill in CI values for targets.
  
  # Perform color transformations.
  df_xref_transform <- data.frame(stringsAsFactors=FALSE) # build a lookup table to get Transforms.
  for (thisModelList in masterModelList) {
    df_xref_transform <- df_xref_transform %>%
      rbind(data.frame(Filter=thisModelList$filter, Transform=thisModelList$transform))
  }
  df_xref_transform$Filter <- as.character(df_xref_transform$Filter)
  df_transformed <- df_predictions %>%
    left_join(df_xref_transform, by="Filter") %>%
    mutate(TransformedMag = UntransformedMag - Transform * CI) %>% # minus is the correct sign: 20160206.
    filter(!is.na(TransformedMag))  # NA often from missing V or I mag, so CI=NA, so Transformed Mag=NA.

  # Compute MagErr (for AAVSO).
  df_xref_modelSigma <- data.frame(stringsAsFactors=FALSE) # build a lookup table to get model mag errs
  for (thisModelList in masterModelList) {
    df_xref_modelSigma <- df_xref_modelSigma %>%
      rbind(data.frame(Filter=thisModelList$filter, ModelSigma=sigma(thisModelList$model)))
  }
  
  # Finish construction of data frame.
  df_xref_modelSigma$Filter <- as.character(df_xref_modelSigma$Filter)
  df_transformed <- df_transformed %>%
    left_join(df_xref_modelSigma, by="Filter") %>%
    mutate(MagErr = pmax(ModelSigma, InstMagSigma)) %>% # pmax = "parallel" element-wise max of 2 vectors
    left_join((df_master %>% select(Serial, FITSfile, Exposure)), by="Serial") %>% # restore a few columns.
    arrange(ModelStarID, JD_num)
  
  # Save df_transformed as .Rdata file (but do not return it).
  path_transformed <- make_safe_path(photometry_folder, "df_transformed.Rdata")
  save(df_transformed, file=path_transformed, precheck=FALSE)
  
  # Make a template-only report_map.txt file if report_map.txt doesn't already exist.
  path_report_map <- make_safe_path(photometry_folder, "report_map.txt")
  if (!file.exists(path_report_map)) {
    lines <- c(
      paste0(";----- This is report_map.txt for AN folder ", AN_rel_folder),
      paste0(";----- Use this file to omit and/or combine target observations from AAVSO report."),
      paste0(";----- Example directive lines:"),
      paste0(";"),
      paste0(";#TARGET  GD Cyg ; to omit this target star altogether from AAVSO report."),
      paste0(";#SERIAL  34 44,129  32  1202 ; to omit these 5 Serial numbers from AAVSO report."),
      paste0(";#COMBINE   80,128 ; to combine (average) these 2 Serial numbers within AAVSO report."),
      paste0(";----- Add your directive lines:"),
      paste0(";\n") 
    )
    writeLines(lines, con=path_report_map)
  }
  nTargetObs <- df_transformed %>% filter(StarType=="Target") %>% nrow()
  nTargets   <- df_transformed %>% filter(StarType=="Target") %>% select(ModelStarID) %>% 
    unique() %>% nrow()
  nCheckObs  <- df_transformed %>% filter(StarType=="Check") %>% nrow()
  Filters    <- df_transformed %>% filter(StarType=="Target") %>% select(Filter) %>% 
    unique()
  cat("PredictAll() yields", nTargetObs, "Target obs for", nTargets, "targets in filters:", 
      paste(Filters,collapse=" ","\n"))
  cat("   and ", nCheckObs, " Check observations.\n")

  cat("Transformed predictions saved to", path_transformed, "\n",
      "   and report_map.txt is ensured available in the same folder\n",
      "Now you are ready to:\n",
      "   1. run targetPlots(),\n",
      "   2. group/select target observations in report_map.txt, repeat 1 & 2 as needed\n",
      "   3. run AAVSO() to make report, and\n",
      "   4. submit report to AAVSO.")
  return(df_transformed)
}

markupReport <- function(AN_top_folder="J:/Astro/Images/C14", AN_rel_folder=NULL) {
  if (is.null(AN_rel_folder)) {stop(">>>>> You must provide a AN_rel_folder, ",
                                    "e.g., AN_rel_folder='20151216'.")}
  source("C:/Dev/Photometry/$Utility.R")
  AN_folder   <- make_safe_path(AN_top_folder, AN_rel_folder)
  photometry_folder <- make_safe_path(AN_folder, "Photometry")
  path_df_predictions <- make_safe_path(photometry_folder, "df_transformed.Rdata")
  load(path_df_predictions)
  
  df_checkStars <- df_predictions %>% filter(StarType=="Check") %>% select(FITSfile, StarID)
  JD_fract <- (df_predictions$JD_num - (df_predictions$JD_num %>% min() %>% floor())) %>% round(digits=4)
  
  df_markupReport <- df_predictions %>%
    mutate(JD_fract=JD_fract) %>%
    filter(StarType=="Target") %>%
    select(Serial, FITSfile, Target=StarID, Filter, Exposure, InstMagSigma, MagErr, JD_fract) %>%
    left_join(df_checkStars, by="FITSfile") %>%
    rename(Check=StarID) %>%
    mutate(Sigma_mMag=round(InstMagSigma*1000)) %>%
    mutate(Err_mMag=round(MagErr*1000)) %>%
    select(-InstMagSigma, -MagErr) %>%
    arrange(Target, FITSfile, Filter, Exposure) %>%
    select(Serial, Target, FITSfile, Filter, Exposure, everything())
  return(df_markupReport)
}

  
AAVSO <- function (AN_top_folder="J:/Astro/Images/C14", AN_rel_folder=NULL, software_version=NULL) {
  ##### Writes AAVSO-ready text file.
  ##### No return value.
  ##### First testing 20160207.
  require(dplyr, quietly=TRUE)
  if (is.null(AN_rel_folder)) {stop(">>>>> You must provide a AN_rel_folder parm, ",
                                  "e.g., AN_rel_folder='20151216'.")}
  if (is.null(software_version)) {stop(">>>>> You must provide a software_version parm, ",
                                    "e.g., software_version='0.12'.")}
  
  # First, make df_report.
  AN_folder   <- make_safe_path(AN_top_folder, AN_rel_folder)
  photometry_folder <- make_safe_path(AN_folder, "Photometry")
  df_report <- make_df_report(photometry_folder)
  
  # Construct report's header lines.
  out <- "#TYPE=EXTENDED" %>%
    c("#OBSCODE=DERA") %>% #       DERA = Eric Dose's observer code @ AAVSO
    c("#SOFTWARE=custom R Scripts, github/edose") %>%
    c("#DELIM=,") %>%
    c("#DATE=JD") %>%
    c("#OBSTYPE=CCD") %>%
    c(paste0("#This report of ", nrow(df_report), " observations was generated ", 
             format(Sys.time(), tz = "GMT"), " UTC.")) %>%
    c(paste0("#Software version for all computation yielding this report = ", software_version)) %>%
    c("#Eric Dose, Bois d'Arc Observatory, Kansas") %>%
    c("#") %>%
    c("#NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,GROUP,CHART,NOTES")
  
    if (nrow(df_report) >= 1) {
  # Final formatting of all observation report fields.
    df_report <- df_report %>%
      mutate(TargetName = TargetName %>% trimws() %>% toupper()) %>%
      mutate(JD         = JD         %>% as.numeric() %>% round(5) %>% format(nsmall=5)) %>%
      mutate(Mag        = Mag        %>% round(3) %>% format()) %>%
      mutate(MagErr     = MagErr     %>% round(3) %>% format()) %>%
      mutate(Filter     = Filter     %>% trimws() %>% toupper()) %>%
      mutate(CompName   = CompName   %>% trimws() %>% toupper()) %>%
      mutate(CompMag    = ifelse(CompMag=="na","na",CompMag      %>% round(3) %>% format())) %>%
      mutate(CheckName  = ifelse(is.na(CheckName),"na",CheckName %>% trimws() %>% toupper())) %>%
      mutate(CheckMag   = ifelse(is.na(CheckMag), "na",CheckMag  %>% round(3) %>% format())) %>%
      mutate(Airmass    = Airmass    %>% round(4) %>% format()) %>%
      mutate(Chart      = Chart      %>% trimws() %>% toupper()) %>%
      mutate(Notes      = ifelse(trimws(Notes)=="","na",trimws(Notes)))
  
    # Append rows of df_report to character vector "out".
    obs_lines <- paste(
      df_report$TargetName,
      df_report$JD,
      df_report$Mag,
      df_report$MagErr,
      df_report$Filter,
      "YES",                   # always Transformed for full model.
      "STD",                   # we use standard comp stars, not "differential" mode.
      df_report$CompName,      # ="ENSEMBLE" when more than one comp star for this observation.
      df_report$CompMag,       # "na" if Ensemble comp, else instrument comp mag (or maybe "na" even so).
      df_report$CheckName,
      df_report$CheckMag,
      df_report$Airmass,
      "na",                   # GROUP, not used here.
      df_report$Chart,
      df_report$Notes,
      sep=",")
    out <- out %>% c(obs_lines)
  } else {
    out <- out %>% c("\n\n\n                    ########## NO OBSERVATIONS TO PRINT ##########\n")
  }
  
  # Now dump the char vector "out" to a text file.
  filename_report <- paste0("AAVSOreport-", AN_rel_folder, ".txt")
  path_report   <- make_safe_path(photometry_folder, filename_report)
  write(out,file=path_report)
  cat(paste0("AAVSO report for AN ",AN_rel_folder," written to: ", path_report, "\n   = ",
    nrow(df_report)," observations."))
}


################################################################################################
##### Below are support-only functions, not called by user. ####################################

make_df_report <- function(photometry_folder) {
  ##### Inputs: df_transformed, masterModelList.
  ##### Returns: df_report for AAVSO().
  ##### In development EVD Feb 6, 2016.
  require(dplyr, quietly=TRUE)
  require(stringi, quietly=TRUE)
  
  path_df_transformed <- make_safe_path(photometry_folder, "df_transformed.Rdata")
  load(path_df_transformed)
  path_masterModelList <- make_safe_path(photometry_folder, "masterModelList.Rdata")
  load(path_masterModelList)
  
  df_report <- df_transformed %>%
    select(Serial, StarType, TargetName=StarID, JD=JD_mid, Mag=TransformedMag, MagErr, Filter) %>%
    mutate(CompName=NA, CompMag=NA, nComps=NA, CheckName=NA, CheckMag=NA) %>%
    mutate(Airmass=df_transformed$Airmass, Chart=df_transformed$Chart) %>%
    mutate(Notes=paste0("Serial=",Serial)) %>%
    filter(StarType=="Target") %>%   # do this late to keep df_report & df_transformed aligned.
    select(-StarType)
           
  # Get report_map.txt file for this AN folder, then parse directive lines.
  report_map_file <- make_safe_path(photometry_folder, "report_map.txt")
  if (!file.exists(report_map_file)) {
    cat("Report Map file", report_map_file, "does not exist.")
    return (NA)
  } else {
    lines <- readLines(report_map_file, warn=FALSE)   # read last line even without EOL character(s).
    for (iLine in 1:length(lines)) {
      lines[iLine] <- lines[iLine] %>% 
        strsplit(";",fixed=TRUE) %>% unlist() %>% first() %>% trimws()  # remove comments
    }
    directiveLines <- lines[stri_detect_regex(lines,'^#')] 
    directiveLines <- directiveLines[!is.na(directiveLines)]
    # detect and collect directive text lines.
  }
  
  ##### Apply report_map.txt omissions before any other manipulations.
  # Apply #TARGET directives.
  omit_lines <- directiveLines[stri_startswith_fixed(directiveLines, "#TARGET")]
  for (thisLine in omit_lines) {
    targetToOmit <- substring(thisLine, first=nchar("#TARGET")+1) %>% trimws()
    df_report <- df_report %>% filter(toupper(TargetName)!=toupper(targetToOmit))
    cat(paste0("Omitted TARGET: ", targetToOmit), "\n")
  }
  # Apply #SERIAL directives.
  omit_lines <- directiveLines[stri_startswith_fixed(directiveLines, "#SERIAL")]
  serialsToOmit <- as.numeric() # we may want to keep this for error msgs, below.
  for (thisLine in omit_lines) {
    words <- stri_extract_all_words(thisLine) %>% unlist() %>% strsplit(",") %>% unlist()
    serials <- words[-1] %>% as.numeric()
    serials <- serials[serials>=1]
    serialsToOmit <- c(serialsToOmit, serials)
  }
  df_report <- df_report %>% filter(!Serial %in% serialsToOmit)
  cat(paste0("Omitting ", length(serialsToOmit), " observations by Serials ",
             paste0(serialsToOmit, collapse=" "), " done.\n"))
  
  # Apply check-star names and mags (lookup from df_transformed, then insert column in-place).
  df_checks <- df_transformed %>% filter(StarType=="Check") %>% select(JD_mid, StarID, TransformedMag)
  df_joined <- df_report      %>% left_join(df_checks, by=c("JD"="JD_mid"))
  df_report <- df_report      %>% mutate(CheckName=df_joined$StarID, CheckMag=df_joined$TransformedMag)
  
  ##### Apply comp-star names and mags.
  # Make df_comps by rbind() all $obs from masterModelList
  df_comps <- data.frame(stringsAsFactors=FALSE)
  for (thisModelList in masterModelList) {
    df_thisComp <- thisModelList$obs %>% filter(StarType=="Comp") %>% select(Serial,JD_mid,StarID,Filter)
    df_comps <- df_comps %>% rbind(df_thisComp)
  }
  # Now insert comp info.
  reportJDs <- (df_report %>% arrange(JD))$JD %>% unique()
  for (thisJD in reportJDs) {
    df_compsThisJD <- df_comps %>% filter(JD_mid==thisJD)
    nCompsThisJD <- nrow(df_compsThisJD)
    df_report$nComps[df_report$JD==thisJD] <- nCompsThisJD
    if (nCompsThisJD<=0) {
      stop(">>>>> No comp stars in model for JD = ", thisJD)
    }
    if(nCompsThisJD > 1) { 
      # the ensemble case.
      df_report$CompName[df_report$JD==thisJD] <- "ENSEMBLE"
      df_report$CompMag[df_report$JD==thisJD]  <- "na"
      df_report$Notes[df_report$JD==thisJD]    <- paste0(df_report$Notes[df_report$JD==thisJD], 
                                                         " (", nCompsThisJD, " comps)")
    } else { 
      # the single-comp-star case.
      df_report$CompName[df_report$JD==thisJD] <- df_compsThisJD$StarID[1]
      df_report$CompMag[df$reportJD==thisJD]  <- "na" # unsure what AAVSO's "Inst Mag of Comp Star" means.
    }
  }

  ###### Lastly, apply report_map.txt #COMBINE directives (careful with the check and comp stars associated).
  allSame <- function (vect) {  # nested function
    return(all(vect==vect[1]))
  }
  combine_lines <- directiveLines[stri_startswith_fixed(directiveLines, "#COMBINE")]
  for (thisLine in combine_lines) {
    words <- stri_extract_all_words(thisLine) %>% unlist() %>% strsplit(",") %>% unlist()
    serials <- words[-1] %>% as.numeric()
    serialsToCombine <- serials[serials>=1]
    df_combine <- df_report %>% filter(Serial %in% serialsToCombine)
    
    # Verify that user-selected combine rows are truly eligible to be combined.
    if (nrow(df_combine)<=1) {
      cat(paste0(">>>>> No combines for line: '", thisLine, "'.\n"))
      next
    }
    if (!allSame(df_combine$TargetName)){
      cat(paste0(">>>>> Unlike Target Names for line: '", thisLine, "'...Combine is skipped.\n"))
      next
    }
    if (!allSame(df_combine$Filter)){
      cat(paste0(">>>>> Unlike Filters for line: '", thisLine, "'...Combine is skipped.\n"))
      next
    }
    if (!allSame(df_combine$CompName)){
      cat(paste0(">>>>> Unlike Comp Names for line: '", thisLine, "'...Combine is skipped.\n"))
      next
    }
    if (!allSame(df_combine$CheckName)){
      cat(paste0(">>>>> Unlike Check Stars for line: '", thisLine, "'...Combine is skipped.\n"))
      next
    }
    if (!allSame(df_combine$Chart)){
      cat(paste0(">>>>> Unlike Chart IDs for line: '", thisLine, "'...Combine is skipped.\n"))
      next
    }
    if (abs(diff(range(as.numeric(df_combine$JD)))) > 1/24){ # greater than 1 hour
      cat(paste0(">>>>> Range of JD times is too large for line: '", thisLine, "'...Combine is skipped.\n"))
      next
    }
    if (abs(diff(range(df_combine$Airmass))) > 0.100){
      cat(paste0(">>>>> Range of Airmasses is too large for line: '", thisLine, 
                 "'...Combine is skipped.\n"))
      next
    }
    
    # Now combine the rows (observations) into the first selected row, delete other rows.
    df_new <- df_combine[1,]
    serialToReplace  <- df_combine$Serial[df_combine$Serial==df_new$Serial]
    serialsToDelete  <- df_combine$Serial[df_combine$Serial!=serialToReplace]
    # Column Serial is simply the first observation's Serial.
    # Columns TargetName, Filter, CompName, CompMag, CheckStar, and Chart: uniform and inherited.
    df_new$JD       <- mean(as.numeric(df_combine$JD))
    df_new$Mag      <- mean(df_combine$Mag)
    df_new$MagErr   <- max( min(df_combine$MagErr), max(df_combine$MagErr/sqrt(nrow(df_combine))))
    df_new$CheckMag <- mean(df_combine$CheckMag)
    df_new$Airmass  <- mean(df_combine$Airmass)
    df_new$Notes    <- paste0("Serials=[", (paste0(df_combine$Serial,collapse=" ")), "] (",
                              min(df_combine$nComps), "-", max(df_combine$nComps), " comps)")
    df_report[df_report$Serial==serialToReplace,] <- df_new[1,]
    df_report <- df_report %>% filter(!Serial %in% serialsToDelete)
    cat(paste0("Combination of Serials ", paste0(df_combine$Serial, collapse=" "), " done.\n"))
  }

  return(df_report) # return without final formatting for AAVSO report (leave that to AAVSO() function).
}

imputeTargetCIs <- function (df_predictions, CI_filters, transforms) {
  # Impute Color Index CI (=true V Mag - true I Mag) from untransformed mags.
  require(dplyr, quietly=TRUE)
  JD_floor <- floor(min(as.numeric(df_predictions$JD_mid)))
  df_predictions <- df_predictions %>% mutate(JD_num=as.numeric(JD_mid)-JD_floor) # avoid scaling problems.
  
  targetStarIDs <- (df_predictions %>%
    filter(StarType=="Target"))$ModelStarID %>%
    unique()
  
  for (thisTargetStarID in targetStarIDs) { # handle one target's ModelStarID at a time.
    df_targetStarID <- df_predictions %>% 
      filter(ModelStarID==thisTargetStarID) %>%
      select(Serial, ModelStarID, Filter, JD_num, CI, UntransformedMag) %>%
      arrange(JD_num)
    df_CI_points <- extractCI_points(df_targetStarID, CI_filters, transforms) %>% 
      arrange(JD_num)

    # Interpolate CI and put it in df_targetStarID
    i <- 0
    if (nrow(df_CI_points) == 0) {
      df_targetStarID$CI <- NA  # because there are no color index values to apply.
      cat(">>>>> ModelStarID=", thisTargetStarID, ": no CI points returned by imputeTargetCIs()\n", sep="")
    }
    if (nrow(df_CI_points) == 1) {
      df_targetStarID$CI <- df_CI_points$CI[1] # set all CI to the same value
    }
    if (nrow(df_CI_points) %in% 2:3) { # linear interpolation (with residuals if 3 points)
      m <- lm (CI ~ JD_num, data=df_CI_points)
      df_targetStarID$CI <- predict.lm(m, data.frame(JD_num=df_targetStarID$JD_num))
      # Enforce no extrapolation.
      df_targetStarID$CI[df_targetStarID$JD_num < df_CI_points$JD_num[1]] <- 
        predict.lm(m,data.frame(JD_num=df_CI_points$JD_num[1]))
      df_targetStarID$CI[df_targetStarID$JD_num > df_CI_points$JD_num[nrow(df_CI_points)]] <-
        predict.lm(m,data.frame(JD_num=df_CI_points$JD_num[nrow(df_CI_points)]))
    }
    if (nrow(df_CI_points) >= 4) { # make and apply smoothing spline (std package "stats").
      degrees_freedom <- round(nrow(df_CI_points)/6) %>% max(4) %>% min(nrow(df_CI_points))
      thisSpline <- smooth.spline(df_CI_points$JD_num, df_CI_points$CI, df=degrees_freedom)
      df_targetStarID$CI <- predict(thisSpline, df_targetStarID$JD_num)$y
      # Enforce no extrapolation.
      df_targetStarID$CI[df_targetStarID$JD_num < df_CI_points$JD_num[1]] <- 
        predict(thisSpline, df_CI_points$JD_num[1])$y
      df_targetStarID$CI[df_targetStarID$JD_num > df_CI_points$JD_num[nrow(df_CI_points)]] <-
        predict(thisSpline, df_CI_points$JD_num[nrow(df_CI_points)])$y
    }

    # Do the CI value replacements.
    df_predictions$CI[match(df_targetStarID$Serial, df_predictions$Serial)] <- df_targetStarID$CI
    i <- 0
  }
  return (df_predictions) # which now includes filled-in CI values and new column JD_num.
}

extractCI_points <- function (df, CI_filters, transforms) {
  ##### df must be a data frame from df_predictions, subset holding one ModelStarID.
  require(dplyr, quietly=TRUE)
  maxDiffJD <- 15 / (24*60) # max time between paired obs, in (Julian) days, typically 15 minutes.
  df <- df %>% filter(Filter %in% CI_filters) %>% arrange(JD_num)
  
  df_CI_point <- data.frame()
  if (nrow(df)>=2) {
    for (iStart in 1:(nrow(df)-1)) {
      if (df$Filter[iStart+1] != df$Filter[iStart]) {
        if (df$JD_num[iStart+1]-df$JD_num[iStart] <= maxDiffJD) {
          this_JD_num <- (df$JD_num[iStart+1] + df$JD_num[iStart]) / 2
          # CI_predicted is color index from naive, predicted, untransformed (not catalog-basis) mags.
          CI_predicted <- ifelse(df$Filter[iStart]==CI_filters[1],
                                 df$UntransformedMag[iStart]-df$UntransformedMag[iStart+1],
                                 df$UntransformedMag[iStart+1]-df$UntransformedMag[iStart])
          this_CI <- CI_predicted / (1 + transforms[1] - transforms[2]) # transforms -> catalog-basis CI.
          df_CI_point <- df_CI_point %>% rbind(data.frame(JD_num=this_JD_num, CI=this_CI))
        } 
      }
    }
} 
  return (df_CI_point)
}


