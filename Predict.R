##### PredictTargets.R   Predict magnitudes of Target (unknown) stars in Astronight's master data frame.
#####    Uses this Astronight's lmer() model object etc from Model.R::modelAll().

predictAll <- function (AN_top_folder="J:/Astro/Images/C14", AN_rel_folder=NULL, 
                        saturatedADU=54000, maxMagUncertainty=0.05, CI_filters=c("V","I")) {
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
  load(path_df_master, verbose=TRUE)
  path_masterModelList <- make_safe_path(photometry_folder, "masterModelList.Rdata")
  load(path_masterModelList, verbose=TRUE)
  
  # Select Target and Check rows for which CatMag to be predicted (all filters).
  # Make df of all needed inputs (newdata) for all filters & all (untransformed) lme4::predict() calls.
  require(dplyr, quietly=TRUE)
  df_predict_input <- df_master %>%
    filter(StarType %in% c("Check","Target")) %>%
    filter(UseInModel==TRUE) %>%
    filter(MaxADU<=saturatedADU) %>%
    filter(MagUncertainty<=maxMagUncertainty) %>%
    mutate(CI=ifelse(is.na(CI),0,CI)) %>%
    select(Serial, ModelStarID, Xpixels, Ypixels, InstMag, MagUncertainty, StarType,
           JD_mid, Filter, Airmass, CI, Vignette, Vignette4) %>%
    mutate(CatMag=0) # arbitrarily chosen, to complete model.

  # Get untransformed predicted Inst Mags (via running predict(), collecting all results).
  df_predictions <- data.frame()
  require(lme4, quietly = TRUE)
  for (thisModelList in masterModelList) {
    df_thisFilter <- df_predict_input %>% filter(Filter==thisModelList$filter)
    predictOutput <- predict(thisModelList$model, newdata=df_thisFilter, re.form=~(1|JD_mid))
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
  df_transformed <- df_predictions %>%
    left_join(df_xref_transform, by="Filter") %>%
    mutate(TransformedMag = UntransformedMag - Transform * CI) # minus is the correct sign: 20160206.

  # Compute MagErr (for AAVSO) as greater of observation's MagUncertainty and model's Sigma.
  df_xref_modelSigma <- data.frame(stringsAsFactors=FALSE) # build a lookup table to get model mag errs
  for (thisModelList in masterModelList) {
    df_xref_modelSigma <- df_xref_modelSigma %>%
      rbind(data.frame(Filter=thisModelList$filter, ModelSigma=sigma(thisModelList$model)))
  }
  df_transformed <- df_transformed %>%
    left_join(df_xref_modelSigma, by="Filter") %>%
    mutate(MagErr = pmax(ModelSigma, MagUncertainty)) # pmax = "parallel" element-wise max of 2 vectors
  
  df_transformed <- df_transformed %>% arrange(ModelStarID, JD_num)
  
  # Save df_transformed as .Rdata file (but do not return the data frame).
  path_transformed <- make_safe_path(photometry_folder, "df_transformed.Rdata")
  save(df_transformed, file=path_transformed, precheck=FALSE)
  
  # Make a template-only report_map.txt file if report_map.txt doesn't already exist.
  path_report_map <- make_safe_path(photometry_folder, "report_map.txt")
  if (!file.exists(path_report_map)) {
    lines <- c(
      paste0(";----- This is report_map.txt for AN folder ", AN_rel_folder),
      paste0(";----- Use this file to combine and/or omit target observations from AAVSO report."),
      paste0(";----- Example directive lines:\n"),
      paste0(";"),
      paste0(";"),
      paste0(";"),
      paste0(";"),
      paste0(";\n;----- Add your directive lines:\n;\n\n")
    )
    writeLines(lines, con=path_report_map)
  }
  nTargetObs <- df_transformed %>% filter(StarType=="Target") %>% nrow()
  nTargets   <- df_transformed %>% filter(StarType=="Target") %>% select(ModelStarID) %>% 
    unique() %>% nrow()
  nCheckObs  <- df_transformed %>% filter(StarType=="Check") %>% nrow()
  Filters    <- df_transformed %>% filter(StarType=="Target") %>% select(Filter) %>% 
    unique()
  cat("PredictAll() yields", nTargetObs, "Target obs for", nTargets, "in filters:", 
      paste(Filters,collapse=" ","\n"))
  cat("   and ", nCheckObs, " Check observations.\n")

  cat("predictAll() has saved this AN's transformed predictions to", path_transformed, "\n",
      "   and has ensured a report_map.txt is available in the same folder\n",
      "Now you are ready to:\n",
      "   1. run targetPlots(),\n",
      "   2. group/select target observations in report_map.txt, repeat 1 & 2 as needed\n",
      "   3. run AAVSO() to make report, and\n",
      "   4. submit report to AAVSO.")
}
  
AAVSO <- function (AN_top_folder="J:/Astro/Images/C14", AN_rel_folder=NULL) {
  #####    Writes AAVSO-ready text file.
  require(dplyr, quietly=TRUE)
  if (is.null(AN_rel_folder)) {stop(">>>>> You must provide a AN_rel_folder parm, ",
                                  "e.g., AN_rel_folder='20151216'.")}
  
  out <- "#TYPE=EXTENDED" %>%
    c("#OBSCODE=DERA") %>% #       DERA = Eric Dose's observer code @ AAVSO
    c("#SOFTWARE=custom R Scripts, github/edose") %>%
    c("#DELIM=,") %>%
    c("#DATE=JD") %>%
    c("#OBSTYPE=CCD") %>%
    c("#NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,GROUP,CHART,NOTES")

  # Read in all required data (should only need df_transformed.Rdata and report_map.txt).
  AN_folder   <- make_safe_path(AN_top_folder, AN_rel_folder)
  photometry_folder <- make_safe_path(AN_folder, "Photometry")
  df_report <- make_df_report(photometry_folder)
  
  # Format rows of df_report as a character vector and append them to character vector "out".
  obs_lines <- paste(
    df_report$target,
    df_report$JD_mid,
    df_report$mag,
    df_report$mag_error,
    df_report$filter,
    df_report$isTransformed, # always TRUE for full model.
    "STD",                   # we use standard comp stars, not "differential" mode.
    df_report$CNAME,         # ="ENSEMBLE" when more than one comp star for this observation.
    "na",                    # CMAG
    df_report$check_star,
    df_report$check_mag,
    df_report$airmass,
    "na",                   # GROUP, not used here.
    df_report$chart,
    df_report$notes,        # usually "na"
    sep=",")
  out <- out %>% c(obs_lines)

  # Now dump the char vector "out" to a text file.
  out_folder <- make_safe_path(AN_folder,"Photometry")
  out_path   <- make_safe_path(out_folder,"AAVSOreport.txt")
  print(out_path)
  write(out,file=out_path)
  print(paste("AAVSO report for folder ",AN_folder," written: ",length(obs_lines)," observations."))
}


################################################################################################
##### Below are support-only functions, not called by user. ####################################

make_df_report <- function() {
  
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
  return (df_CI_point)
}


