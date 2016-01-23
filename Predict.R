##### PredictTargets.R   Predict magnitudes of Target (unknown) stars in Astronight's master data frame.
#####    Uses this Astronight's lmer() model object etc from Model.R::modelAll().

predictAll <- function (AN_top_folder="J:/Astro/Images/C14", AN_rel_folder=NA, 
                        maxMagUncertainty=0.05, CI_filters=c("V","I")) {
  ##### Performs ALL steps to get transformed magnitudes for Check and Target observations in df_master.
  ##### Returns big data frame df_predict (which is transformed).
  ##### Requires that AN folder contain R files df_master.Rdata (from Input.R::make_df_master()) and 
  #####    masterModelList.Rdata (from Model.R::make_masterModelList()).
  ##### Typical usage: df_predict <- predictAll(AN_rel_folder="20151216")
  if (is.na(AN_rel_folder)) {stop(">>>>> You must provide a AN_rel_folder, ",
                                  "e.g., AN_rel_folder='20151216'.")}
  source("C:/Dev/Photometry/$Utility.R")
  AN_folder   <- make_safe_path(AN_top_folder, AN_rel_folder)
  photometry_folder <- make_safe_path(AN_folder, "Photometry")
  df_master_path <- make_safe_path(photometry_folder, "df_master.Rdata")
  load(df_master_path, verbose=TRUE)
  masterModelList_path <- make_safe_path(photometry_folder, "masterModelList.Rdata")
  load(masterModelList_path, verbose=TRUE)
  
  # Select Target and Check rows for which CatMag to be predicted (all filters).
  require(dplyr, quietly=TRUE)
  df_targets <- df_master %>%
    filter(StarType %in% c("Check","Target")) %>%
    filter(UseInModel==TRUE) %>%
    filter(MaxADU<=saturatedADU) %>%
    filter(MagUncertainty<=maxMagUncertainty) %>%
    mutate(CI=ifelse(is.na(CI),0,CI))
  
  # Make df of all needed inputs (newdata) for all filters & all (untransformed) lme4::predict() calls.
  df_predict_input <- df_targets %>%
    select(Serial, ModelStarID, XPixels, YPixels, InstMag, MagUncertainty, StarType,
           JD_mid, Filter, Airmass, CI, Vignette, Vignette4) %>%
    mutate(CatMag=0) # arbitrarily chosen, to complete model.

  # Run predict() & collect all results, which are untransformed predicted instrument magnitudes.
  df_predictions <- data.frame()
  require(lme4, quietly = TRUE)
  for (thisModelList in masterModelList) {
    df_newdata <- df_predict_input %>%
      filter(Filter==thisModelList$filter)
    predictOutput <- predict(thisModelList$model, newdata=df_newdata, re.form=~(1|JD_mid))
    # Correct for extinction if extinction was not fit ("Airmass") in model.
    if (!thisModelList$fit_extinction) {
      predictOutput <- predictOutput + extinction * df_newdata$Airmass
    }
    this_df_predictions <- df_newdata %>%
      mutate(PredictedInstMag <- predictOutput)
    df_predictions <- rbind(df_predictions, this_df_predictions)
  }

  # Get UntransformedMags from PredictedInstMag, InstMag (and CatMag=0).
  df_predictions <- df_predictions %>%
    mutate(UntransformedMag = InstMag - PredictedInstMag)
  
  # Impute Color Index values for target observations by interpolation.
  df_predictions <- imputeTargetCIs(df_predictions, CI_filters) # Fill in CI values for target obs.
  
  # Perform transformations.
  df_xref_transform <- data.frame() # build a lookup table to fill in Transforms for each filter.
  for (thisModelList in masterModelList) {
    df_xref_transform <- df_xref_transform %>%
      rbind(data.frame(Filter=modelList$filter, Transform=modelList$transform))
  }
  df_predictions <- df_predictions %>%
    left_join(df_xref_transform) %>%
    mutate(TransformedMag = UntransformedMag - Transform * CI) # unsure of sign, here.
  
  # Save df_predictions as .Rdata file (but do not return the data frame).
  predictions_path <- make_safe_path(photometry_folder, "df_predictions.Rdata")
  save(df_predictions, file=predictions_path, precheck=FALSE)
  cat("predictAll() has saved this AN's predictions to", predictions_path, 
      "\n   Now ready to:\n   (1) predictionPlots(),\n   (2) group/select target observations to report,\n",
      "   (3) run AAVSO() to make report,\n   and (4) submit report to AAVSO.")
}

predictOneFilter <- function (filterModelList, df_master, saturatedADU=54000, maxMagUncertainty) {
  require(dplyr, quietly=TRUE)
  # Unpack input model list.
  model      <- filterModelList$model
  transform  <- filterModelList$transform
  extinction <- filterModelList$extinction
  JD_mid_model_levels <- model.frame(model) %>% select(JD_mid) %>% unique() %>% unlist() %>% levels()
  
  # Make working data frame and perform raw prediction.
  df <- df_master %>%
    filter(StarType %in% c("Check","Target")) %>%
    filter(Filter==filter) %>%
    filter(UseInModel==TRUE) %>%
    filter(MaxADU<=saturatedADU) %>%
    filter(MagUncertainty<=maxMagUncertainty) %>%
    filter(JD_mid %in% JD_mid_model_levels) %>% # remove images not represented in this filter's model.
    mutate(CI=ifelse(is.na(CI),0,CI)) %>%     # zero Targets' (but not Checks') color index (adjust later).
    mutate(CatMagSaved=CatMag) %>%
    mutate(CatMag=0, PredictedMag=NA, estimError=NA) # CatMag must be zero to do prediction properly.
  modelMag <- predict(model, df, re.form=~(1|JD_mid)) # re.form to omit modelStarID from formula.
  df <- df %>%
    mutate(CatMag=CatMagSaved) %>%   # restore CatMag (which needed to be zero to do predictions).
    select(-CatMagSaved)            
  
  # Adjust for terms not included in model's formula.
  if (!is.na(extinction)) {
    modelMag <- modelMag + extinction * df$Airmass
  }

  # Compute and store predicted magnitude.
  df <- df %>% mutate(PredictedMag = df$InstMag - modelMag)
  
  ##### Compute & store estimated errors.

  return (df)
}

writeAAVSO <- function (AN_top_folder="J:/Astro/Images/C14", AN_rel_folder=NA) {
  #####    Writes AAVSO-ready text file.
  require(dplyr, quietly=TRUE)
  if (is.na(AN_rel_folder)) {stop(">>>>> You must provide a AN_rel_folder parm, ",
                                  "e.g., AN_rel_folder='20151216'.")}
  
  out <- "#TYPE=EXTENDED" %>%
    c("#OBSCODE=DERA") %>% #       DERA = Eric Dose's observer code @ AAVSO
    c("#SOFTWARE=custom R Scripts, github/edose") %>%
    c("#DELIM=,") %>%
    c("#DATE=JD") %>%
    c("#OBSTYPE=CCD") %>%
    c("#NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,GROUP,CHART,NOTES")

  # Read in all df_target_V etc data frames.
  
  # Create data frame df_report.
  
  # Format each row of df_report as a text line.
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


