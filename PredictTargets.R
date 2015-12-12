##### PredictTargets.R   Predict magnitudes of Target (unknown) stars in Astronight's master data frame.
#####    Uses this Astronight's lmer() model object etc from Model.R::modelAll().

predictAll <- function (masterModelList, df_master, AN_folder, maxMagUncertainty=0.05) {
  require(dplyr)
  filtersModelList <- names(masterModelList)
  filtersDfMaster  <- df_master$Filter %>% unique()
  filters <- intersect(filtersModellist, filtersDfMaster)
  
  df_targets <- data.frame()
  for (filter in filters) {
    df_targets <- df_targets %>%
      rbind(predictOneFilter(modelList, df_master, filter, maxMagUncertainty))
  }
  return (df_targets)
}


predictOneFilter <- function (modelList, df_master, filter, maxMagUncertainty) {
  require(dplyr)
  # Unpack input model list.
  model      <- modelList$model
  filter     <- modelList$filter
  transform  <- modelList$transform
  extinction <- modelList$extinction
  
  # Make working data frame and perform raw prediction.
  df <- df_master %>%
    filter(StarType %in% c("Check","Target")) %>%
    filter(Filter==filter) %>%
    filter(UseInModel==TRUE) %>%
    filter(Saturated==FALSE) %>%
    filter(MagUncertainty<=maxMagUncertainty) %>%
    mutate(CI=ifelse(is.na(CI),0,CI)) %>%     # zero Targets' (but not Checks') color index (adjust later).
    mutate(CatMag=0, PredictedMag=NA, estimError=NA)
  modelMag <- predict(model,df,re.form=~(1|JD_mid)) # re.form to omit modelStarID from formula.
  
  # Adjust for terms not included in model's formula.
  if (!is.na(transform)) {
    modelMag <- modelMag + transform * df$CI
  }
  if (!is.na(extinction)) {
    modelMag <- modelMag + extinction * Airmass
  }

  # Compute and store predicted magnitude.
  df <- df %>% mutate(PredictedMag = df$InstMag - modelMag)
  
  ##### TODO : compute & store estimated errors.

  return (df)
}


write_AAVSO_reportfile <- function (AN_folder) {
  #####    Writes AAVSO-ready text file.
  require(dplyr)
  out <- "#TYPE=EXTENDED" %>%
    c("#OBSCODE=DERA") %>% #       DERA = Eric Dose's observer code @ AAVSO
    c("#SOFTWARE=custom R Scripts, github/edose") %>%
    c("#DELIM=,") %>%
    c("#DATE=JD") %>%
    c("#OBSTYPE=CCD") %>%
    c("#NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,GROUP,CHART,NOTES")

  # TODO : (1) read in all df_target_V etc data frames, (2) create df_report.
  
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
