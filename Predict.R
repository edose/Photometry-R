##### PredictTargets.R   Predict magnitudes of Target (unknown) stars in Astronight's master data frame.
#####    Uses this Astronight's lmer() model object etc from Model.R::modelAll().

predictAll <- function (AN_top_folder="J:/Astro/Images/C14", AN_rel_folder=NA, maxMagUncertainty=0.05) {
  ##### Calls predictOneFilter for each modelList in the masterModelList stored in AN's Photometry folder.
  ##### Returns big data frame df_predict (which is ready for transformation by applyTransforms()).
  ##### Typical usage: df_predict <- predictAll(AN_rel_folder="20151216")
  if (is.na(AN_rel_folder)) {stop(">>>>> You must provide a AN_rel_folder parm, ",
                                  "e.g., AN_rel_folder='20151216'.")}
  require(dplyr, quietly=TRUE)
  source("C:/Dev/Photometry/$Utility.R")
  AN_folder   <- make_safe_path(AN_top_folder, AN_rel_folder)
  photometry_folder <- make_safe_path(AN_folder, "Photometry")
  df_master_path <- make_safe_path(photometry_folder, "df_master.Rdata")
  load(df_master_path, verbose=TRUE)
  masterModelList_path <- make_safe_path(photometry_folder, "masterModelList.Rdata")
  load(masterModelList_path, verbose=TRUE)
  
  df_predict <- data.frame()
  for (modelList in masterModelList) {
    df_predict <- df_predict %>%
      rbind(predictOneFilter(modelList=modelList, df_master=df_master, maxMagUncertainty=maxMagUncertainty))
  }
  return (df_targets)
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
  if (!is.na(transform)) {
    modelMag <- modelMag + transform * df$CI
  }
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


