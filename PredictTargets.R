##### PredictTargets.R   Predict magnitudes of Target (unknown) stars in Astronight's master data frame.
#####    Writes AAVSO-ready text file.
#####    Uses this Astronight's lmer() model object from Model.R::model().

write_AAVSO_reportfile <- function (df_master, AN_folder) {
  require(dplyr)
  out <- "#TYPE=EXTENDED"
  out <- out %>% c("#OBSCODE=DERA") # DERA = Eric Dose
  out <- out %>% c("#SOFTWARE=custom R Scripts")
  out <- out %>% c("#DELIM=,")
  out <- out %>% c("#DATE=JD")
  out <- out %>% c("#OBSTYPE=CCD")
  out <- out %>% c("#NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,GROUP,CHART,NOTES")
  
  # TODO : curate df_report from df_in.
  
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

