##### PredictTargets.R   Predict magnitudes of Target (unknown) stars in Astronight's master data frame.
#####    Uses this Astronight's lmer() model object etc from Model.R::modelAll().
##### UNDER REWORK from JUNE 25 2016, to convert to new Hybrid Model/Ensemble approach as in SAS 2016 talk.
#####
##### Typical workflow:
#####    Ensure masterModelList is ready to go from ListV, etc via Model.R::make_masterModelList().
#####    df_predictions <- predictAll(AN_rel_folder="20151216")
#####    -- if eclipsers: eclipserComps(df=df_predictions, fov="ST Tri", starID="ST Tri", this_filter="V")
#####    -- if curating eclipser comps, re-run predictAll() (no change in call signature)
#####    eclipserPlot(starID="ST Tri") 
#####    df_markupReport <- markupReport(AN_rel_folder="20151216")
#####    -- examine markup report, esp for COMBINES and poor check star agreement.
#####    -- edit report_map.txt as needed for #SERIAL & #COMBINE directives.
#####    AAVSO(AN_rel_folder="", software_version="1.2.1")
#####    -- examine AAVSO report, re-edit report_map.txt if needed, rerun AAVSO().
#####    Submit/upload AAVSOreport-nnnnnnnn.txt to AAVSO; check for proper upload.
#####    Set all /Photometry files to read only (in Windows).

predictAll <- function (AN_top_folder="J:/Astro/Images/C14", AN_rel_folder=NULL, 
                        saturatedADU=54000, maxInstMagSigma=0.05, CI_filters=c("V","I")) {
  ##### Performs ALL steps to get transformed magnitudes for Check and Target observations in df_master.
  ##### Returns big data frame df_predict (which is transformed).
  ##### Requires that AN folder contain R files df_master.Rdata (from Input.R::make_df_master()) and 
  #####    masterModelList.Rdata (from Model.R::make_masterModelList()).
  ##### Uses masterModelList and df_master to predict, applying:
  #####    CatMag = 0, Color Index=0, and cirrus effect (random effect "JD_mid") = 0
  #####    to the models to get very raw predictions, which then are treated to restore the above 3 effects.
  #####
  ##### Typical usage: df_predictions <- predictAll(AN_rel_folder="20151216")
  if (is.null(AN_rel_folder)) {stop(">>>>> You must provide a AN_rel_folder, ",
                                  "e.g., AN_rel_folder='20151216'.")}
  source("C:/Dev/Photometry/$Utility.R")
  require(dplyr, quietly=TRUE)
  AN_folder   <- make_safe_path(AN_top_folder, AN_rel_folder)
  photometry_folder <- make_safe_path(AN_folder, "Photometry")
  path_df_master <- make_safe_path(photometry_folder, "df_master.Rdata")
  load(path_df_master)
  source("C:/Dev/Photometry/Model.R")
  df_filtered <- omitObs(AN_top_folder=AN_top_folder, AN_rel_folder=AN_rel_folder)
  df_filtered <- curateEclipserComps(AN_top_folder=AN_top_folder, AN_rel_folder = AN_rel_folder, 
                                     df_filtered_master = df_filtered)
  path_masterModelList <- make_safe_path(photometry_folder, "masterModelList.Rdata")
  load(path_masterModelList)
  
  # COMP STARS: Get raw, predicted mag estimates for ALL eligible Comp Star observations in all filters.
  df_input_comps <- df_filtered %>%
    filter(StarType == "Comp") %>%
    filter(MaxADU_Ur<=saturatedADU) %>%
    filter(InstMagSigma<=maxInstMagSigma) %>%
    filter(!is.na(CatMag)) %>%
    filter(!is.na(CI)) %>%  
    mutate(CatMagSaved=CatMag) %>%  # because we'll need this after running predict().
    mutate(CatMag=0)                # zero is required in predict(), which predicts instrum mags.
  images_input_comps <- df_input_comps$FITSfile %>% 
    unique()
  images_with_targets <- (df_filtered %>% 
    filter(StarType=="Target"))$FITSfile %>% 
    unique() %>% 
    intersect(images_input_comps)  # these are the images we need to handle.
  df_estimates_comps <- data.frame(stringsAsFactors = FALSE)
  for (thisModelList in masterModelList) {
    df_predict_input <- df_input_comps %>%
      filter(Filter == thisModelList$filter) %>%
      filter(FITSfile %in% images_with_targets)
    predictOutput <- predict(thisModelList$model, newdata=df_predict_input, re.form= ~0) # a vector
    if (!thisModelList$fit_extinction) {
      predictOutput <- predictOutput + thisModelList$extinction * df_predict_input$Airmass
    }
    ######################################
    if (!thisModelList$fit_transform) {
      predictOutput <- predictOutput + thisModelList$transform * df_predict_input$CI
    }
    #####################################
    columns_post_predict <- c("Serial", "ModelStarID", "FITSfile", "StarID", "Chart", 
                             "Xcentroid", "Ycentroid", "InstMag", "InstMagSigma", "StarType", 
                             "CatMag", "CatMagSaved", "CatMagError", "Exposure",
                             "JD_mid", "Filter", "Airmass", "CI", "SkyBias", "Vignette", "LogADU")
    df_estimates_thisFilter <- df_predict_input %>% 
        select(one_of(columns_post_predict)) %>%
        mutate(EstimatedMag = InstMag - predictOutput)
    df_estimates_comps <- rbind(df_estimates_comps, df_estimates_thisFilter)
  }  
  df_estimates_comps <- df_estimates_comps %>%
    mutate(CatMag=CatMagSaved) %>%  # restoring CatMag, now that predict() has been run
    select(-CatMagSaved) %>%        # we no longer need this
    mutate(UseInEnsemble=TRUE)
  df_input_comps <- NULL            # we no longer need this; raise an error if it's used again.
  
  # Here, df_estimates_comps accounts for airmass(extinction) and color (transforms).
  # Now, we will derive per-image cirrus-effect values from them (for all filters together).
  df_cirrus_effect <- data.frame(Image=images_with_targets, CirrusEffect=NA_real_, CirrusSigma=NA_real_, 
                                 Criterion1=NA_real_, Criterion2=NA_real_, 
                                 NumCompsUsed = NA_integer_, CompIDsUsed=NA_character_, 
                                 NumCompsRemoved=NA_integer_,
                                 stringsAsFactors = FALSE)
  for (image in df_cirrus_effect$Image) {
    df_estimates_this_image <- df_estimates_comps %>% filter(FITSfile==image)
    # cat(paste0(image, "  n=", nrow(df_estimates_this_image), "\n"))
    cirrus_effect_per_comp <- df_estimates_this_image$EstimatedMag - df_estimates_this_image$CatMag
    # Computes comp-star weights for mean & sigma.
    sigma2 <- df_estimates_this_image$CatMagError^2 + df_estimates_this_image$InstMagSigma^2
    least_allowed_sigma <- 0.01
    raw_weights <- 1/(pmax(sigma2, least_allowed_sigma^2))
    normalized_weights <- raw_weights / sum(raw_weights)  # normalized to sum to one (thus v1 must = 1).
    v2 <- sum(normalized_weights^2)  # var(wt.mean) / wt.var(points)
    # Calculate cirrus effect and error for this image.
    cirrus_effect_this_image <- weighted.mean(cirrus_effect_per_comp, normalized_weights)
    resid2 <- (cirrus_effect_per_comp - cirrus_effect_this_image)^2
    cirrus_sigma_this_image <- ifelse(nrow(df_estimates_this_image)==1, 
                                      df_estimates_this_image[1,"CatMagError"], 
                                      sqrt(v2*sum(normalized_weights * resid2)))
    compIDs_used <- df_estimates_this_image$StarID   # default; kept if all comps good.
    num_comps_used <- nrow(df_estimates_this_image)  #   "
    num_comps_removed <- 0                           #   "
    # Reject this image's worst comp stars and recalculate, if:
    #   (at least 4 comp stars) AND (criterion1 >= 16  OR  criterion2 >= 20).
    criterion1 <- NA  # default
    criterion2 <- NA  # default
    if (nrow(df_estimates_this_image) >= 4) {
      # criterion 1 = how many times worse is the worst comp vs the average of the other comps.
      x <- normalized_weights * resid2
      criterion1 <- max(x) / ((sum(x)-max(x))/(length(x)-1)) 
      # criterion2 = square of worst comp's effective t-value (relative to CatMagError)
      y <- resid2 * raw_weights
      criterion2 <- max(y)
      # cat(paste0(image, "  n_comps=", length(x), 
      #            "  criterion=", criterion,  "  criterion2=", criterion2, 
      #            "  cirrus sigma=", cirrus_sigma_this_image, "\n"))
      selection <- rep(TRUE, nrow(df_estimates_this_image))
      if ((criterion1 >= 16) | (criterion2 >= 20)) {
        # Here, we will have to remove 1 or more points, up to 1/4 of them.
        c1 <- x / ((sum(x)-x)/(length(x)-1)) 
        c2 <- y
        score <- pmax(c1/16, c2/20)              # which >= 1 for any bad comp.
        max_to_remove <- floor(length(score)/4)  # but do not remove more than 1/4 of comps.
        to_remove <- (score >= 1) & 
          (row_number(score) > (length(score)-max_to_remove))  # TRUE iff one of worst comps.
        # Now, set weights to zero for worst comps, then recalculate cirrus effect and sigma.
        raw_weights[to_remove] <- 0
        normalized_weights <- raw_weights / sum(raw_weights)  # normalized to sum to one.
        v2 <- sum(normalized_weights^2)  # var(wt.mean) / wt.var(points)
        cirrus_effect_this_image <- weighted.mean(cirrus_effect_per_comp, normalized_weights)
        resid2 <- (cirrus_effect_per_comp - cirrus_effect_this_image)^2
        cirrus_sigma_this_image <- sqrt(v2*sum(normalized_weights * resid2))
        compIDs_used <- df_estimates_this_image$StarID[!to_remove]
        num_comps_used    <- length(compIDs_used)
        num_comps_removed <- nrow(df_estimates_this_image) - num_comps_used
        # Record these omitted comp stars in master comp data frame.
        removed_serials <- df_estimates_this_image[to_remove,"Serial"]
        df_estimates_comps[df_estimates_comps$Serial %in% removed_serials,"UseInEnsemble"] <- FALSE
      }
    }
    
    # Insert results into this image's row in df_cirrus_effect.
    cirrus_row_this_image <- df_cirrus_effect$Image==image
    df_cirrus_effect[cirrus_row_this_image,"CirrusEffect"]    <- cirrus_effect_this_image
    df_cirrus_effect[cirrus_row_this_image,"CirrusSigma"]     <- cirrus_sigma_this_image
    df_cirrus_effect[cirrus_row_this_image,"Criterion1"]      <- criterion1
    df_cirrus_effect[cirrus_row_this_image,"Criterion2"]      <- criterion2
    df_cirrus_effect[cirrus_row_this_image,"NumCompsUsed"]    <- num_comps_used
    df_cirrus_effect[cirrus_row_this_image,"CompIDsUsed"]     <- compIDs_used %>% paste(collapse=",")
    df_cirrus_effect[cirrus_row_this_image,"NumCompsRemoved"] <- num_comps_removed
  }
  
  # TARGET and CHECK STARS: # Get best mag estimates for ALL such observations in all filters.
  df_input_checks_targets <- df_filtered %>%
    filter(StarType %in% c("Check","Target")) %>%
    filter(MaxADU_Ur<=saturatedADU) %>%
    filter(InstMagSigma<=maxInstMagSigma) %>%
    mutate(CI=0) %>%    # because they are presumed unknown; corrected later.
    mutate(CatMagSaved=CatMag) %>%  # because we could want this after running predict().
    mutate(CatMag=0)    # estimated below by imputation from predicted mag.
  df_estimates_checks_targets <- data.frame(stringsAsFactors = FALSE)
  for (thisModelList in masterModelList) {
    df_input_this_filter <- df_input_checks_targets %>% 
      filter(Filter==thisModelList$filter) %>%
      filter(FITSfile %in% images_with_targets)
    predictOutput <- predict(thisModelList$model, newdata=df_input_this_filter, re.form= ~0) # a vector
    # Add extinction*airmass term if not included in model (as formula term "Airmass"):
    if (!thisModelList$fit_extinction) {
      predictOutput <- predictOutput + thisModelList$extinction * df_input_this_filter$Airmass
    }
    # Add transform*[Color Index] term if not included in model (as formula term "CI"):
    if (!thisModelList$fit_transform) {
      predictOutput <- predictOutput + thisModelList$transform * df_input_this_filter$CI
    }
    df_estimates_this_filter <- df_input_this_filter %>%
      select(one_of(columns_post_predict)) %>%             # defined well above
      mutate(PredictedMag = InstMag - predictOutput)
    df_estimates_checks_targets <- rbind(df_estimates_checks_targets, df_estimates_this_filter)
  } 
  df_estimates_checks_targets <- df_estimates_checks_targets %>%
    mutate(CatMag=CatMagSaved) %>%  # restoring CatMag, now that predict() has been run
    select(-CatMagSaved) %>%        # we no longer need this
    mutate(UseInEnsemble=NA)
  df_input_checks_targets <- NULL   # we no longer need this; raise an error if it's used again.
  # Here, df_estimates_checks_targets accounts for airmass(extinction), but not for color (transform),
  #    or for cirrus effect.
  
  # CIRRUS CORRRECTION: Apply per-image cirrus-effect to checks and targets (for all filters together).
  df_predictions_checks_targets <- left_join(df_estimates_checks_targets, df_cirrus_effect, 
                                             by=c("FITSfile" = "Image")) %>%
    mutate(UntransformedMag = PredictedMag - CirrusEffect) # = after cirrus correction, before transform.
  
  # COLOR CORRECTION: 
  transforms <- c(masterModelList[[CI_filters[1]]]$transform, masterModelList[[CI_filters[2]]]$transform)
  # Impute each targets' Color Index values by interpolation amongst its observations, fill them ALL in.
  df_predictions_checks_targets <- imputeTargetCIs(df_predictions_checks_targets, CI_filters, transforms)
  # Build a lookup table to get Transforms.
  df_xref_transform <- data.frame(stringsAsFactors=FALSE) 
  for (thisModelList in masterModelList) {
    df_xref_transform <- df_xref_transform %>%
      rbind(data.frame(Filter=thisModelList$filter, Transform=thisModelList$transform, 
                       stringsAsFactors=FALSE))
  }
  # Perform color transformations using the interpolated Color Index values.
  df_transformed <- df_predictions_checks_targets %>%
    left_join(df_xref_transform, by="Filter") %>%
    mutate(TransformedMag = UntransformedMag - Transform * CI) %>% # minus is the correct sign: 20160206.
    filter(!is.na(TransformedMag))  # NA often from missing V or I mag, so CI=NA, so Transformed Mag=NA.
  
  # COMPUTE TotalSigma for each target and check star in each image.
  # There will be 3 contributions to TotalSigma for each Target and Check observation:
  #   (1) model_sigma: from mixed-model regression (same for all observations in this filter).
  #   (2) cirrus_sigma: from variance in ensemble comp stars (same for all obs in this image).
  #   (3) inst_mag_sigma: from shot noise etc in observation (unique to each observation).
  df_transformed$ModelSigma  <- NA   # new column.
  df_transformed$CirrusSigma <- NA  # new column.
  df_transformed$TotalSigma  <- NA   # new column.
  for (thisModelList in masterModelList) {
    this_filter <- thisModelList$filter
    images_comps <- (df_estimates_comps %>% 
      filter(Filter==this_filter))$FITSfile %>%
      unique()
    for (image in images_comps) {
      n <- df_estimates_comps %>% filter(FITSfile==image) %>% filter(UseInEnsemble==TRUE) %>%
        nrow() %>% max(1)
      model_sigma <- summary(thisModelList$model)$sigma
      cirrus_sigma <- df_cirrus_effect[df_cirrus_effect$Image == image, "CirrusSigma"]
      df_targets_checks <- df_estimates_checks_targets %>% 
        filter(FITSfile==image) %>%
        select(Serial, InstMagSigma)
      for (serial in df_targets_checks$Serial) {
        inst_mag_sigma <- df_targets_checks[df_targets_checks$Serial == serial,"InstMagSigma"]
        total_sigma <- sqrt((model_sigma^2)/n + cirrus_sigma^2 + inst_mag_sigma^2) # add in quadrature.
        thisRow <- df_transformed$Serial==serial
        df_transformed[thisRow,"ModelSigma"]  <- model_sigma
        df_transformed[thisRow,"CirrusSigma"] <- cirrus_sigma
        df_transformed[thisRow,"TotalSigma"]  <- total_sigma
      }
    }
  }

  # Clean up (sort only) of df_estimates_comps (all data for comp stars).
  df_estimates_comps <- df_estimates_comps %>%
    arrange(ModelStarID, JD_mid)
  
  # Clean up df_transformed (all data for targets & checks).
  df_transformed <- df_transformed %>%
    select(-PredictedMag, -Criterion1, -Criterion2, -UntransformedMag, -Transform) %>%
    left_join((df_master %>% select(Serial, FOV, MaxADU_Ur, FWHM, SkyADU, SkySigma)), 
              by="Serial") %>%
    arrange(ModelStarID, JD_mid)
  
  # Make a template-only report_map.txt file if report_map.txt doesn't already exist.
  path_report_map <- make_safe_path(photometry_folder, "report_map.txt")
  if (!file.exists(path_report_map)) {
    lines <- c(
      paste0(";----- This is report_map.txt for AN folder ", AN_rel_folder),
      paste0(";----- Use this file to omit and/or combine target observations from AAVSO report."),
      paste0(";----- Example directive lines:"),
      paste0(";"),
      paste0(";#TARGET  GD Cyg ; to omit this target star altogether from AAVSO report."),
      paste0(";#JD  0.233 0.311 ; to omit this JD (fractional) range from AAVSO report."),
      paste0(";#SERIAL  34 44,129  32  1202 ; to omit these 5 Serial numbers from AAVSO report."),
      paste0(";#COMBINE   80,128 ; to combine (average) these 2 Serial numbers within AAVSO report."),
      paste0(";----- Add your directive lines:"),
      paste0(";\n") 
    )
    writeLines(lines, con=path_report_map)
  }
  
  # Save df_estimates_comps and df_transformed as .Rdata files.
  path_estimates_comps <- make_safe_path(photometry_folder, "df_estimates_comps.Rdata")
  save(df_estimates_comps, file=path_estimates_comps, precheck=FALSE)
  path_transformed <- make_safe_path(photometry_folder, "df_transformed.Rdata")
  save(df_transformed, file=path_transformed, precheck=FALSE)
  
  # Write some summary lines to R console.
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
      "   1. run eclipserPlot(starID='???') if any eclipsers,\n",
      "   2. run markup Report and group/select target observations in report_map.txt,\n",
      "   3. repeat 1 & 2 as needed\n",
      "   4. run AAVSO() to make report, and\n",
      "   5. submit report to AAVSO.")

  # Save df_estimates_comps and df_transformed as .Rdata files.
  path_estimates_comps <- make_safe_path(photometry_folder, "df_estimates_comps.Rdata")
  save(df_estimates_comps, file=path_estimates_comps, precheck=FALSE)
  path_transformed <- make_safe_path(photometry_folder, "df_transformed.Rdata")
  save(df_transformed, file=path_transformed, precheck=FALSE)
  return(df_transformed)
}

eclipserComps <- function (df_in=df_predictions, fov=NULL, starID=NULL, this_filter=NULL) {
  require(dplyr, quietly=TRUE)
  df <- df_in %>% 
    filter(FOV==fov) %>% 
    filter(StarID==starID) %>%
    filter(Filter==this_filter) %>%
    select(Serial, FITSfile, CompIDsUsed)
  set <- df$CompIDsUsed %>% strsplit(",") %>% unlist() %>% unique()
  df_comps <- data.frame(CompID=set, Included=FALSE, stringsAsFactors = FALSE)

  cat("EDIT file 'pre-predict.txt' with one of the following lines:\n")
  for (numToTest in 1:nrow(df_comps)) {
    base <- df_comps$CompID[df_comps$Included]
    test <- df_comps$CompID[!df_comps$Included]
    maxImagesQualifying <- 0
    bestSet <- "(none)"
    for (thisTest in test) {
      numImagesQualifying <- 0
      testComps <- c(base, thisTest)
      for (thisRow in 1:nrow(df)) {
        thisSet <- df$CompIDsUsed[thisRow] %>% strsplit(",") %>% unlist()
        if (all(thisTest %in% thisSet)) {
          numImagesQualifying <- numImagesQualifying + 1
        }
      }
      if(numImagesQualifying > maxImagesQualifying) {
        bestTest <- thisTest
        maxImagesQualifying <- numImagesQualifying
      }
    }
    df_comps$Included[df_comps$CompID == bestTest] <- TRUE
    cat(paste("    ", sum(df_comps$Included), " comp -> ",
              maxImagesQualifying," images qualify,",
              "-->   #COMPS ", fov, ",", this_filter, ",",
              paste0(df_comps$CompID[df_comps$Included],collapse=","), "\n"))
  }
}

markupReport <- function(AN_top_folder="J:/Astro/Images/C14", AN_rel_folder=NULL) {
  # Usage: df_markupReport <- markupReport(AN_rel_folder="AN20151216")
  #    then print df_markupReport (probably in Word).
  if (is.null(AN_rel_folder)) {stop(">>>>> You must provide a AN_rel_folder, ",
                                    "e.g., AN_rel_folder='20151216'.")}
  source("C:/Dev/Photometry/$Utility.R")
  require(dplyr,quietly=TRUE)
  AN_folder   <- make_safe_path(AN_top_folder, AN_rel_folder)
  photometry_folder <- make_safe_path(AN_folder, "Photometry")
  path_df_predictions <- make_safe_path(photometry_folder, "df_transformed.Rdata")
  load(path_df_predictions)
  
  df_checkStars <- df_predictions %>% 
    filter(StarType=="Check") %>% 
    select(FITSfile, Check=StarID, CkMag=TransformedMag, CkCat=CatMag)
  JD_fract <- (df_predictions$JD_num - (df_predictions$JD_num %>% min() %>% floor())) %>% round(digits=4)
  
  df_markupReport <- df_predictions %>%
    mutate(JD_fract=JD_fract) %>%
    filter(StarType=="Target") %>%
    select(Serial, FITSfile, Target=StarID, Filter, Exp=Exposure, Mag=TransformedMag, InstMagSigma, 
           ModelSigma, CirrusSigma, TotalSigma, 
           MaxADU=MaxADU_Ur, FWHM, JD_fract) %>%
    left_join(df_checkStars, by="FITSfile") %>%
    mutate(FITSfile=strsplit(FITSfile,'.fts') %>% unlist()) %>%
    mutate(Mag=round(Mag,3)) %>%
    mutate(FWHM=round(FWHM,2)) %>%
    mutate(CkMag=round(CkMag,3)) %>%
    mutate(Inst=round(InstMagSigma*1000)) %>%
    mutate(Model=round(ModelSigma*1000)) %>%
    mutate(Cirr=round(CirrusSigma*1000)) %>%
    mutate(Err=round(TotalSigma*1000)) %>%
    select(-InstMagSigma, -ModelSigma, -CirrusSigma, -TotalSigma) %>%
    arrange(Target, FITSfile, Filter, Exp) %>%
    select(Serial, Target, FITSfile, Filter, Exp, Mag, MaxADU, FWHM, everything())
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
  df_report <- make_df_report(photometry_folder)  # holds all the data required for AAVSO report.
  
  # REPORT HEADER lines: Construct.
  out <- "#TYPE=EXTENDED" %>%
    c("#OBSCODE=DERA") %>% #       DERA = Eric Dose's observer code @ AAVSO
    c("#SOFTWARE=custom R Scripts, github/edose") %>%
    c("#DELIM=,") %>%
    c("#DATE=JD") %>%
    c("#OBSTYPE=CCD") %>%
    c(paste0("#This report of ", nrow(df_report), " observations was generated ", 
             format(Sys.time(), tz = "GMT"), " UTC", 
             " from raw data in folder ",AN_rel_folder,".")) %>%
    c(paste0("#Software version for all computation yielding this report = ", software_version)) %>%
    c("#Eric Dose, New Mexico Mira Project, ABQ, NM") %>%
    c("#") %>%
    c("#NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,GROUP,CHART,NOTES")
  
    if (nrow(df_report) >= 1) {
  # FINAL FORMATTING of all observation report fields.
    df_formatted <- df_report %>%
      mutate(TargetName = TargetName %>% trimws() %>% toupper()) %>%
      mutate(JD         = JD         %>% as.numeric() %>% round(5) %>% format(nsmall=5)) %>%
      mutate(Mag        = Mag        %>% round(3) %>% format()) %>%
      mutate(MagErr     = TotalSigma %>% round(3) %>% format()) %>%
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
      df_formatted$TargetName,
      df_formatted$JD,
      df_formatted$Mag,
      df_formatted$MagErr,
      df_formatted$Filter,
      "YES",                   # always Transformed for full model.
      "STD",                   # we use standard comp stars, not "differential" mode.
      df_formatted$CompName,   # ="ENSEMBLE" when more than one comp star for this observation.
      df_formatted$CompMag,    # "na" if Ensemble comp, else instrument comp mag (or maybe "na" even so).
      df_formatted$CheckName,
      df_formatted$CheckMag,
      df_formatted$Airmass,
      "na",                   # GROUP, not used here.
      df_formatted$Chart,
      df_formatted$Notes,
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
    nrow(df_formatted)," observations."))
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
  path_df_estimates_comps <- make_safe_path(photometry_folder, "df_estimates_comps.Rdata")
  load(path_df_estimates_comps)

  df_report <- df_transformed %>%
    select(Serial, StarType, TargetName=StarID, JD=JD_mid, Mag=TransformedMag, 
           TotalSigma, InstMagSigma, ModelSigma, CirrusSigma, Filter) %>%
    mutate(CompName=NA_character_, CompMag=NA_real_, nComps=NA_integer_, 
           CheckName=NA_character_, CheckMag=NA_real_) %>%
    mutate(Airmass=df_transformed$Airmass, Chart=df_transformed$Chart) %>%
    mutate(Notes=paste0("obs#",Serial)) %>%
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
    badSerials <- serials[!(serials %in% df_report$Serial)]
    if (length(badSerials) >= 1) {
      cat(paste0(">>>>> Bad Serial numbers (", badSerials, ") to omit: '", thisLine, "'.\n"))
      serials <- serials[!(serials %in% badSerials)]
    }
    serialsToOmit <- c(serialsToOmit, serials)
  }
  df_report <- df_report %>% filter(!Serial %in% serialsToOmit)
  cat(paste0("Omitting ", length(serialsToOmit), " observations by Serials ",
             paste0(serialsToOmit, collapse=" "), " done.\n"))
  
  # Apply #JD directives.
  omit_lines <- directiveLines[stri_startswith_fixed(directiveLines, "#JD")]
  for (thisLine in omit_lines) {
    parms <- stri_extract_all_words(thisLine) %>% unlist() %>% strsplit(",") %>% unlist()
    if (length(parms)==3) {
      minJD <- as.numeric(parms[2]) # fractional part of JD expected
      maxJD <- as.numeric(parms[3]) #   "
      if (minJD>=0 | maxJD < 2) {
        floorJD <- floor(min(as.numeric(df_report$JD)))
        df_report <- df_report %>% 
          filter((JD <= minJD+floorJD) | (JD >= maxJD+floorJD))
        cat(paste0("Omitting JD_fract range",  minJD, " to ", maxJD, "\n"))
      } else {
        cat(paste(">>>>> minJD and/or maxJD values suspect: ", thisLine))
      }
    } else {
      cat(paste(">>>>> Can't parse line: ", thisLine))
    }
  }

  # Apply check-star names and mags (lookup from df_transformed, then insert column in-place).
  df_checks <- df_transformed %>% filter(StarType=="Check") %>% select(JD_mid, StarID, TransformedMag)
  df_joined <- df_report      %>% left_join(df_checks, by=c("JD"="JD_mid"))
  df_report <- df_report      %>% mutate(CheckName=df_joined$StarID, CheckMag=df_joined$TransformedMag)
  
  ##### Apply comp-star names and mags.
  # Make df_comps by rbind() all $obs from masterModelList
  df_comps <- df_estimates_comps %>%
    select(Serial, JD_mid, StarID, ObsMag=EstimatedMag, Filter)
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
      df_report$CompName[df_report$JD==thisJD] <- "ENSEMBLE" # & leave CompMag as NA
      # df_report$CompMag[df_report$JD==thisJD]  <- NA
      df_report$Notes[df_report$JD==thisJD]    <- paste0(df_report$Notes[df_report$JD==thisJD], 
                                                         "/", nCompsThisJD, " comps")
    } else { 
      # the single-comp-star case.
      df_report$CompName[df_report$JD==thisJD] <- df_compsThisJD$StarID[1]
      df_report$CompMag[df_report$JD==thisJD]  <- df_compsThisJD$ObsMag[1]
      df_report$Notes[df_report$JD==thisJD]    <- paste0(df_report$Notes[df_report$JD==thisJD], 
                                                         "/", nCompsThisJD, " comp")
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
      cat(paste0(">>>>> No combine possible for line: '", thisLine, "'.\n"))
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
    realCheckNames <- df_combine$CheckName[!is.na(df_combine$CheckName)] # Check Stars could be absent.
    if(length(realCheckNames) >= 2) {
      if (!allSame(realCheckNames)){
        cat(paste0(">>>>> Unlike Check Stars for line: '", thisLine, "'...Combine is skipped.\n"))
        next
      }
    }
    if (!allSame(df_combine$Chart)){
      cat(paste0(">>>>> Unlike Chart IDs for line: '", thisLine, "'...Combine is skipped.\n"))
      next
    }
    if (abs(diff(range(as.numeric(df_combine$JD)))) > 1/24){ # greater than 1 hour
      cat(paste0(">>>>> Range of JD times is too large for line: '", thisLine, "'...Combine is skipped.\n"))
      next
    }
    if (abs(diff(range(df_combine$Airmass))) > 0.400){
      cat(paste0(">>>>> Range of Airmasses (", 
                 diff(range(df_combine$Airmass)), 
                 ") is too large for line: '", thisLine, "'...Combine is skipped.\n"))
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
    # Error estimation got more sophisticated with v 1.0 (July 2016).
    require(magrittr, quietly=TRUE)
    instMagSigma <- df_combine$InstMagSigma %>% pmax(0.001) %>% raise_to_power(-2) %>% 
      sum() %>% raise_to_power(-0.5)  # combined in quadrature
    modelSigma   <- df_combine$ModelSigma[1] # not combined; same for all in this image (i.e., this filter).
    nModelSigma  <- mean(df_combine$nComp)
    cirrusSigma  <- df_combine$CirrusSigma %>% pmax(0.001) %>% raise_to_power(-2) %>% 
      sum() %>% raise_to_power(-0.5)  # combined in quadrature
    df_new$TotalSigma <- sqrt((modelSigma^2)/nModelSigma + cirrusSigma^2 + instMagSigma^2) # add in quadr.
    df_new$CheckMag <- mean(df_combine$CheckMag)
    df_new$Airmass  <- mean(df_combine$Airmass)
    # df_new$Notes    <- paste0("obs#[", (paste0(df_combine$Serial,collapse=" ")), "]/",
    #                           min(df_combine$nComps), "+ comps")    
    df_new$Notes    <- paste0(length(df_combine$Serial), " obs  >=", min(df_combine$nComps), " comps")
    df_report[df_report$Serial==serialToReplace,] <- df_new[1,]
    df_report <- df_report %>% filter(!Serial %in% serialsToDelete)
    cat(paste0("Combination of Serials ", paste0(df_combine$Serial, collapse=" "), " done.\n"))
  }

  return(df_report) # return without final formatting for AAVSO report (leave that to AAVSO() function).
}

curateEclipserComps <- function(AN_top_folder="J:/Astro/Images/C14", AN_rel_folder=NULL, 
                                df_filtered_master=NULL) {
  require(dplyr, quietly=TRUE)
  source('C:/Dev/Photometry/$Utility.R')
  
  AN_folder <- make_safe_path(AN_top_folder, AN_rel_folder)
  photometry_folder <- make_safe_path(AN_folder, "Photometry")
  
  # Get pre-predict file for this AN folder, then parse directive lines.
  pre_predict_file <- make_safe_path(photometry_folder, "pre-predict.txt")
  if (!file.exists(pre_predict_file)) {
    cat("\n>>>>> WARNING: Pre-predict file", omit_file, "does not exist.\n\n")
    return (df_filtered_master)
  } else {
    lines <- readLines(pre_predict_file, warn=FALSE)   # read last line even without EOL character(s).
    for (iLine in 1:length(lines)) {
      lines[iLine] <- lines[iLine] %>% 
        strsplit(";",fixed=TRUE) %>% unlist() %>% first() %>% trimws()  # remove comments
    }
    # Detect and collect directive text lines.
    directiveLines <- lines[stri_detect_regex(lines,'^#')]
    directiveLines <- directiveLines[!is.na(directiveLines)]
  }
  
  # Process directive lines to curate df_omitted (i.e., to keep only the comp stars that user wants).
  for (thisLine in directiveLines) {
    directive <- thisLine %>% strsplit("[ \t]",fixed=FALSE) %>% unlist() %>% 
      first() %>% trimws() %>% toupper()
    value     <- thisLine %>% substring(nchar(directive)+1) %>% trimws()  # all but the directive
    if (is.na(directive)) {directive <- ""}
    if (directive=="#COMPS") {
      parms <- value %>% strsplit(",",fixed=TRUE) %>% unlist() %>% trimws() # all parms in a vector
      if (length(parms) >= 3) {
        fov     <- parms[1]
        filter  <- parms[2]
        # compIDs is a vector of the comp star StarID values TO KEEP.
        compIDs <- parms[-1:-2] %>% strsplit(" ", fixed=TRUE) %>% unlist() %>% trimws()
        # comparison star(s) to remove for this FOV and filter (comp stars user did not specify to keep)
        to_remove <- df_filtered_master %>%
          filter(StarType=="Comp") %>%
          filter(FOV==fov) %>%
          filter(Filter==filter) %>%
          filter(!StarID %in% compIDs) %>%  # note the ! (not) operator, here
          select(Serial, ModelStarID, Filter) 
        # Remove unwanted comp star observations (i.e., to_remove$Serial)
        df_filtered_master <- df_filtered_master %>% filter(!Serial %in% to_remove$Serial)
        # TODO: Remove images with too few remaining comp star observations (maybe later).
        
        
      } else {
        cat(paste(">>>>> Can't parse line: ", thisLine))
      } # if (length(parms)>=3)
    } # if (directive==#COMPS)
  } # for thisLine
  
  # Now, (optionally) remove observations with too few 
  
  return (df_filtered_master)
}

imputeTargetCIs <- function (df_predictions, CI_filters, transforms) {
  # Impute Color Index CI (=true V Mag - true I Mag) from untransformed mags.
  # REPLACES CI for Target and Check stars (nb: CI had been 0 for Targets but Catalog CI for Checks)
  require(dplyr, quietly=TRUE)
  JD_floor <- floor(min(as.numeric(df_predictions$JD_mid)))
  df_predictions <- df_predictions %>% mutate(JD_num=as.numeric(JD_mid)-JD_floor) # avoid scaling problems.
  
  # We're only interested in images with Target stars.
  StarIDs_targets_checks <- 
    (df_predictions %>% 
       filter(StarType %in% c("Target", "Check")))$ModelStarID %>%
    unique()
  
  for (thisStarID_targets_checks in StarIDs_targets_checks) { # handle one target's ModelStarID at a time.
    df_StarID_targets_checks <- df_predictions %>% 
      filter(ModelStarID==thisStarID_targets_checks) %>%
      select(Serial, ModelStarID, Filter, JD_num, CI, UntransformedMag) %>%
      arrange(JD_num)
    df_CI_points <- extractCI_points(df_StarID_targets_checks, CI_filters, transforms) %>% 
      arrange(JD_num)

    # Interpolate CI and put it in df_StarID_targets_checks
    i <- 0
    if (nrow(df_CI_points) == 0) {
      df_StarID_targets_checks$CI <- NA  # because there are no color index values to apply.
      cat(">>>>> ModelStarID=", thisStarID_targets_checks, 
          ": no CI points returned by imputeTargetCIs()\n", sep="")
    }
    if (nrow(df_CI_points) == 1) {
      df_StarID_targets_checks$CI <- df_CI_points$CI[1] # set all CI to the same value
    }
    if (nrow(df_CI_points) %in% 2:3) { # linear interpolation (with residuals if 3 points)
      m <- lm (CI ~ JD_num, data=df_CI_points)
      df_StarID_targets_checks$CI <- predict.lm(m, data.frame(JD_num=df_StarID_targets_checks$JD_num))
      # Enforce no extrapolation.
      df_StarID_targets_checks$CI[df_StarID_targets_checks$JD_num < df_CI_points$JD_num[1]] <- 
        predict.lm(m,data.frame(JD_num=df_CI_points$JD_num[1]))
      df_StarID_targets_checks$CI[df_StarID_targets_checks$JD_num > 
                                  df_CI_points$JD_num[nrow(df_CI_points)]] <-
        predict.lm(m,data.frame(JD_num=df_CI_points$JD_num[nrow(df_CI_points)]))
    }
    if (nrow(df_CI_points) >= 4) { # make and apply smoothing spline (std package "stats").
      degrees_freedom <- round(nrow(df_CI_points)/6) %>% max(4) %>% min(nrow(df_CI_points))
      thisSpline <- smooth.spline(df_CI_points$JD_num, df_CI_points$CI, df=degrees_freedom)
      df_StarID_targets_checks$CI <- predict(thisSpline, df_StarID_targets_checks$JD_num)$y
      # Enforce no extrapolation.
      df_StarID_targets_checks$CI[df_StarID_targets_checks$JD_num < df_CI_points$JD_num[1]] <- 
        predict(thisSpline, df_CI_points$JD_num[1])$y
      df_StarID_targets_checks$CI[df_StarID_targets_checks$JD_num > 
                                  df_CI_points$JD_num[nrow(df_CI_points)]] <-
        predict(thisSpline, df_CI_points$JD_num[nrow(df_CI_points)])$y
    }

    # Do the CI value replacements.
    df_predictions$CI[match(df_StarID_targets_checks$Serial, 
                            df_predictions$Serial)] <- df_StarID_targets_checks$CI
    i <- 0
  }
  return (df_predictions) # which now includes filled-in CI values and new column JD_num.
}

extractCI_points <- function (df, CI_filters, transforms) {
  ##### df must be a data frame from df_predictions, subset holding one ModelStarID.
  
  require(dplyr, quietly=TRUE)
  maxDiffJD <- 60 / (24*60) # max time between paired obs, in (Julian) days, (60 minutes as of 20160926).
  df <- df %>% filter(Filter %in% CI_filters) %>% arrange(JD_num)
  
  df_CI_point <- data.frame(stringsAsFactors = FALSE)
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
          df_CI_point <- df_CI_point %>% rbind(data.frame(JD_num=this_JD_num, CI=this_CI, 
                                                          stringsAsFactors = FALSE))
        } 
      }
    }
} 
  return (df_CI_point)
}


