#####  Input.R -- for general photometry, all the functions needed to get from 
#####               (1) a folder of one night's FITS files generated via Bob Denny's ACP 
#####                      and a filtered optical rig, &
#####               (2) pre-prepared field-of-view text file founded on AAVSO sequences (from VPhot).
#####                      no longer calling out to APT (CalTech's software) to find and measure stars,
#####                      as of March 2016.
#####             The above delivers a master data frame of that Astronight's photometric data. 
#####             This master data frame is to be used by Model.R to build a mixed-model regression.
#####             and finally to predict unknown 
#####
##### Typical sequence will be (starting with AN folder copied directly from obs laptop/ACP):
#####    checkFOVs(AN_rel_folder="20151216")
#####    renameObject(AN_rel_folder="20151216", oldObject="XXX", newObject="YYY") probably rarely.
#####    beforeCal(AN_rel_folder="20151216")
#####    (get any missing Masters from prev ANs),
#####    In MaxIm: (1) 'Set Calibration' to this /CalibrationMasters, 'Replace w/Masters'
#####              (2) Edit/File Batch and Convert, Select All /Uncalibrated, Destination Path=/Calibrated, 
#####                     check the 'Perform Calibration' box, click 'OK'.
#####    finishFITS(AN_rel_folder="200151216")
#####    df_master <- make_df_master(AN_rel_folder="20151216")
#####    df_image  <- images(AN_rel_folder="20151216")
#####    ...then start modeling with Model.R functions.

checkFOVs <- function(AN_top_folder="J:/Astro/Images/C14", AN_rel_folder=NULL) {
  require(dplyr, quietly=TRUE)
  require(stringi, quietly=TRUE)
  source("C:/Dev/Photometry/$Utility.R")
  require(FITSio, quietly=TRUE)
  get_header_value <- function(header, key) {  # nested function.
    value <- header[which(header==key)+1]
    if (length(value)==0) value <- NA
    trimws(value)
  }  
  
  AN_folder <- make_safe_path(AN_top_folder, AN_rel_folder)
  
  # Collect file names for all relevant FITS files.
  df <- data.frame(RelPath=list.files(AN_folder, full.names=FALSE, recursive=TRUE, include.dirs=FALSE),
                   stringsAsFactors = FALSE) %>%
    filter(!stri_startswith_fixed(RelPath,"AutoFlat/")) %>%
    filter(!stri_startswith_fixed(RelPath,"Calibration/")) %>%
    filter(!stri_startswith_fixed(RelPath,"Ur/")) %>%
    filter(!stri_startswith_fixed(RelPath,"Exclude/")) %>%
    filter(!stri_startswith_fixed(RelPath,"Excluded/")) %>%
    filter(!stri_startswith_fixed(RelPath,"CalibrationMasters/")) %>%
    filter(!stri_startswith_fixed(RelPath,"Photometry/")) %>%
    mutate(Filename="", Object="")
  
  # Collect header data from FITS files.
  anyError <- FALSE
  cat("Testing: ")
  for (iRow in 1:nrow(df)) {
    relPath <- df$RelPath[iRow]
    fullPath <- make_safe_path(AN_folder, relPath)
    filename <- strsplit(relPath,"/")[[1]] %>% last()
    
    pattern<- "^(.+)-[[:digit:]S][[:digit:]]{3}-" # gets both old and new formats of FITS filenames.
    objectFromFilename <- regmatches(filename,regexec(pattern,filename)) %>% unlist() %>% nth(2)
    fileHandle <- file(description=fullPath, open="rb")
    header <- parseHdr(readFITSheader(fileHandle))
    close(fileHandle)

    objectFromFITS <- get_header_value(header, "OBJECT")
    errorThisFile <- FALSE
    #cat(paste0("testing: ",relPath,": ", objectFromFITS, " vs ", objectFromFilename, "\n"))
    cat(paste0(iRow," "))
    if (objectFromFilename != objectFromFITS) {
      cat(paste0("\n>>>>> ", fullPath,": Object mismatch, ", objectFromFilename, 
                " vs ", objectFromFITS, sep=""))
      errorThisFile <- TRUE
    }
    if (!errorThisFile) {
      df$Filename[iRow]   <- filename
      df$Object[iRow]     <- objectFromFITS
    } else {
      anyError <- TRUE
    }
    this_FOV_list <- read_FOV_file(objectFromFITS)
    FOV_file_absent <- (this_FOV_list %>% is.na() %>% first())
    if (FOV_file_absent) {
      cat(paste0("\n>>>>> No FOV file >", objectFromFITS, "< for FITS file >", relPath, "<\n"))
      anyError <- TRUE
    }
  }
  cat("\n")
  if (anyError) {
    cat(">>>>> Check the above listing for errors.\n")
  } else {
    cat("+++++ All OK: No FITS object errors, and all required FOV files are in place.\n")
  }
}

renameObject <- function(AN_top_folder="J:/Astro/Images/C14", AN_rel_folder=NULL, 
                         oldObject, newObject) {
  ## Tests OK 20151220.
  ##### For occasional use in renaming objects (FITS header *and* FITS file name), typically as first step.
  ##### Handles files in subdirectories OK. 
  ##### Typical usage:  
  #####      renameObject(AN_rel_folder="20151216", oldObject="Landolt_SA98", newObject="Std_SA98")
  require(dplyr, quietly=TRUE)
  require(stringi, quietly=TRUE)
  source("C:/Dev/Photometry/$Utility.R")
  require(FITSio, quietly=TRUE)
  get_header_value <- function(header, key) {  # nested function.
    value <- header[which(header==key)+1]
    if (length(value)==0) value <- NA
    trimws(value)
  }  
  
  if (is.null(AN_rel_folder)) {stop(">>>>> You must provide a AN_rel_folder parm, ",
                                    "e.g., AN_rel_folder='20151216'.")}
    AN_folder <- make_safe_path(AN_top_folder, AN_rel_folder)

  # Make list of all FITS files with oldObject at beginning of file name.
  df <- data.frame(RelPath=list.files(AN_folder, full.names=FALSE, recursive=TRUE, include.dirs=FALSE),
                   stringsAsFactors = FALSE) %>%
    filter(!stri_startswith_fixed(RelPath,"AutoFlat/")) %>%
    filter(!stri_startswith_fixed(RelPath,"Calibration/")) %>%
    filter(!stri_startswith_fixed(RelPath,"Ur/")) %>%
    mutate(RelDir="", OldFilename="", NewFilename="")
  if (nrow(df)<=0) {
    cat("Stopping: folder'", AN_folder, "' does not exist or is empty.\n", sep="")
    return()
  }
  
  for (iRow in 1:nrow(df)) {
    df$OldFilename[iRow]=strsplit(df$RelPath[iRow],"/")[[1]] %>% last() # doesn't work with %>%filter()
  }
  df <- df %>%
    filter(stri_startswith_fixed(OldFilename,paste(oldObject,"-",sep=""))) %>%
    mutate(OldObject=oldObject) %>%
    mutate(Stub=substring(OldFilename, nchar(oldObject)+1)) %>%
    mutate(NewFilename=paste(newObject, Stub, sep="")) %>%
    mutate(RelDir=substring(RelPath, 1, nchar(RelPath)-nchar(OldFilename))) %>%
    mutate(newRelPath=make_safe_path(RelDir,NewFilename)) %>%
    mutate(oldFullPath=make_safe_path(AN_folder,RelPath), newFullPath=make_safe_path(AN_folder, newRelPath))
  if (nrow(df)<=0) {
    cat("Stopping: no file object '", oldObject, "' in folder '", AN_folder, "'.\n", sep="")
    return()
  }
  
  # Update FITS Object in header, write new file, & if successful then delete old file.
  cat("Renaming",nrow(df),"files: ")
  for (iRow in 1:nrow(df)) {
    fullPath <- make_safe_path(AN_folder, df$RelPath[iRow])
    fits <- readFITS(fullPath)
    fits$header <- modVal("OBJECT", newObject, paste("renamed: was '",oldObject,"'", sep=""), 
                          headerName=fits$header)
    writeFITSim16i(fits$imDat, file=df$newFullPath[iRow], axDat=fits$axDat, header=fits$header) # 16-bit
    if (file.exists(df$newFullPath[iRow])) {
      if (file.size(df$newFullPath[iRow]) >= 0.9 * file.size(df$oldFullPath[iRow])) {
        file.remove(df$oldFullPath[iRow])
        cat(iRow,"")
      }
    }
  }
  cat("Done.")
}

beforeCal <- function(AN_top_folder="J:/Astro/Images/C14", AN_rel_folder) {
  ##### Calls everything (except renameObject()) that's needed before image calibration.
  ##### Typical usage: beforeCal(AN_rel_folder="20151216")
  copyToUr(AN_rel_folder=AN_rel_folder)
  renameACP(AN_rel_folder=AN_rel_folder)
  prepareForCal(AN_rel_folder=AN_rel_folder)
}

finishFITS <- function(AN_top_folder="J:/Astro/Images/C14", AN_rel_folder) {
  ##### Tests OK 20151220; mods of 20160121 test OK.
  ##### Run this after using MaxIm to make Calibration masters and using Batch Save & Convert to
  #####    calibrate all FITS from \Uncalibrated to \Calibrated.
  ##### Typical usage:  finishFITS(AN_rel_folder="20151216")
  source("C:/Dev/Photometry/$Utility.R")
  require(dplyr, quietly=TRUE)
  AN_folder                 <- make_safe_path(AN_top_folder, AN_rel_folder)
  CalibrationFolder         <- make_safe_path(AN_folder, "Calibration")
  CalibrationMastersFolder  <- make_safe_path(AN_folder, "CalibrationMasters")
  CalibratedFolder          <- make_safe_path(AN_folder, "Calibrated")
  UncalibratedFolder        <- make_safe_path(AN_folder, "Uncalibrated")
  UrFolder                  <- make_safe_path(AN_folder, "Ur")
  
  # Delete non-master calibration files from \CalibrationMasters.
  allCalibrationMasterFiles <- list.files(CalibrationMastersFolder, full.names=FALSE)
  nonMasterFiles <- allCalibrationMasterFiles[!substr(allCalibrationMasterFiles,1,7)=="Master_"]
  if (length(nonMasterFiles)>0) {
    file.remove(make_safe_path(CalibrationMastersFolder, nonMasterFiles))
  }
  
  # Verify (via FITS headers) all files in \Uncalibrated were calibrated, then rename folder to \Calibrated.
  require(FITSio, quietly=TRUE)
  get_header_value <- function(header, key) {  # nested function.
    value <- header[which(header==key)+1]
    if (length(value)==0) value <- NA
    trimws(value)
  }
  
  allegedlyCalibrated <- setdiff(list.files(CalibratedFolder, full.names=TRUE), 
                                 list.dirs(CalibratedFolder))
  nUncalibrated <- 0
  for (thisFile in allegedlyCalibrated) {
    fileHandle <- file(description=thisFile, open="rb")
    header <- parseHdr(readFITSheader(fileHandle))
    close(fileHandle)
    calstatValue <- get_header_value(header, "CALSTAT")
    if (is.na(calstatValue)) {
      nUncalibrated <- nUncalibrated + 1
    } else {
      if(calstatValue != "BDF") {
        nUncalibrated <- nUncalibrated + 1
      }
    }
  }
  if (nUncalibrated!=0)  {
    stop(paste(">>>>> STOP: of", length(allegedlyCalibrated), "FITS files,", 
               nUncalibrated, "have not been calibrated."))
  }

  preCalibration <- setdiff(list.files(UncalibratedFolder, full.names=TRUE), 
                             list.dirs(UncalibratedFolder))
  if (!(length(allegedlyCalibrated) == length(preCalibration))) {
    cat(paste(">>>>>", length(allegedlyCalibrated), "FITS in /Calibrated folder, but",
               length(preCalibration), "FITS in /Uncalibrated. (OK if some moved to /Excluded)\n"))
    answerYES <- "Y" == (cat("Proceed? (y/n)") %>% readline() %>% trimws() %>% toupper())
    if (!answerYES) {
      stop("STOPPING at user request.")
    }
  }
  
  
  # Now, change .fit file extension from MaxIm calibration back to .fts.
  require(stringi, quietly=TRUE)
  pathsCalibrated <- setdiff(list.files(CalibratedFolder, full.names=TRUE), 
                             list.dirs(CalibratedFolder))
  fts_paths <- pathsCalibrated[stri_endswith_fixed(pathsCalibrated,".fts")]
  fit_paths <- pathsCalibrated[stri_endswith_fixed(pathsCalibrated,".fit")]
  new_fts_paths <- stri_replace_last_fixed(fit_paths, ".fit", ".fts")
  RenamedOK <- file.rename(fit_paths, new_fts_paths)
  if(all(RenamedOK)) {
    cat("In /Calibrated, renamed", length(fit_paths), ".fit file extensions to .fts extensions",
      paste0("(",ifelse(length(fts_paths)==0,"no",length(fts_paths))), ".fts file extensions already existed).")
  } else {
    stop(paste(">>>>> Problem renaming some or all of .fit files to .fts, in folder /Calibrated."))
  }
  # Now it's safe to delete /Uncalibrated and all files in it.
  file.remove(preCalibration)
  unlink(UncalibratedFolder,recursive=TRUE)
}

make_df_master <- function(AN_top_folder="J:/Astro/Images/C14", AN_rel_folder,
                        CCDcenterX=1536, CCDcenterY=1024) {
  ##### Tests OK.
  ##### Process all FITS files in folder (w/o APT as of March 2016) to get (X,Y) for FOV's (RA,Dec) pairs, 
  #####    then build & return master data frame.
  ##### Typical usage:  df_master <- make_df_master(AN_rel_folder="20151216")
  require(dplyr, quietly=TRUE)
  source("C:/Dev/Photometry/$Utility.R")
  source("C:/Dev/Photometry/Aperture.R")
  
  AN_folder   <- make_safe_path(AN_top_folder, AN_rel_folder)

  # make list of all FITS (FITS in /Calibrated folder ONLY; hide files in /Exclude to remove from workflow).
  FITS_folder <- make_safe_path(AN_folder,"Calibrated")
  if(!dir.exists(FITS_folder)) {
    stop(paste0(">>>>> folder '", FITS_folder), "' does not exist.")
  }
  FITS_files  <- trimws(list.files(FITS_folder, pattern=".fts$", full.names=FALSE, 
                                    recursive=FALSE, ignore.case=TRUE))
  FITS_paths  <- trimws(list.files(FITS_folder, pattern=".fts$", full.names=TRUE, 
                                    recursive=FALSE, ignore.case=TRUE))
  cat("FITS folder = ",FITS_folder,"\n")
  
  # make data frame with all FOV names for this AN, one row per FITS file.
  pattern<- "^([^-]+)-{1}[[:digit:]]{4}-{1}[[:alpha:]]" # must use dash alone as FOV/Object terminator, 
                                                        # as we don't use ACP-format names.
  substrings <- regmatches(FITS_files,regexec(pattern,FITS_files)) %>% 
    data.frame(stringsAsFactors=FALSE) %>% 
    t()
  rownames(substrings)<- NULL
  FOVdf <- data.frame(substrings[,2],stringsAsFactors=FALSE)
  colnames(FOVdf) <- "FOVname"
  FOVdf$FITS_file <- FITS_files
  FOVdf$FITS_path <- FITS_paths
  
  # List FOVs and counts in target FITS files, ask to continue.
  photometry_folder <- make_safe_path(AN_folder,"Photometry")
  if (!dir.exists(photometry_folder)) {
    dir.create(photometry_folder)
  }
  ask_df <- FOVdf %>% 
    group_by(FOVname) %>%
    summarize(nFITS=n()) %>% 
    as.data.frame() %>%
    mutate(FOV_file_exists=FALSE, FOV_file_exists_string="*NA*", NCheckStars=0, CheckMsg="")
  for (iFOV in 1:nrow(ask_df)) {
    this_FOV_file <- read_FOV_file(ask_df$FOVname[iFOV])
    ask_df$FOV_file_exists[iFOV] <- !(this_FOV_file %>% is.na() %>% first())
    ask_df$FOV_file_exists_string[iFOV] <- ask_df$FOV_file_exists[iFOV] %>% 
        ifelse("OK","MISSING")
    if (ask_df$FOV_file_exists[iFOV]) {
      ask_df$NCheckStars[iFOV] <- sum(this_FOV_file$star_data$StarType=="Check")
      ask_df$CheckMsg[iFOV] <- ifelse(ask_df$NCheckStars[iFOV]==1, "OK", 
                                paste("WARNING: ", ask_df$NCheckStars[iFOV], 
                                       "Check Stars (must be 1 for Target FOVs)."))
    }
  }
  ask_df %>%
    select(FOVname, nFITS, FOV_file=FOV_file_exists_string, CheckStar=CheckMsg) %>% 
    print(right=FALSE)
  if (any(!ask_df$FOV_file_exists)) {
    stop("STOPPING: at least one FOV file is missing.")
  }
  answerYES <- "Y" == (cat("Proceed? (y/n)") %>% readline() %>% trimws() %>% toupper())
  if (!answerYES) {
    stop("STOPPING at user request.")
  }

  # For each FOV, (1) read FOV file, (2) for each FITS file, each star, apply aperture and measure star,
  #    (3) add each FITS file's star data to master df.
  df_master <- data.frame()  # start with a null data frame and later add rows to it.

  FOVs <- FOVdf$FOVname %>% unique() %>% sort()
  for (thisFOV in FOVs) {
    cat("FOV >", thisFOV, "<\n")
    FOV_list <- read_FOV_file(thisFOV)  # all static info defining this FOV
    FOV_FITS_paths <- FOVdf %>%         # list of all full FITS paths under this FOV
      filter(FOVname==thisFOV) %>% 
      select(FITS_path) %>% 
      unlist()

    df_star_data_numbered <- FOV_list$star_data %>% 
      mutate(Number=1:nrow(FOV_list$star_data)) # this is needed later.
    df_RADec <- df_star_data_numbered %>% 
      select(Number, StarID, degRA, degDec)
    Rdisc  <-  8
    Rinner <- 12
    Router <- 18
    df_punch <- FOV_list$punch

        # For each FITS file derived from this FOV, run APT to make APT output file.
    for (thisFITS_path in FOV_FITS_paths) {
      df_apertures <- getXYfromWCS(df_RADec, thisFITS_path)
      df_apertures <- df_apertures %>%
        mutate(RawADUMag=NA_real_, InstMagSigma=NA_real_, 
               SkyADU=NA_real_, SkySigma=NA_real_, 
               DiscRadius=Rdisc, SkyRadiusInner=Rinner, SkyRadiusOuter=Router,
               FWHM=NA_real_)
      thisImage <- getImageMatrix(thisFITS_path) # matrix of numbers. 
      Xsize <- (dim(thisImage))[1] # In FITS, X is 1st dim, Y the 2nd.
      Ysize <- (dim(thisImage))[2]
      if (nrow(df_apertures) >= 1) { # if we have any apertures at all...
      # Remove any rows that will lie outside or too near this image's boundaries.
        marginPixels <- Router + 3
        df_apertures <- df_apertures %>%
          filter((Xcentroid - marginPixels) >= 1) %>%
          filter((Ycentroid - marginPixels) >= 1) %>%
          filter((Xcentroid + marginPixels) <= Xsize) %>%
          filter((Ycentroid + marginPixels) <= Ysize)
      }
      if (nrow(df_apertures) >= 1) { # if we have apertures within image boundaries...
        df_apertures$FITSfile  <- substring(thisFITS_path, nchar(FITS_folder)+2)
        df_apertures$MaxADU_Ur <- getMaxADU_Ur(thisFITS_path, df_apertures, Rdisc)
        source("C:/Dev/Photometry/Aperture.R")
        for (i_ap in 1:nrow(df_apertures)) {
          # Make "aperture" object from this image and this aperture center (X,Y).
          Xcenter <- df_apertures$Xcentroid[i_ap]
          Ycenter <- df_apertures$Ycentroid[i_ap]
          thisAp <- makeRawAperture(thisImage, Xcenter, Ycenter, Rdisc, Rinner, Router)
          
          # Do two cycles of centroid refinement (without annulus punch) before accepting centroids.
          # (If far off from center, it may take the first cycle just to get anywhere close.)
          nCycles <- 2
          for (iCycle in 1:nCycles) {
            ev <- evalAperture(thisAp, evalSkyFunction=evalSky005)
            Xcenter <- ev$Xcentroid
            Ycenter <- ev$Ycentroid
            thisAp <- makeRawAperture(thisImage, Xcenter, Ycenter, Rdisc, Rinner, Router)
          }
          
          # Centroid position is refined once), so now prepare punch list for this star in this image.
          this_df_punch <- FOV_list$punch %>%
            filter(StarID==df_apertures$StarID[i_ap]) %>%
            mutate(dX=NA_real_, dY=NA_real_) # in pixels
          if (nrow(this_df_punch) >=1) {
            hdr <- getFITSheaderValues(thisFITS_path, c("PA", "CDELT1", "CDELT2"))
            pa_FITS <- as.numeric(hdr$PA)
            cdelt1_FITS <- as.numeric(hdr$CDELT1)
            cdelt2_FITS <- as.numeric(hdr$CDELT2)
            cat("punch:")
            for (i_punch in 1:nrow(this_df_punch)) {  # might be able to do this as vectors rather than loop.
              skyAngle <- atan2(this_df_punch$DEast[i_punch], this_df_punch$DNorth[i_punch])
              skyToPlateRotation <- (pi/180) * (360-pa_FITS)
              plateAngle <- skyAngle + skyToPlateRotation
              arcsecPerPixel <- 3600 * mean(abs(cdelt1_FITS), abs(cdelt2_FITS))
              distArcsec <- sqrt(this_df_punch$DEast[i_punch]^2 + this_df_punch$DNorth[i_punch]^2)
              distPixels <- distArcsec / arcsecPerPixel
              this_df_punch$dX[i_punch] <- distPixels * sin(-plateAngle)
              this_df_punch$dY[i_punch] <- -distPixels * cos(+plateAngle) # minus because +Y means *down*.
              cat(" ",this_df_punch$StarID[i_punch])
            }
            cat("\n")
          }
          # Apply punch list to raw aperture
          thisApPunched <- punch(thisAp, this_df_punch) 
          
          # for PRESENTATION: plot the image & punched aperture
          if (nrow(this_df_punch)>=1) { 
            end_FITS_name = substring(df_apertures$FITSfile[i_ap],
                                      nchar(df_apertures$FITSfile[i_ap])-5, 
                                      nchar(df_apertures$FITSfile[i_ap]))
            if (end_FITS_name == "-V.fts") {
              plotAperturePresentation(thisApPunched, 
                                       title = paste0(df_apertures$FITSfile[i_ap], 
                                                      ":    target ", df_apertures$StarID[i_ap]))
              iiii <- 1  # dummy stmt for bookmark during debugging.
              }
            }
          cat(paste0(df_apertures$FITSfile[i_ap], ":    target ", df_apertures$StarID[i_ap]), " --> ",
              nrow(this_df_punch), " punches.\n")
          
          # For TESTING: plot image & punched aperture
          #plotAperture(thisApPunched, title = paste0(df_apertures$FITSfile[i_ap], ":    target ",
          #  df_apertures$StarID[i_ap]))
          
          # Get values from this aperture and image matrix.
          ev <- evalAperture(thisApPunched, evalSkyFunction=evalSky005)
          
          # Replace values in this row of df_apertures.
          df_apertures$Xcentroid[i_ap]      <- ev$Xcentroid # update centroid.
          df_apertures$Ycentroid[i_ap]      <- ev$Ycentroid #   " "
          df_apertures$RawADUMag[i_ap]      <- -2.5 * log10(ev$netFlux) # NB: ev$netFlux=NA for undet star.
          df_apertures$InstMagSigma[i_ap]   <- (2.5/log(10)) * (ev$netFluxSigma / ev$netFlux)
          df_apertures$SkyADU[i_ap]         <- ev$skyADU
          df_apertures$SkySigma[i_ap]       <- ev$skySigma
          df_apertures$FWHM[i_ap]           <- ev$FWHM
        }
        df_apertures <- df_apertures %>% filter(!is.na(RawADUMag)) # remove rows for undetected stars.
      
        # Combine different data sets to make one data frame for this image.
        df_FITSheader <- getFITSheaderInfo(thisFITS_path)
        df_master_thisFITS <- make_df_master_thisFITS(df_apertures, df_star_data_numbered, df_FITSheader, 
                                                      FOV_list$FOV_data)
        # Append all this image's data to master data frame.
        df_master <- rbind(df_master, df_master_thisFITS)
        
        cat(paste(thisFITS_path %>% strsplit("/",fixed=TRUE) %>% unlist() %>% last(), 
                    " --> df_master now has", nrow(df_master), "rows\n"))
      } else {
        cat(paste(thisFITS_path %>% strsplit("/",fixed=TRUE) %>% unlist() %>% last(), 
                    " --> NO ROWS from this file ... df_master now has", nrow(df_master), "rows\n"))
      }
    }
  }

  # Construct Vignette variable (stars' squared distance in pixels from image center)
  #    and X1024 & Y1024 (well-scaled pixel distances from CCD center), for all stars incl check and target.
  df_master <- df_master %>%
    mutate(X1024=(Xcentroid-CCDcenterX)/1024, Y1024=(Ycentroid-CCDcenterY)/1024) %>%
    mutate(Vignette =(X1024^2 + Y1024^2))
  
  # Add SkyBias column.
  df_master <- addBiasTerm(df_master, biasType="sigma")
  
  # Write out entire APT text log (all runs).
  # write(allAPTstdout, APTstdout_path)
  
  # Add column for old (Ur) filenames.
  df_rename <- read.table(make_safe_path(AN_folder, "Photometry/File-renaming.txt"), 
                          header=TRUE, sep=";", stringsAsFactors=FALSE, strip.white=TRUE) %>%
    select(NewFilename, RelPath) %>%
    rename(FITSfile=NewFilename, UrFITSfile=RelPath)
  df_master <- left_join(df_master, df_rename, by="FITSfile")
  
  # Sort by JD_mid then Number, add master Serial number column & move it to left end.
  # (Serial numbers mostly used to omit specific data points from models.)
  df_master <- df_master %>%
    arrange(JD_mid, Number) %>%
    mutate(Serial=1:nrow(df_master)) %>%
    select(Serial, ModelStarID, FITSfile, everything()) %>%
    select(-Number)
  df_master_path <- make_safe_path(photometry_folder, "df_master.Rdata")
  save(df_master, file=df_master_path, precheck=FALSE) # recover via: load_df_master(AN_rel_folder="201..").
  cat("make_df_master() has saved df_master to", df_master_path, "\n   now returning df_master.\n")
  
  # Make a template-only omit.txt file if omit.txt doesn't already exist.
  PhotometryFolder <- make_safe_path(AN_folder, "Photometry")
  omitPath <- make_safe_path(PhotometryFolder, "omit.txt")
  if (!file.exists(omitPath)) {
    lines <- c(
      paste0(";----- This is omit.txt for AN folder ", AN_rel_folder),
      paste0(";----- Use this file to omit observations from input to modelOneFilter()."),
      paste0(";----- Example directive lines:\n"),
      paste0(";#OBS   Obj-0000-V, 132 ; to omit star 132 from FITS image Obj-0000-V.fts"),
      paste0(";#STAR  Obj, 132, V     ; to omit star 132 from all FITS with object Obj and filter V"),
      paste0(";#STAR  Obj, 132        ; to omit star 132 from all FITS with object Obj and ALL filters"),
      paste0(";#IMAGE Obj-0000-V      ; to omit FITS image Obj-0000-V specifically"),
      paste0(";#JD    0.72, 1         ; to omit fractional JD from 0.72 through 1"),
      paste0(";\n;----- Add your directive lines:\n;\n\n")
    )
    writeLines(lines, con=omitPath)
  }
  return(df_master)  # one row per observation, for every FITS file.
}

images <- function(AN_top_folder="J:/Astro/Images/C14", AN_rel_folder=NULL) {
  source("C:/Dev/Photometry/$Utility.R")
  AN_folder   <- make_safe_path(AN_top_folder, AN_rel_folder)
  photometry_folder <- make_safe_path(AN_folder, "Photometry")
  df_master_path <- make_safe_path(photometry_folder, "df_master.Rdata")
  load(df_master_path)
  
  maxInstMagSigma <- 0.03
  maxColorIndex <- 2.5
  saturatedADU <- 54000
  
  df <- df_master %>% 
    filter(StarType=="Comp") %>% 
    filter(InstMagSigma<=maxInstMagSigma) %>%
    filter(!is.na(CI)) %>%
    filter(CI<=maxColorIndex) %>%
    filter(MaxADU_Ur<=saturatedADU) %>%
    filter(!is.na(Airmass)) %>%
    filter(!is.na(CatMag))
  df_image <- df %>%
    group_by(FITSfile) %>% 
    summarize(nComp=n()) %>% 
    as.data.frame() %>%
    left_join(df_master %>% select(FOV, FITSfile, Filter, Exposure, Airmass)) %>% 
    unique() %>%
    select(FOV, Filter, Exposure, Airmass, FITSfile, nComp) %>%
    arrange(FOV, Filter, desc(Exposure), FITSfile)
  return (df_image)
}

load_df_master <- function(AN_top_folder="J:/Astro/Images/C14", AN_rel_folder) {
  source("C:/Dev/Photometry/$Utility.R")
  AN_folder   <- make_safe_path(AN_top_folder, AN_rel_folder)
  photometry_folder <- make_safe_path(AN_folder, "Photometry")
  df_master_path <- make_safe_path(photometry_folder, "df_master.Rdata")
  load(df_master_path, envir=globalenv())
  cat(AN_folder, "::df_master now available in current workspace (global envir.).\n",sep="")
}


##################################################################################################
##### The following can be called directly, but normally call instead: beforeCal().

copyToUr <- function(AN_top_folder="J:/Astro/Images/C14", AN_rel_folder=NULL){
  ##### Run this before anything, except possibly after renameObject().
  ##### Typical usage:  copyToUr(AN_rel_folder="20151216")
  ##### Tests OK 20151220.
  require(dplyr, quietly=TRUE)
  require(stringi, quietly=TRUE)
  source("C:/Dev/Photometry/$Utility.R")
  
  if (is.null(AN_rel_folder)) { stop("AN_rel_folder may not be null.") }
  AN_folder <- make_safe_path(AN_top_folder, AN_rel_folder)
  UrFolder <- make_safe_path(AN_folder, "Ur")
  if (!dir.exists(UrFolder)) { dir.create(UrFolder) }
  
  df <- data.frame(RelPath=list.files(AN_folder, full.names=FALSE, recursive=TRUE, include.dirs=FALSE),
                   stringsAsFactors = FALSE) %>%
    filter(!stri_startswith_fixed(RelPath,"AutoFlat/")) %>%
    filter(!stri_startswith_fixed(RelPath,"Calibration/")) %>%
    filter(!stri_startswith_fixed(RelPath,"Ur/")) %>%
    filter(!stri_startswith_fixed(RelPath,"Exclude/")) %>%
    filter(!stri_startswith_fixed(RelPath,"Excluded/")) %>%
    mutate(OldFullPath=RelPath %>% make_safe_path(AN_folder,.)) %>%
    mutate(NewFullPath=UrFolder %>% make_safe_path(RelPath)) %>%
    mutate(NewFullDir="", RelDir="")
  # Make new subdirectories if any.
  for (iRow in 1:nrow(df)) {
    filename <- strsplit(df$RelPath[iRow],"/")[[1]] %>% unlist() %>% last() # doesn't work with %>%filter()
    df$NewFullDir[iRow] <- substring(df$NewFullPath[iRow],1,nchar(df$NewFullPath[iRow])-nchar(filename)-1)
    df$RelDir[iRow] <- substring(df$RelPath[iRow],1,nchar(df$RelPath[iRow])-nchar(filename)-1)
  }
  subDirs <- df %>% filter(RelDir!="") %>% select(RelDir) %>% unique() %>% unlist()
  if (length(subDirs>=1)) {
      dir.create(make_safe_path(UrFolder,subDirs))
  }
  CopiedOK <- file.copy(df$OldFullPath, df$NewFullPath)
  cat("CopyToUr() has copied",sum(CopiedOK),"of",length(CopiedOK),"files to",UrFolder,"\n")
}

renameACP <- function(AN_top_folder="J:/Astro/Images/C14", AN_rel_folder) {
  ##### Tests OK 20151220.
  ##### Rename all FITS (including in subfolders) from ACP names to serial names.
  ##### Typical Usage:  renameACP(AN_rel_folder="20151216")
  
  require(dplyr, quietly=TRUE)
  require(stringi, quietly=TRUE)
  source("C:/Dev/Photometry/$Utility.R")
  require(FITSio, quietly=TRUE)
  get_header_value <- function(header, key) {  # nested function.
    value <- header[which(header==key)+1]
    if (length(value)==0) value <- NA
    trimws(value)
  }  
  
  AN_folder <- make_safe_path(AN_top_folder, AN_rel_folder)
  
  # Collect file names for all relevant FITS files.
  df <- data.frame(RelPath=list.files(AN_folder, full.names=FALSE, recursive=TRUE, include.dirs=FALSE),
                   stringsAsFactors = FALSE) %>%
    filter(!stri_startswith_fixed(RelPath,"AutoFlat/")) %>%
    filter(!stri_startswith_fixed(RelPath,"Calibration/")) %>%
    filter(!stri_startswith_fixed(RelPath,"Ur/")) %>%
    filter(!stri_startswith_fixed(RelPath,"Exclude/")) %>%
    filter(!stri_startswith_fixed(RelPath,"Excluded/")) %>%
    mutate(OldRelDir="", NewFilename="", Object="", Filter="", JD_start=NA, Airmass=NA)

  # Collect header data from FITS files.
  nErrors <- 0
  for (iRow in 1:nrow(df)) {
    relPath <- df$RelPath[iRow]
    fullPath <- make_safe_path(AN_folder, relPath)
    oldFilename <- strsplit(relPath,"/")[[1]] %>% last()
    
    pattern<- "^(.+)-S[[:digit:]]{3}-R"
    objectFromFilename <- regmatches(oldFilename,regexec(pattern,oldFilename)) %>% unlist() %>% nth(2)
    pattern <- "-R[[:digit:]]{3}-C[[:digit:]]{3}-([[:alpha:]]{1,6})[._]"
    filterFromFilename <- regmatches(oldFilename,regexec(pattern,oldFilename)) %>% unlist() %>% nth(2)
    
    fileHandle <- file(description=fullPath, open="rb")
    header <- parseHdr(readFITSheader(fileHandle))
    close(fileHandle)
    objectFromFITS <- get_header_value(header, "OBJECT")
    filterFromFITS <- get_header_value(header, "FILTER")
    JD_start       <- get_header_value(header, "JD")
    airmass        <- get_header_value(header, "AIRMASS")
    errorThisFile <- FALSE
    if (objectFromFilename != objectFromFITS) {
      cat(paste(fullPath,": Object mismatch, ", objectFromFilename, " vs ", objectFromFITS, sep=""))
      nErrors <- nErrors + 1
      errorThisFile <- TRUE
    }
    if (filterFromFilename!=filterFromFITS) {
      cat(paste(fullPath,": Filter mismatch, ", filterFromFilename, " vs ", filterFromFITS, sep=""))
      nErrors <- nErrors + 1
      errorThisFile <- TRUE
    }
    if (!errorThisFile) {
      df$OldRelDir[iRow]  <- substring(relPath, 1, nchar(relPath)-nchar(oldFilename))
      df$Object[iRow]     <- objectFromFITS
      df$Filter[iRow]     <- filterFromFITS
      df$JD_start[iRow]   <- JD_start
      df$Airmass[iRow]    <- airmass  # handy for selecting FITS files e.g. for VPhot.
    }
    if (nErrors>0) {
      stop(paste("There were",nErrors,"FITS/filename mismatches in Object or Filter IDs",sep=""))
    }
  }

  # Make new filenames.
  df <- df %>% arrange(Object, as.numeric(JD_start))
  iObjectFile <- 1
  for (iRow in 1:nrow(df)) {
    if (iRow >= 2) {
      if (df$Object[iRow] != df$Object[iRow-1]) {
        iObjectFile <- 1
      }
    }
    df$NewFilename[iRow] <- paste(df$Object[iRow], "-", stri_pad_left(iObjectFile,4,"0"), "-", 
                         df$Filter[iRow],".fts",sep="")
    iObjectFile <- iObjectFile + 1
  }
  
  # Rename all FITS files, moving any files in subfolders to main AN_folder.
  oldFullPath <- make_safe_path(AN_folder,df$RelPath)
  newFullPath <- make_safe_path(AN_folder,df$NewFilename) # all files end up in AN_folder
  file.rename(oldFullPath,newFullPath)

  # Remove all old subdirectories (must now be emptied by previous file.rename() call).
  dirsToRemove <- df %>% select(OldRelDir) %>% filter(nchar(OldRelDir)>0) %>% unique() %>% unlist()
  dirsToRemove <- ifelse(stri_endswith_fixed(dirsToRemove,"/"),
                         substring(dirsToRemove,1,nchar(dirsToRemove)-1),
                         dirsToRemove) # remove any trailing "/"
  if (length(dirsToRemove>=1)) {
    cat("Deleting", length(dirsToRemove), "directories:\n")
    write.table(dirsToRemove, file="", row.names=FALSE, col.names=FALSE) 
    unlink(make_safe_path(AN_folder,dirsToRemove), recursive=TRUE, force=TRUE)
  }
  
  # write data frame as text file and as R object, return data frame.
  PhotometryFolder <- make_safe_path(AN_folder, "Photometry")
  if (!dir.exists(PhotometryFolder)) { dir.create(PhotometryFolder) }
  txtFilename <- make_safe_path(PhotometryFolder, "File-renaming.txt")
  write.table(df, file=txtFilename, quote=FALSE, row.names=FALSE, sep=" ; ")  # save as text file
  R_Filename <- make_safe_path(PhotometryFolder, "File-renaming.RData")
  save(df, file=R_Filename)                                                   # save as .RData file
  cat("RenameACP() has renamed",nrow(df), "files, all of which now reside in top folder", AN_folder, "\n")
}

prepareForCal <- function(AN_top_folder="J:/Astro/Images/C14", AN_rel_folder=NULL) {
  ##### Tests OK 20151220.
  #####   ALWAYS run renameACP() (or similar renaming fn) on a AN folder before running this.
  ##### Typical usage:  prepareForCal(AN_rel_folder="20151216")
  
  if (is.null(AN_rel_folder)) { stop("AN_rel_folder may not be null.") }
  require(dplyr, quietly=TRUE)
  source("C:/Dev/Photometry/$Utility.R")
  AN_folder                 <- make_safe_path(AN_top_folder, AN_rel_folder)
  CalibrationFolder         <- make_safe_path(AN_folder, "Calibration")
  CalibrationMastersFolder  <- make_safe_path(AN_folder, "CalibrationMasters")
  AutoFlatFolder            <- make_safe_path(AN_folder, "AutoFlat")
  UncalibratedFolder        <- make_safe_path(AN_folder, "Uncalibrated")
  PhotometryFolder          <- make_safe_path(AN_folder, "Photometry")
  CalibratedFolder          <- make_safe_path(AN_folder, "Calibrated")
  scriptCompletedOK <- TRUE  # default until negated by problem.  
  
  # Make new folders as needed.
  if (!dir.exists(CalibrationFolder)) {
    if (dir.exists(AutoFlatFolder)) {
      dir.create(CalibrationFolder)   # create Calibration folder if either flats or darks exist.
    }
  }
  if (!dir.exists(CalibrationMastersFolder)) { dir.create(CalibrationMastersFolder) }
  if (!dir.exists(UncalibratedFolder))       { dir.create(UncalibratedFolder) }
  if (!dir.exists(PhotometryFolder))         { dir.create(PhotometryFolder) }
  if (!dir.exists(CalibratedFolder)) { 
    dir.create(CalibratedFolder)
  } else {
    fileList <- list.files(path=CalibratedFolder, all.files=TRUE, full.names=TRUE)
    if (length(fileList)>=1) {
      file.remove(fileList)
    }
  }
  
  # Move flats to \Calibration, remove \AutoFlat.
  if (dir.exists(AutoFlatFolder)) {
    allAutoFlatFiles <- list.files(AutoFlatFolder, full.names=TRUE)
    CopiedOK <- file.copy(allAutoFlatFiles, CalibrationFolder, copy.date=TRUE)
    if (all(CopiedOK)) {
      file.remove(allAutoFlatFiles)
      unlink(AutoFlatFolder,recursive=TRUE)
      cat(paste("Flats:", length(allAutoFlatFiles), "moved to Calibration folder OK.\n" ))
    } else {
      cat(paste(">>>>> Problem moving", sum(!CopiedOK), "Flats to Calibration folder.\n"))
      scriptCompletedOK <- FALSE
    }
  }
  
  # Copy all files (flats+darks+bias) from \Calibration to \CalibrationMasters.
  allCalibrationFiles <- list.files(CalibrationFolder, full.names=TRUE)
  CopiedOK <- file.copy(allCalibrationFiles, CalibrationMastersFolder, copy.date=TRUE)
  if (all(CopiedOK)) {
    cat(paste("Darks+bias+flats:", length(allCalibrationFiles), "moved to CalibrationMasters folder OK.\n"))
  } else {
    cat(paste(">>>>> Problem moving", sum(!CopiedOK), "Darks+bias+flats to CalibrationMasters folder.\n"))
    scriptCompletedOK <- FALSE
  }

  # Move all target files to \Uncalibrated (later renamed \Calibrated).
  allTargetFiles <- setdiff(list.files(AN_folder, full.names=TRUE),
                            list.dirs(AN_folder))  # files only, not subfolders.
  CopiedOK <- file.copy(allTargetFiles, UncalibratedFolder, copy.date=TRUE)
  if (all(CopiedOK)) {
    file.remove(allTargetFiles)
    cat("prepareForCal() has moved", length(allTargetFiles), "to Uncalibrated folder.\n")
  } else {
    cat(paste(">>>>> Problem moving", sum(!CopiedOK), "target FITS to Uncalibrated folder.\n"))
    scriptCompletedOK <- FALSE
  }
  
  if (scriptCompletedOK) {
    cat("prepareForCal() completed OK.\n")
    cat("Now --> ensure all needed flats+darks (or Masters) are in /CalibrationMasters\n")
    cat("   (get any missing Masters from prev ANs),\n")
    cat("   then in MaxIm: 'Set Calibration' to this /CalibrationMasters\n")
    cat("   then in MaxIm: 'Replace w/Masters'\n")
    cat("   then in MaxIm: Edit/File Batch and Convert, select all in /Uncalibrated,\n",
        "        Destination Path=/Calibrated, check the Perform Calibration box, click OK.\n", sep="")
    cat("   Then in R: finishFITS(), then df_master<-make_df_master().\n")
  } else {
    cat(paste(">>>>> Problem completing script #1. Check above warnings to correct.\n"))
  }
}




################################################################################################
##### Below are support-only functions, not called by user. ####################################

getMaxADU_Ur <- function(CalFITS_path, df_XY, Rdisc, saturatedADU=54000) {
  # TESTED 11/29/2015.
  # Open Ur (not Calibrated) FITS and determine max ADU within the disc radius.
  # We want to use the original (Ur) FITS file for true checking of pixel saturation --
  #    unfortunately that means we have to look up the old filename.
  # Note: FITS (& APT's) first pixel is numbered "1", MaxIm's is numbered "0".
  # Later, use Aperture.R functions rather than ~ duplicating the code here.

  # make UrFITS_path (original file name, prob ACP-format) from CalFITS_path & new file name.
  pattern <- "^(.+)/Calibrated/(.+)"
  substrings <- unlist(regmatches(CalFITS_path, regexec(pattern, CalFITS_path)))
  AN_folder <- substrings[2]
  renamedFilename <- substrings[3]
  df_rename <- read.table(make_safe_path(AN_folder, "Photometry/File-renaming.txt"),
                          header=TRUE, sep=";", stringsAsFactors=FALSE, strip.white=TRUE)
  UrFilename <- df_rename %>% filter(NewFilename==renamedFilename) %>% select(RelPath) %>% unlist()
  UrFITS_path <- paste(AN_folder, "/Ur/", UrFilename, sep="")
  
  # grab image matrix from Ur FITS file (not from Calibrated FITS file)
  zz <- file(description=UrFITS_path, open="rb")
  require(FITSio, quietly=TRUE)
  header <- readFITSheader(zz)
  D <- readFITSarray(zz, header)
  close(zz)
  image <- D$imDat
  Xsize <- D$axDat$len[1]
  Ysize <- D$axDat$len[2]
  MaxADU_Ur <- rep(NA_real_,nrow(df_XY)) # default vector to begin with.

  # now, get max ADU from Ur file, within source aperture.
  for (iObj in 1:nrow(df_XY)) {
    Xcenter     <- df_XY$Xcentroid[iObj]
    Ycenter     <- df_XY$Ycentroid[iObj]
    #testRadius  <- Rdisc + 0.5
    #Xlow  <- max(1, floor(Xcenter-testRadius))
    #Xhigh <- min(Xsize, ceiling(Xcenter+testRadius))
    #Ylow  <- max(1, floor(Ycenter-testRadius))
    #Yhigh <- min(Ysize, ceiling(Ycenter+testRadius))
    #for (X in Xlow:Xhigh) {
    #  for (Y in Ylow:Yhigh) {
    #    ADU <- image[X,Y]
    #    if ((X-Xcenter)^2+(Y-Ycenter)^2 <= testRadius^2) {
    #      if (is.na(MaxADU_Ur[iObj])) { MaxADU_Ur[iObj] <- ADU }
    #      if (ADU > MaxADU_Ur[iObj])  { MaxADU_Ur[iObj] <- ADU }
    #    } # if
    #  } # for Y
    #} # for X
    Rinner <- Rdisc + 2
    Router <- Rinner + 2
    aperture <- makeRawAperture(image=image, Xcenter=Xcenter, Ycenter=Ycenter, Rdisc=Rdisc,
                                Rinner=Rinner, Router=Router) # Rinner and Router are not used here.
    MaxADU_Ur[iObj] <- evalAperture(aperture)$maxADU
  } # for iObj
  return(MaxADU_Ur)
}

getFITSheaderInfo <- function (FITS_path) {
  source('C:/Dev/Photometry/$Utility.R')
  require(FITSio, quietly=TRUE)
  require(dplyr, quietly=TRUE)
  get_header_value <- function(header, key) {  # nested function.
    value <- header[which(header==key)+1]
    if (length(value)==0) value <- NA
    trimws(value)
  }
  fileHandle <- file(description=FITS_path, open="rb")
  header <- parseHdr(readFITSheader(fileHandle))
  close(fileHandle)

  df <- data.frame(Object=NA)
  df$Object    <- get_header_value(header, "OBJECT")
  df$JD_start  <- get_header_value(header, "JD")
  df$UTC_start <- get_header_value(header, "DATE-OBS")
  df$Exposure  <- as.numeric(get_header_value(header, "EXPOSURE"))
  df$JD_mid    <- as.character(as.numeric(df$JD_start) + ((df$Exposure/2) / (24*3600)))
  df$Filter    <- get_header_value(header, "FILTER")
  df$Airmass   <- as.numeric(get_header_value(header, "AIRMASS"))
  return (df)
}

addBiasTerm <- function (df_master, biasType="sigma") {
  # IN TESTING March 25 2016, to add SkyBias term to improve ap. measurements of very low-flux sources.
  # biasType must = "" or "sigma", for now.
  # Later, add biasTypes of:
  #    "medianOfSigma": which uses a median of skySigma rather than the individual levels.
  #    "regularized"  : which uses an average of the "sigma" and "medianOfSigma" values.
  require(dplyr)
  if ("SkyBias" %in% colnames(df_master)) {
    df_master <- df_master %>% select(-SkyBias) # remove any pre-existing column named "SkyBias"
  }
  
  if(is.null(biasType) | is.na(biasType) | biasType=="") {
    df_master$SkyBias <- 0
  } 
  else {
    biasType      <- biasType %>% trimws() %>% tolower()
    if (biasType=="sigma") {
      netDiscFlux   <- 10 ^ (-df_master$InstMag/2.5) * df_master$Exposure
      nDiscPixels   <- pi * df_master$DiscRadius^2
      skyFluxUncert <- nDiscPixels * df_master$SkySigma
      # This next line is correct--should NOT be a log10(), as (1+y/x) already approximates ln(y/x) for y<<x.
      # Using the identity log10(z) = ln(z) / ln(10): [R's log() == natural log (i.e., ln())]
      magUncert     <- abs(-2.5 * (skyFluxUncert / netDiscFlux) / log(10)) 
      df_master$SkyBias <- magUncert
    }
    else {
      df_master$SkyBias <- 0
    }
  }
  return(df_master)
}


make_df_master_thisFITS <- function (df_apertures, df_star_data_numbered, df_FITSheader, FOV_data) {
  # Combine all relevant data for all APT-detected stars for this FITS file.
  # Inputs:
  #   df_apertures = parsed output from APT run on this FITS file.
  #   df_star_data_numbered = FOV star data IN ORDER given to APT (so that data can be joined herein).
  #   df_FITSheader = data from FITS file header (Filter used, JD of exposure, etc; uniform for all rows)
  #   FOV_data = data about the FOV (uniform across all rows)

  # Join APT, FOV star, and APT data  to give master data frame.
  #   CI is color index V-I; RawADUMag is ADU-based Instrument Mag not yet corrected for exposure time.
  
  df_star_data_numbered <- df_star_data_numbered %>%
    mutate(ModelStarID=paste0(FOV_data$Sequence,"_",df_star_data_numbered$StarID))  # as "ST Tri_137"

  df <- left_join(df_apertures, df_star_data_numbered, by=c("Number", "StarID")) %>%
    cbind(df_FITSheader) %>%
    mutate(FOV=FOV_data$Sequence, 
           Chart=FOV_data$Chart, 
           FOV_date=FOV_data$Date,
           CI=MagV-MagI) %>%
    mutate(InstMag = RawADUMag + 2.5 * log10(Exposure)) %>%
    select(-Object,-RawADUMag)
  
  # Make new column "CatMag" from MagX column where X is FILTER (from FITS header).
  magColumnName <- paste("Mag",df_FITSheader$Filter[1],sep="")
  columnIndex <- match(magColumnName, colnames(df))
  df$CatMag <- df[,columnIndex]  # new column
  df <- df[,-which(colnames(df) %in% c("MagU", "MagB", "MagV", "MagR", "MagI"))] # remove columns.

  return(df)
}