#####  Input.R -- for general photometry, all the functions needed to get from 
#####               (1) a folder of one night's FITS files generated via Bob Denny's ACP 
#####                      and a filtered optical rig, &
#####               (2) pre-prepared field-of-view text file founded on AAVSO sequences (from VPhot)
#####             through use of APT (Cal Tech's aperture photometry software), to get a master
#####             data frame (AN) of that Astronight's photometric data. 
#####             This master data frame is to be used by Model.R to build a mixed-model regression.
#####             and finally to predict unknown 
#####
##### Typical sequence will be (starting with AN folder copied directly from obs laptop/ACP):
#####    renameObject(AN_rel_folder="20151216", oldObject="XXX", newObject="YYY") probably rarely.
#####    beforeCal(AN_rel_folder="20151216")
#####    Make calibration masters and calibrate all FITS in /Uncalibrated with MaxIm DL.
#####    finishFITS(AN_rel_folder="200151216")
#####    df_master <- make_df_master(AN_rel_folder="200151216")
##### ...then start modeling with Model.R functions.

renameObject <- function(AN_top_folder="J:/Astro/Images/C14", AN_rel_folder, 
                         oldObject, newObject) {
  ## Tests OK 20151220.
  ##### For occasional use in renaming objects (FITS header *and* FITS file name), typically as first step.
  ##### Handles files in subdirectories OK. 
  ##### Typical usage:  
  #####      renameObject(AN_rel_folder="20151216", oldObject="Landolt_SA98", newObject="Std_SA98")
  require(dplyr)
  require(stringi)
  source("C:/Dev/Photometry/$Utility.R")
  require(FITSio)
  get_header_value <- function(header, key) {  # nested function.
    value <- header[which(header==key)+1]
    if (length(value)==0) value <- NA
    trimws(value)
  }  
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
  ##### Tests OK 20151220.
  ##### Run this after using MaxIm to make Calibration masters and calibrating FITS in \Uncalibrated.
  ##### Typical usage:  finishFITS(AN_rel_folder="20151216")
  source("C:/Dev/Photometry/$Utility.R")
  require(dplyr)
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
  require(FITSio)
  get_header_value <- function(header, key) {  # nested function.
    value <- header[which(header==key)+1]
    if (length(value)==0) value <- NA
    trimws(value)
  }
  
  allegedlyCalibrated <- setdiff(list.files(UncalibratedFolder, full.names=TRUE), 
                                 list.dirs(UncalibratedFolder))
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
    
  # Here, all FITS were verified to be fully calibrated, e.g., in MaxIm.
  if (!dir.exists(CalibratedFolder)) { dir.create(CalibratedFolder) }
  CopiedOK <- file.copy(allegedlyCalibrated, CalibratedFolder, copy.date=TRUE)
  if (all(CopiedOK)) {
    file.remove(allegedlyCalibrated)
    unlink(UncalibratedFolder,recursive=TRUE)
    print(paste("Calibrated FITS:", length(allegedlyCalibrated), "moved to new Calibrated folder OK."))
  } else {
    print(paste(">>>>> Problem moving", sum(!CopiedOK), "calibrated FITS to Calibrated folder."))
  }
}

make_df_master <- function(AN_top_folder="J:/Astro/Images/C14", AN_rel_folder,
                        APT_preferences_path="C:/Dev/Photometry/APT/APT-C14.pref") {
  ##### Tested OK, except new columns Vignette4 & MaxADU (from markSaturatedObs()) need verification.
  ##### Process all FITS files in folder through APT, build & return master df.
  ##### Typical usage:  df_master <- run_APT_all(AN_rel_folder="20151216")
  require(dplyr)
  source("C:/Dev/Photometry/$Utility.R")
  AN_folder   <- make_safe_path(AN_top_folder, AN_rel_folder)
  
  # make list of all FITS.
  FITS_folder <- make_safe_path(AN_folder,"Calibrated")
  FITS_files  <- trimws(list.files(FITS_folder, pattern=".fts$", full.names=FALSE, 
                                    recursive=FALSE, ignore.case=TRUE))
  FITS_paths  <- trimws(list.files(FITS_folder, pattern=".fts$", full.names=TRUE, 
                                    recursive=FALSE, ignore.case=TRUE))
  cat(">>>>> FITS folder= ",FITS_folder,"\n")
  
  # make data frame with all FOV names for this AN, one row per FITS file.
  pattern<- "^([^-]+)-{1}[[:digit:]]{4}-{1}[[:alpha:]]" # must use dash alone as FOV/Object terminator, as we prob don't use ACP-format names.
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
    mutate(FOV_file_exists=FALSE, FOV_file_exists_string="*NA*")
  for (iFOV in 1:nrow(ask_df)) {
    ask_df$FOV_file_exists[iFOV] <- !read_FOV_file(ask_df$FOVname[iFOV]) %>% is.na() %>% first()
    ask_df$FOV_file_exists_string[iFOV] <- ask_df$FOV_file_exists[iFOV] %>% 
        ifelse("OK","MISSING")
  }
  ask_df %>%
    select(FOVname,nFITS,FOV_file=FOV_file_exists_string) %>% 
    print(right=FALSE)
  if (all(ask_df$FOV_file_exists)) {
    isYES <- "Y" == (cat("Proceed? (y/n)") %>% readline() %>% trimws() %>% toupper())
  } else {
    stop("STOPPING: at least one FOV file is missing.")
  }

  APTstdout_path  <- make_safe_path(photometry_folder, "APTstdout.txt")
  APTsourcelist_path <- make_safe_path(photometry_folder,"APTsource.txt")
  APToutput_path  <- make_safe_path(photometry_folder, "APToutput.txt")
  unlink(APTstdout_path, force=TRUE) # delete any old file before we start appending.
  allAPTstdout <- character(0)
  
  # For each FOV, (1) read FOV file, (2) make APT source file,
  #    (3) run_APT_oneFITS() on each FITS file & add to master df.
  df_master <- data.frame()  # start with a null data frame and later add rows to it.

  FOVs <- FOVdf$FOVname %>% unique() %>% sort()
  for (thisFOV in FOVs) {
    cat("FOV >", thisFOV, "<\n")
    FOV_list <- read_FOV_file(thisFOV)  # all static info defining this FOV
    FOV_FITS_paths <- FOVdf %>%         # list of all full FITS paths under this FOV
      filter(FOVname==thisFOV) %>% 
      select(FITS_path) %>% 
      unlist()

    # write APT's sourcelist file (one only) for all FITS files using this FOV.
    APT_star_data <- write_APTsourcelist_file(FOV_list$star_data, APTsourcelist_path)
    
    # for each FITS file derived from this FOV, run APT to make APT output file.
    for (thisFITS_path in FOV_FITS_paths) {
      df_FITSheader <- getFITSheaderInfo(thisFITS_path)
      APTstdout <- run_APT_oneFITS(thisFITS_path, APTsourcelist_path, 
                                   APTpreferences_path="C:/Dev/Photometry/APT/APT-C14.pref",
                                   APToutput_path)
      df_APT <- parse_APToutput(APToutput_path) %>% markSaturatedObs()
      df_APT$FITSfile <- substring(df_APT$FITSpath, nchar(FITS_folder)+2)
      df_master_thisFITS <- make_df_master_thisFITS(df_APT, APT_star_data, df_FITSheader, FOV_list$FOV_data)
      df_master <- rbind(df_master, df_master_thisFITS)         # append master rows to master data frame.
      allAPTstdout <- c(allAPTstdout, APTstdout)
      print(paste(thisFITS_path %>% strsplit("/",fixed=TRUE) %>% unlist() %>% last(), 
                  " --> df_master now has", nrow(df_master), "rows"), quote=FALSE)
    }
  }

  # Construct Vignette variables (stars' squared or 4th-power distance in pixels from image center)
  # (Do it here, not in Model.R, so that they are available for predicting check and target stars.)
  Xcenter <- 3*1024/2  # for 3K x 2K chip
  Ycenter <- 2*1024/2
  XY2_corner <- Xcenter^2 + Ycenter^2  # squared distance in pixels at image corner, for normalization.
  XY4_corner <- Xcenter^4 + Ycenter^4  # 4-power of distance in pixels at image corner, for normalization.
  df_master <- df_master %>%
    mutate(Vignette =((Xpixels-Xcenter)^2+(Ypixels-Ycenter)^2) / XY2_corner) %>%
    mutate(Vignette4=((Xpixels-Xcenter)^4+(Ypixels-Ycenter)^4) / XY4_corner)
  
  # Write out entire APT text log (all runs).
  write(allAPTstdout, APTstdout_path)
  
  # Add column for old (Ur) filenames.
  df_rename <- read.table(make_safe_path(AN_folder, "Photometry/File-renaming.txt"), 
                          header=TRUE, sep=";", stringsAsFactors=FALSE, strip.white=TRUE) %>%
    select(NewFilename, RelPath) %>%
    rename(FITSfile=NewFilename, UrFITSfile=RelPath)
  df_master <- left_join(df_master, df_rename, by="FITSfile")
  
  # Sort by JD_mid then Number, add Serial number column & move it to left end.
  # Serial numbers mostly used to omit specific data points from models.
  df_master <- df_master %>%
    arrange(JD_mid, Number) %>%
    mutate(Serial=1:nrow(df_master)) %>%
    select(Serial, UseInModel, ModelStarID, FITSfile, everything())
  df_master_path <- make_safe_path(photometry_folder, "df_master.Rdata")
  save(df_master, file=df_master_path, precheck=FALSE) # recover via: load_df_master(AN_rel_folder="201..").
  cat("make_df_master() has saved df_master to", df_master_path, "\n   now returning df_master.\n")
  return(df_master)  # one row per observation, for every FITS file.
}

load_df_master <- function(AN_top_folder="J:/Astro/Images/C14", AN_rel_folder) {
  source("C:/Dev/Photometry/$Utility.R")
  AN_folder   <- make_safe_path(AN_top_folder, AN_rel_folder)
  photometry_folder <- make_safe_path(AN_folder, "Photometry")
  df_master_path <- make_safe_path(photometry_folder, "df_master.Rdata")
  load(df_master_path, envir=globalenv(), verbose=TRUE)
  cat(AN_folder, "::df_master now available in current workspace (global envir.).\n",sep="")
}



##################################################################################################
##### The following can be called directly, but normally call instead: beforeCal().

copyToUr <- function(AN_top_folder="J:/Astro/Images/C14", AN_rel_folder){
  ##### Run this before anything, except possibly after renameObject().
  ##### Typical usage:  copyToUr(AN_rel_folder="20151216")
  ##### Tests OK 20151220.
  require(dplyr)
  require(stringi)
  source("C:/Dev/Photometry/$Utility.R")
  
  AN_folder <- make_safe_path(AN_top_folder, AN_rel_folder)
  UrFolder <- make_safe_path(AN_folder, "Ur")
  if (!dir.exists(UrFolder)) { dir.create(UrFolder) }
  
  df <- data.frame(RelPath=list.files(AN_folder, full.names=FALSE, recursive=TRUE, include.dirs=FALSE),
                   stringsAsFactors = FALSE) %>%
    filter(!stri_startswith_fixed(RelPath,"AutoFlat/")) %>%
    filter(!stri_startswith_fixed(RelPath,"Calibration/")) %>%
    filter(!stri_startswith_fixed(RelPath,"Ur/")) %>%
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
  
  require(dplyr)
  require(stringi)
  source("C:/Dev/Photometry/$Utility.R")
  require(FITSio)
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

prepareForCal <- function(AN_top_folder="J:/Astro/Images/C14", AN_rel_folder) {
  ##### Tests OK 20151220.
  #####   ALWAYS run renameACP() (or similar renaming fn) on a AN folder before running this.
  ##### Typical usage:  prepareForCal(AN_rel_folder="20151216")
  
  require(dplyr)
  source("C:/Dev/Photometry/$Utility.R")
  AN_folder                 <- make_safe_path(AN_top_folder, AN_rel_folder)
  CalibrationFolder         <- make_safe_path(AN_folder, "Calibration")
  CalibrationMastersFolder  <- make_safe_path(AN_folder, "CalibrationMasters")
  AutoFlatFolder            <- make_safe_path(AN_folder, "AutoFlat")
  UncalibratedFolder        <- make_safe_path(AN_folder, "Uncalibrated")
  PhotometryFolder          <- make_safe_path(AN_folder, "Photometry")
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
    cat("Script #1 completed OK.\n")
    cat("Now --> ensure all needed flats+darks (or Masters) are in /CalibrationMasters\n")
    cat("   (get any missing Masters from prev ANs),\n")
    cat("   then in MaxIm: 'Set Calibration' to this /CalibrationMasters\n")
    cat("   then in MaxIm: 'Replace w/Masters'\n")
    cat("   then in MaxIm: load all FITS from /Uncalibrated, 'Calibrate All', 'Close All/Save All'.\n")
    cat("   Then in R: finishFITS(), then df_master<-make_df_master().\n")
  } else {
    cat(paste(">>>>> Problem completing script #1. Check above warnings to correct.\n"))
  }
}




################################################################################################
##### Below are support-only functions, not called by user. ####################################


run_APT_oneFITS <- function (thisFITS_path=NULL, 
                             APTsourcelist_path, APTpreferences_path, APToutput_path) {
# Process one FITS file through APT. 
#    Requires calibrated FITS file. 
#    Also requires APT source-list file (previously constructed from VPhot sequence file).
#    Writes output file "APToutput.txt" into the AN directory.

  APTprogram_path <- "C:/Programs/APT/APT_v2.5.8/APT.jar"
  unlink(APToutput_path, force=TRUE)   # else APT might append to old data (bad). 

  # Construct arguments.
  APT_arguments <- c("-Duser.language=en",
                     "-Duser.region=US",
                     "-mx1024M",
                     "-jar", APTprogram_path,
                     "-i", thisFITS_path,
                     "-s", APTsourcelist_path,
                     "-p", APTpreferences_path,
                     "-o", APToutput_path)  
                     # prefs file is default (set up in APT GUI)
  
  # Run APT to generate output file.
  stdOutput <- system2("java", shQuote(APT_arguments), stdout=TRUE, wait=TRUE) # Run APT program.
  return(stdOutput)      # return text output as character vector
}

write_APTsourcelist_file <- function (star_data, APTsourcelist_path) {
  APT_star_data <- star_data %>% mutate(Number=1:nrow(star_data))
  
  APT_star_data %>%
    arrange(Number) %>%          # ensure stars are in correct order in APT source file.
    select(degRA,degDec) %>%     # APT source file needs only degRA & degDec
    format(nsmall=6) %>% 
    write.table(APTsourcelist_path, quote=FALSE, row.names=FALSE, col.names=FALSE) # no header line.
  
  return(APT_star_data) # with Number to ensure join works correctly.
}

parse_APToutput <- function (APToutput_path) {
  # parse APT text output file, return data frame for this one file (one image).
  
  lines <- readLines(APToutput_path)
  
  # Parse APT output header, set up a column range for each key.
  require(stringi)
  header <- lines[3]
  headerKeys <- c("Number", "CentroidX", "CentroidY",
                  "Magnitude", "MagUncertainty", "SkyMedian", "SkySigma",
                  "RadiusCentroid", "SkyRadiusInner", "SkyRadiusOuter", "ApertureNumRejected",
                  "RadialProfileFWHM") # key "Image" is left-justified; handle separately below.
  fieldWidths <- c(6,16,16, 16,16,16,16,  12,12,12,12,  12)
  rightmostColumns <- NULL
  for (key in headerKeys) {
    rightmostColumns <- c(rightmostColumns,
                          stri_locate_first_fixed(header, paste(" ",key," ",sep=""))[2]-1) # append
  }
  leftmostColumns <- rightmostColumns - fieldWidths + 1
  FITSpath_leftColumn <- (stri_locate_first_fixed(header, " Image "))[1]+1
  
  # For each star, parse all values and add a row to data frame.
  stopString <- "End of "
  df_APT <- data.frame()    # empty data frame to begin.
  pathnames <- vector()     # empty vector to begin.
  for (line in lines[-1:-3]) {    # i.e., skip header lines (1st 3 lines)
    if (substring(line,1,nchar(stopString))==stopString) 
      break
    FITSpath <- trimws(substring(line,FITSpath_leftColumn))
    line <- line %>% substring(1,FITSpath_leftColumn-1) %>% paste(" ")
    values <- rep(NA,length(headerKeys))
    for (iKey in 1:length(headerKeys)) {
      values[iKey] <- as.numeric(substring(line,leftmostColumns[iKey],rightmostColumns[iKey]))
    }
    df_APT   <- rbind(df_APT, values)
    pathnames <- c(pathnames,FITSpath) # build separate vector of full path names.
  }
  df_APT <- cbind(df_APT, pathnames)   # add vectors of pathnames as new columns
  colnames(df_APT) <- c("Number", "Xpixels", "Ypixels",
                        "RawMagAPT", "MagUncertainty", "SkyMedian", "SkySigma", 
                        "Radius", "SkyRadiusInner", "SkyRadiusOuter", "ApNumRej", 
                        "FWHMpixels", "FITSpath")
  df_APT$FITSpath <- as.character(df_APT$FITSpath) # without as.character() it ends up as a factor (bad).
  return(df_APT %>% arrange(Number))
  }

markSaturatedObs <- function(df_APT, saturatedADU=54000) {
  # TESTED 11/29/2015.
  # Open Ur (not Calibrated) FITS and for any observation showing saturated pixels,
  # set Saturated=TRUE for that row in df_FITS and return.
  # We want to use the original (Ur) FITS file for true checking of pixel saturation --
  #    unfortunately that means we have to look up the old filename.
  # Note: APT's first pixel is numbered "1", MaxIm's is numbered "0".
  CalFITS_path <- as.character(df_APT$FITSpath[1])

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
  require(FITSio)
  header <- readFITSheader(zz)
  D <- readFITSarray(zz, header)
  close(zz)
  image <- D$imDat
  Xsize <- D$axDat$len[1]
  Ysize <- D$axDat$len[2]
  df_APT$Saturated <- FALSE  # add column to data frame with initial values.
  df_APT$MaxADU <- NA        #  "

  # now, mark Saturated=TRUE for any row (Observation) with a saturated pixel within aperture.
  for (iObj in 1:nrow(df_APT)) {
    n <- df_APT$Number[iObj]
    Xcenter <- df_APT$Xpixels[iObj]
    Ycenter <- df_APT$Ypixels[iObj]
    testRadius  <- df_APT$Radius[iObj] + 0.5
    Xlow  <- max(1, floor(Xcenter-testRadius))
    Xhigh <- min(Xsize, ceiling(Xcenter+testRadius))
    Ylow  <- max(1, floor(Ycenter-testRadius))
    Yhigh <- min(Ysize, ceiling(Ycenter+testRadius))
    for (X in Xlow:Xhigh) {
      for (Y in Ylow:Yhigh) {
        ADU <- image[X,Y]
        if ((X-Xcenter)^2+(Y-Ycenter)^2 <= testRadius^2) {
          if (ADU > saturatedADU)         { df_APT$Saturated[iObj] <- TRUE }
          if (is.na(df_APT$MaxADU[iObj])) { df_APT$MaxADU[iObj] <- ADU }
          if (ADU > df_APT$MaxADU[iObj])  { df_APT$MaxADU[iObj] <- ADU }
        } # if
      } # for Y
    } # for X
  } # for iObj
  return(df_APT)
}

getFITSheaderInfo <- function (FITS_path) {
  source('C:/Dev/Photometry/$Utility.R')
  require(FITSio)
  require(dplyr)
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

make_df_master_thisFITS <- function (df_APT, APT_star_data, df_FITSheader, FOV_data) {
  # Combine all relevant data for all APT-detected stars for this FITS file.
  # Inputs:
  #   df_APT = parsed output from APT run on this FITS file.
  #   APT_star_data = FOV star data IN ORDER given to APT (so that data can be joined herein).
  #   df_FITSheader = data from FITS file header (Filter used, JD of exposure, etc; uniform for all rows)
  #   FOV_data = data about the FOV (uniform across all rows)

  # Join APT, FOV star, and APT data  to give master data frame.
  #   CI is color index V-I; RawMagAPT is 
  
  APT_star_data <- APT_star_data %>%
    mutate(ModelStarID=paste0(FOV_data$Sequence,"_",APT_star_data$StarID)) %>%  # as "ST Tri_137"
    mutate(UseInModel=TRUE)  # default value only
  
  df <- left_join(df_APT, APT_star_data, by="Number") %>%
    cbind(df_FITSheader) %>%
    mutate(Sequence=FOV_data$Sequence, 
           Chart=FOV_data$Chart, 
           FOV_date=FOV_data$Date,
           CI=MagV-MagI) %>%
    mutate(RawMagAPT = RawMagAPT + 2.5 * log10(Exposure)) %>%
    rename(InstMag=RawMagAPT)
  
  # Make new column "CatMag" from MagX column where X is FILTER (from FITS header).
  magColumnName <- paste("Mag",df_FITSheader$Filter[1],sep="")
  columnIndex <- match(magColumnName, colnames(df))
  df$CatMag <- df[,columnIndex]  # new column
  df <- df[,-which(colnames(df) %in% c("MagU", "MagB", "MagV", "MagR", "MagI"))] # remove columns.

  return(df)
}