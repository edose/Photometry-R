#####  Input.R -- for general photometry, all the functions needed to get from 
#####               (1) a folder of one night's FITS files generated via Bob Denny's ACP 
#####                      and a filtered optical rig, &
#####               (2) pre-prepared field-of-view text file founded on AAVSO sequences (from VPhot)
#####             through use of APT (Cal Tech's aperture photometry software), to get a master
#####             data frame (AN) of that Astronight's photometric data. 
#####             This master data frame is to be used by Model.R to build a mixed-model regression.
#####             and finally to predict unknown 
#####      APT was adopted to replace PI (PixInsight), which has proven unreliable in its
#####         identification of source signals, and is in any case too tied to catalogs, 
#####         when all we want is fast, accurate AP of our comp, check, and target stars.

renumberACP <- function (folder="J:/Astro/Images/C11/2015/20150825/Renaming/") {
  # Renames FITS files from ACP naming ("T Cep-S001-R001-C001-I.fts")
  # to sequential naming ("T Cep-001-I.fts").
  # Also returns data frame of FITS files' metadata, including new file names,
  # and writes this data to a file $catalog.csv stored with the FITS files.
  df <- data.frame(ACP_name=list.files(path=folder,pattern="[.]*.f[[:alpha:]]*t"), 
                   target=NA, stringsAsFactors=FALSE)
  
  # Parse ACP filenames for all FITS files in this folder.
  pattern <- paste("(.+)-S[[:digit:]]{3}-R[[:digit:]]{3}-C[[:digit:]]{3}-",
                   "([[:alpha:]][[:alnum:]]*)[_]*(.*)f[[:alpha:]]+", sep="")
  filesDetected <- nrow(df)
  for (iRow in 1:filesDetected) {
    filename <- df$ACP_name[iRow]
    substrings <- regmatches(filename, regexec(pattern, filename))[[1]]
    df$target[iRow] <- substrings[2]
  }
  df <- df[!is.na(df$target),]  # keep only files with ACP-style names.
  filesToRename <- nrow(df)
  cat(filesDetected, "FITS files detected, of which", filesToRename, "shall be renamed.\n")
  if (filesToRename <= 0) {
    stop(cat("STOPPING: no FITS files to rename in folder", folder, 
             "(have you already renamed them?)"))
  }
  
  # Extract multiple metadata from FITS files' headers.  
  source('C:/Dev/Photometry/$Utility.R')
  require(FITSio)
  require(dplyr)
  get_header_value <- function(header, key) {  # nested function.
    value <- header[which(header==key)+1]
    if (length(value)==0) value <- NA
    value
  }
  df <- cbind(df, Object=NA, JD_start=NA, JD_mid=NA, Filter=NA, Airmass=NA, FWHM=NA, Exp_secs=NA)
  for (iRow in 1:nrow(df)) {  
    fullfilename <- make_safe_path(folder, df$ACP_name[iRow])
    fileHandle <- file(description=fullfilename, open="rb")
    header <- parseHdr(readFITSheader(fileHandle))
    close(fileHandle)
    df$Object[iRow]   <- get_header_value(header, "OBJECT")
    df$JD_start[iRow] <- get_header_value(header, "JD")
    df$Filter[iRow]   <- get_header_value(header, "FILTER")
    df$Airmass[iRow]  <- get_header_value(header, "AIRMASS")
    df$FWHM[iRow]     <- get_header_value(header, "FWHM")
    df$Exp_secs[iRow] <- get_header_value(header, "EXPTIME")
  }
  JD_half_duration <- (as.numeric(df$Exp_secs)/2) / (24*3600) # convert seconds to JD days.
  df$JD_mid <- as.character(as.numeric(df$JD_start) + JD_half_duration)
  
  # Report any mismatches between ACP name and FITS-header Object.
  mismatches <- df %>%
    select(ACP_name,target,Object) %>%
    filter(target!=Object)
  cat ("Object-target mismatches (s/b zero) =",nrow(mismatches), "\n")
  if (nrow(mismatches) > 0) {
    write.table(mismatches,file="",quote=FALSE,row.names=FALSE)
  }
  
  # Construct time-based index for files *within* each target, then construct new names.
  df2 <- df %>% 
    group_by(Object) %>% 
    arrange(JD_start) %>% 
    mutate(one=1) %>% 
    mutate(index=trimws(as.character(cumsum(one)))) %>%
    mutate(index=substring(paste("000",index,sep=""),nchar(index)+1,nchar(index)+3)) %>%
    select(-target,-one)
  df2$newName <- paste(df2$Object,"-",df2$index,"-",df2$Filter,".fts",sep="") # (dplyr mutate didn't work)
  
  # Rename all FITS files (after confirmation).
  userConfirmedRename <- "Y"==toupper(trimws(readline(
    cat("Do you really want to rename all",nrow(df2), "FITS files? (y/n)"))))
  if (userConfirmedRename) {
    df_rename <- df2 %>%
      select(ACP_name,newName) %>%
      mutate(oldName=make_safe_path(folder,ACP_name)) %>%
      mutate(newName=make_safe_path(folder,newName))
    outcome <- file.rename(df_rename$oldName, df_rename$newName)
    cat(sum(outcome),"files have been renamed...")
    write.csv(df2,file=make_safe_path(folder,"$catalog.csv"),row.names=FALSE)
    cat("and file $catalog.csv has been written.\n")
  }
  return(df2)
}

#####    TODO: new function renameObject(): changes FITS filename and FITS header field to a user-specified 
#####       new object name after checking that name absent from other FITS files of that date.

run_APT_all <- function(AN_rel_folder="20151101-test") {
# Process all FITS files in folder through APT, build master df.
  require(dplyr)
  # make list of all FITS.
  AN_folder   <- make_safe_path("J:/Astro/Images/C14",AN_rel_folder)
  FITS_folder <- make_safe_path(AN_folder,"Calibrated")
  FITS_files  <- trimws(list.files(FITS_folder, pattern=".fts$", full.names=FALSE, 
                                    recursive=FALSE, ignore.case=TRUE))
  FITS_paths  <- trimws(list.files(FITS_folder, pattern=".fts$", full.names=TRUE, 
                                    recursive=FALSE, ignore.case=TRUE))
  
  # make data frame with all FOV names for this AN, one row per FITS file.
  pattern<- "^(.+)-S[[:digit:]]{3}-R"
  substrings <- regmatches(FITS_files,regexec(pattern,FITS_files))
  FOVdf <- data.frame(matrix(unlist(substrings), nrow=length(substrings), byrow=T),
                      stringsAsFactors=FALSE) %>%
    select(-1)
  colnames(FOVdf) <- "FOVname"
  FOVdf$FITS_file <- FITS_files
  FOVdf$FITS_path <- FITS_paths
  
  # for each FOV, (1) read FOV file, (2) make APT source file,
  #    (3) run_APT_oneFITS() on each FITS file & add to master df.
  photometry_folder <- make_safe_path(AN_folder,"Photometry")
  if (!dir.exists(photometry_folder)) {
    dir.create(photometry_folder)
  }
  master_df <- data.frame()  # start with a null data frame and later add to it.
  
  FOVs <- unique(FOVdf$FOVname)
  for (thisFOV in FOVs) {
    FOV_list    <- read_FOV_file(thisFOV)

    # write APT's sourcelist file for this FOV.
    APTsourcelist_path <- make_safe_path(AN_folder,"APTsource.txt")
    write_APTsourcelist_file(FOV_list$star_data, APTsourcelist_path)
    
    # for each FITS file derived from this FOV, run APT to make APT output file.
    FOV_FITS_paths <- FOVdf %>% 
      filter(FOVname==thisFOV) %>% 
      select(FITS_path)
    for (thisFITS_path in FOV_FITS_paths) {
      APToutput_path <- run_APT_oneFITS(AN_folder, thisFITS_path, APTsourcelist_path)
      df_FITS <- parse_APToutput(APToutput_path)
      df_FITS <- removeSaturatedObs (df_FITS)
      # TODO : remove rows for obs saturated near obs (check *Ur* FITS, not Calibrated).
      # TODO : add columns for FOV, FITS etc info to df_FITS.
      master_df <- rbind (master_df, df_FITS)  
    }
  }
    
  return(master_df)  # one row per observation.
}


################################################################################################
##### Below are support-only functions, not called by user. ####################################



run_APT_oneFITS <- function (AN_folder=NULL, thisFITS_path=NULL, APT_sourcelist_path) {
# Process one FITS file through APT. 
#    Requires calibrated FITS file. 
#    Also requires APT source-list file (previously constructed from VPhot sequence file).
#    Writes output file "APToutput.txt" into the AN directory.

  # Construct arguments.
  APTprogram_path <- "C:/Programs/APT/APT_v2.5.8/APT.jar"
  APToutput_path  <- make_safe_path(AN_folder, "APToutput.txt")
  
  APT_arguments <- c("-Duser.language=en",
                     "-Duser.region=US",
                     "-mx1024M",
                     "-jar", APTprogram_path,
                     "-i", thisFITS_path,
                     "-s", APTsourcelist_path,
                     "-o", APToutput_path)  
                     # prefs file is default (set up in APT GUI)
  
  # Run APT to generate output file.
  errorCode <- system2("java", shQuote(APT_arguments), wait=TRUE) # Run APT program.
  if (errorCode == 0) return(APToutput_path) else return("")      # return output path.
}

write_APTsourcelist_file <- function (star_data, APTsourcelist_path) {
  FOV_list$star_data %>% 
    select(degRA,degDec) %>% 
    format(nsmall=6) %>% 
    write.table(APTsourcelist_path, quote=FALSE, row.names=FALSE)
}

parse_APToutput <- function (APToutput_path) {
  # parse APT text output file, return data frame for this one file (one image).
  lines <- readLines(APToutput_path)
  
  # Parse header, assign column range to each key.
  header <- lines[3]
  headerKeys <- c("Number", "CentroidX", "CentroidY", "SourceIntensity", "SourceUncertainty",
                  "Magnitude", "MagUncertainty", "SkyMedian", "SkySigma",
                  "RadiusCentroid", "SkyRadiusInner", "SkyRadiusOuter", "ApertureNumRejected",
                  "RadialProfileFWHM") # key "Image" is left-justified; handle separately.
  fieldWidths <- c(6,16,16,16,16,  16,16,16,16,  12,12,12,12,  12)
  rightmostColumns <- NULL
  for (key in headerKeys) {
    rightmostColumns <- c(rightmostColumns,
                          stri_locate_first_fixed(header, paste(" ",key," ",sep=""))[2]-1) # append
  }
  leftmostColumns <- rightmostColumns - fieldWidths + 1
  FITSpath_leftColumn <- stri_locate_first_fixed(header, paste(" Image ",sep=""))[1]+1
  stopString <- "End of "
  df_APT <- as.data.frame(NULL)   # empty data frame to begin.
  filenames <- as.vector(NULL)     # empty vector to begin.
  # For each star, parse all values and add a row to data frame.
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
    filenames <- c(filenames,FITSpath)
  }
  df_APT <- cbind(df_APT, filenames)
  colnames(df_APT) <- c("Number","Xpixels","Ypixels","Intensity","Uncertainty","Mag","MagUncertainty",
                        "SkyMedian","SkySigma","Radius","SkyRadiusInner","SkyRadiusOuter",
                        "ApNumRej","FWHMpixels","FITSpath")
  return(df_APT %>% arrange(Number))
  }

removeSaturatedObs <- function(df_FITS) {
  # Open Ur (not Calibrated) FITS and for any observation showing saturated pixels,
  # remove that row from df_FITS and return.
  CalFITS_path <- df_FITS$FITS_path
  
  
  
  return(df_FITS)
}

