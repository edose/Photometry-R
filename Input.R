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

##### Always run renumberACP() and renameTarget() BEFORE calibrating the FITS files.

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
  AN_folder   <- make_safe_path("J:/Astro/Images/C14",AN_folder)
  FITS_folder <- make_safe_path(AN_folder,"Calibrated")
  filenames   <- trimws(list.files(FITS_folder, pattern=".fts$", full.names=FALSE, 
                                   recursive=FALSE, ignore.case=TRUE))
  
  # make data frame with all FOV names for this AN
  pattern<- "^(.+)-S[[:digit:]]{3}-R"
  substrings <- regmatches(filenames,regexec(pattern,filenames))
  FOVdf <- data.frame(matrix(unlist(substrings), nrow=length(substrings), byrow=T),
                      stringsAsFactors=FALSE)
  colnames(FOVdf) <- c("fragment","FOV")
  FOVdf$filename <- filenames
  FOVdf$fragment <- NULL
  
  # for each FOV, (1) read FOV file, (2) make APT source file,
  #    (3) run_APT_oneFITS() on each FITS file & add to master df.
  photometry_folder <- make_safe_path(AN_folder,"Photometry")
  if (!dir.exists(photometry_folder)) {
    dir.create(photometry_folder)
  }
  FOVs <- unique(FOVdf$target)
  for (thisFOV in FOVs) {
    # get FOV data from FOV file.
    FOV_data <- read_FOV_file(thisFOV)
    sequence_df <- FOV_data$sequence  # data frame of all objects in VPhot sequence
    chart       <- FOV_data$chart     # ID of photometry (not graphical) AAVSO chart
    
    # write APT's sourcelist file for this FOV.
    APTsourcelist_path <- make_safe_path(photometry_folder,"APTsource.txt")
    write_APTsourcelist_file(APTsourcelist_path, sequence_df)
    
    # for each FITS file derived from this FOV, run APT to make APT output file.
    FOVfiles <- FOVdf %>% 
      filter(FOV==thisFOV) %>% 
      select(filename)
    for (FOVfile in FOVfiles) {
      # TODO: set up APT and run_APT_oneFITS()
      
      
    }
    
    
  }
  
  return(df)
}


################################################################################################
##### Below are support-only functions, not called by user. ####################################



run_APT_oneFITS <- function (AN_folder=NULL, FITS_file=NULL) {
# Process one FITS file through APT. 
#    Requires calibrated FITS file. 
#    Also requires APT source-list file (previously constructed from VPhot sequence file).
#    Writes output file "APT-[FITS filename].txt" into FITS file's directory.
  APTprogram_folder <- "C:/Programs/APT/APT_v2.5.8/"
  
  # Construct arguments.
  AN_folder            <- make_safe_path("J:/Astro/Images/C14",AN_folder)
  FITS_folder          <- make_safe_path(AN_folder,"Calibrated")
  FITS_path            <- make_safe_path(FITS_folder,FITS_file)
  APTsourcelist_folder <- make_safe_path(AN_folder,"APT")
  APTsourcelist_path   <- make_safe_path(APTsourcelist_folder,)
  APToutput_path       <- make_safe_path(FITSfolder, paste("APT-",FITSfile,sep=""))
  
  APT_arguments <- c("-Duser.language=en",
                     "-Duser.region=US",
                     "-mx1024M",
                     "-jar", make_safe_path(APTprogram_folder, "APT.jar"),
                     "-i", FITS_path,
                     "-s", APTsourceListFile,
                     "-o", outputFile)
  
  # Run APT to generate output file.
  errorCode <- system2("java", shQuote(APT_arguments), wait=TRUE)
  if (errorCode == 0) return(outputFile) else return("")
}
