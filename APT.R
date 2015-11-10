#####  APT.R -- functions to control & extract data from APT (Cal Tech aperture photometry software).
#####      This is to replace previous use of PI (PixInsight), which has proven unreliable in its
#####         identification of source signals, and is in any case too tied to catalogs, 
#####         when all we want is fast, accurate AP of our targets.

# run_APT_all(): process all FITS files in folder through APT, build master df.
run_APT_all <- function(AN_rel_folder="20151101-test") {
  require(dplyr)
  # make list of all FITS.
  AN_folder   <- make_safe_path("J:/Astro/Images/C14",AN_folder)
  FITS_folder <- make_safe_path(AN_folder,"Calibrated")
  filenames   <- trimws(list.files(FITS_folder, pattern=".fts$", full.names=FALSE, 
                                 recursive=FALSE, ignore.case=TRUE))
  
  # make data frame with FOV names
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

# run_APT_oneFITS(): process one FITS file through APT. 
#    Requires calibrated FITS file. 
#    Also requires APT source-list file (previously constructed from VPhot sequence file).
#    Writes output file "APT-[FITS filename].txt" into FITS file's directory.
run_APT_oneFITS <- function (AN_folder=NULL, FITS_file=NULL) {
  
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


