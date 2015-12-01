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

#####    TODO: new function renameObject(): changes FITS filename and FITS header field to a user-specified 
#####       new object name after checking that name absent from other FITS files of that date.
#####       (Useful when an ACP plan didn't use the best Object name.)

run_APT_all <- function(AN_top_folder="J:/Astro/Images/C14", AN_rel_folder="20151101-test") {
# Process all FITS files in folder through APT, build master df.
  require(dplyr)
  # make list of all FITS.
  AN_folder   <- make_safe_path(AN_top_folder, AN_rel_folder)
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
  df_master <- data.frame()  # start with a null data frame and later add rows to it.
  FOVdf %>% 
    select(FOVname) %>% 
    group_by(FOVname) %>%
    summarize(nFITS=n()) %>% 
    as.data.frame() %>% 
    print()
  isYES <- "Y"==toupper(trimws(readline(cat("Proceed? (y/n)"))))
  print(isYES)
  if (!isYES) stop("Stopped at user request.")
  
  FOVs <- unique(FOVdf$FOVname)
  for (thisFOV in FOVs) {
    print(paste("Top of outer loop, thisFOV >", thisFOV, "<", sep=""))
    FOV_list <- read_FOV_file(thisFOV)  # all static info defining this FOV
    FOV_FITS_paths <- FOVdf %>%         # list of all full FITS paths under this FOV
      filter(FOVname==thisFOV) %>% 
      select(FITS_path) %>% 
      unlist()

    # write APT's sourcelist file (one only) for all FITS files using this FOV.
    APTsourcelist_path <- make_safe_path(AN_folder,"APTsource.txt")
    APT_star_data <- write_APTsourcelist_file(FOV_list$star_data, APTsourcelist_path)
    
    # for each FITS file derived from this FOV, run APT to make APT output file.
    for (thisFITS_path in FOV_FITS_paths) {
      df_FITSheader <- getFITSheaderInfo(thisFITS_path)
      APToutput_path <- run_APT_oneFITS(AN_folder, thisFITS_path, APTsourcelist_path) # runs APT.
      df_APT <- parse_APToutput(APToutput_path) %>%
        markSaturatedObs()
      df_master_thisFITS <- make_df_master_thisFITS(df_APT, APT_star_data, df_FITSheader, FOV_list$FOV_data)
      df_master <- rbind(df_master, df_master_thisFITS)         # append master rows to master data frame.
      print(paste("df_master: ",nrow(df_master), " rows",sep=""))
    }
  }
  save(df_master, 
       file=make_safe_path(AN_folder, "df_master.Rdata")) # may later recover via df <- load(df_master_path).
  return(df_master)  # one row per observation, for every FITS file.
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
  unlink(APToutput_path, force=TRUE)   # else APT might append to old data (bad). 
  
  APT_arguments <- c("-Duser.language=en",
                     "-Duser.region=US",
                     "-mx1024M",
                     "-jar", APTprogram_path,
                     "-i", thisFITS_path,
                     "-s", APTsourcelist_path,
                     "-o", APToutput_path)  
                     # prefs file is default (set up in APT GUI)
  
  # Run APT to generate output file.
  print("################## Call APT here.")
  errorCode <- system2("java", shQuote(APT_arguments), wait=TRUE) # Run APT program.
  if (errorCode == 0) return(APToutput_path) else return("")      # return output path.
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
  FITSpath_leftColumn <- (stri_locate_first_fixed(header, " Image "))[1]+1
  
  # For each star, parse all values and add a row to data frame.
  stopString <- "End of "
  df_APT <- data.frame()    # empty data frame to begin.
  filenames <- vector()     # empty vector to begin.
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
  df_APT$FITSpath <- as.character(df_APT$FITSpath) # else it ends up as a factor.
  return(df_APT %>% arrange(Number))
  }

markSaturatedObs <- function(df_APT, saturatedADU=55000) {
  # TESTED 11/29/2015.
  # Open Ur (not Calibrated) FITS and for any observation showing saturated pixels,
  # set Saturated=TRUE for that row in df_FITS and return.
  # We want to use the original (Ur) FITS file for true checking of pixel saturation.
  # Note: APT's first pixel is numbered "1", MaxIm's is numbered "0".
  CalFITS_path <- as.character(df_APT$FITSpath[1])

  # make UrFITS_path from CalFITS_path
  pattern <- "^(.+)/Calibrated/(.+)"
  substrings <- unlist(regmatches(CalFITS_path, regexec(pattern, CalFITS_path)))
  UrFITS_path <- paste(substrings[2],"/Ur/",substrings[3],sep="")
  
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
        if ((X-Xcenter)^2+(Y-Ycenter)^2 <= testRadius^2) {
          if (image[X,Y] > saturatedADU) {
            df_APT$Saturated[iObj] <- TRUE
          } # if
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
  df$Object   <- get_header_value(header, "OBJECT")
  df$JD_start <- as.numeric(get_header_value(header, "JD"))
  df$Exposure <- as.numeric(get_header_value(header, "EXPOSURE"))
  df$JD_mid   <- df$JD_start + (df$Exposure / (24*3600) / 2)
  df$Filter   <- get_header_value(header, "FILTER")
  df$Airmass  <- as.numeric(get_header_value(header, "AIRMASS"))
  return (df)
}

make_df_master_thisFITS <- function (df_APT, APT_star_data, df_FITSheader, FOV_data) {
  # Combine all relevant data for all APT-detected stars for this FITS file.
  # Inputs:
  #   df_APT = parsed output from APT run on this FITS file.
  #   APT_star_data = FOV star data IN ORDER given to APT (so that data can be joined herein).
  #   df_FITSheader = data from FITS file header (Filter used, JD of exposure, etc; uniform for all rows)
  #   FOV_data = data about the FOV (uniform across all rows)

  df <- left_join(df_APT, APT_star_data, by="Number") %>%
    cbind(df_FITSheader) %>%
    mutate(Sequence=FOV_data$Sequence, Chart=FOV_data$Chart, FOV_date=FOV_data$Date)
  
  # Make new column "CatMag" from MagX column where X is FILTER (from FITS header).
  magColumnName <- paste("Mag",df_FITSheader$Filter[1],sep="")
  columnIndex <- match(magColumnName, colnames(df))
  df$CatMag <- df[,columnIndex]  # new column
  df <- df[,-which(colnames(df) %in% c("MagU", "MagB", "MagV", "MagR", "MagI"))]
  return(df)
}