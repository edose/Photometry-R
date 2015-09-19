##### Utility.R, Input/output files for Photometry
##### Eric Dose, Bois d'Arc Observatory, Kansas, USA -- begun September 18 2015.

##### This version depends on: 
#        csv files from PixInxight (AperturePhotometryEVD.js script), one per image,
#        sequence files exported from VPHOT, one per RA-Dec field, and
#        cross-reference between catalog & AAVSO star names (file manually generated).

##### Reads a folder of PixInsight photometry output (.csv) files, 
#        aggregates to one master data frame and curates.
make_PI_master_df <- function (PI_folder="C:\\") {
  ##### Argument "folder" must be a folder in which every .txt file is a transform VPHOT file.
  filenames <- trimws(list.files(PI_folder, pattern=".csv$", full.names=TRUE, 
                                 recursive=FALSE, ignore.case=TRUE))
  df <- data.frame()
  for (filename in filenames){
    df <- rbind(df, get_one_PI_csv(filename)) # get next raw data frame and append it.
  }
  return(df)
}

read_one_PI_csv <- function (filepathname=
                            "C:/Dev/Photometry/PixInsight/csv/NGC 7790-S001-R001-C002-V.csv",
                            minSNR=10) {
  t <- read.table(filepathname, header=TRUE, sep=";", skip=4, strip.white=TRUE, stringsAsFactors = FALSE)
  lines <- readLines(filepathname)
  key <- "Aperture;"
  ap_line <- lines[substring(lines,1,nchar(key))==key][1]
  ap <- trimws(unlist(strsplit(ap_line, ";", fixed=TRUE)))[2]
  FluxName <- paste("FLUX", trimws(as.character(ap)), sep="")
  SNRName  <- paste("SNR",  trimws(as.character(ap)), sep="")
  df <- data.frame(JD=t$DATE_OBS, CatalogStar=t$NAME, Filter=t$FILTER, 
                   ImageRA=t$IMGRA, ImageDec=t$IMGDEC, X=t$IMGX, Y=t$IMGY,
                   Bkgd=t$BKGROUND, Flux=t[,FluxName], SNR=t[,SNRName], 
                   Flag=as.integer(as.hexmode(t$FLAG)), Aperture=as.integer(ap),
                   stringsAsFactors = FALSE)
  keep_rows <- (df$SNR >= minSNR) & (df$Flag < 16)
  return(df[keep_rows,])
}

#get_PI_csv_aperture <- function(){
#  line <- lines[substring(lines,1,nchar(key))==key][1]
#  trimws(substring(line,nchar(key)+1))
#}

get_one_AAVSO_sequence <- function(sequence=c("CF Cas sparse"),
                                   folder="C:/Dev/Photometry/VPhot/") {
  filepathname <- make_safe_path(folder, sequence, ".txt")
  df <- read.table(filepathname, header=FALSE, sep="\t", skip=0, fill=TRUE, strip.white=TRUE, 
                   stringsAsFactors = FALSE)
  df <- df[,1:5]
  colnames(df) <- c("StarID", "degRA", "degDec", "Mags", "StarType")
  df$MagU <- NA
  df$MagB <- NA
  df$MagV <- NA
  df$MagR <- NA
  df$MagI <- NA
  df$Sequence <- trimws(sequence)
  mag_xref <- data.frame(passband=c("MagU","MagB","MagV","MagR","MagI"), stringsAsFactors = FALSE)
  rownames(mag_xref) <- c("1024","1","2","4","8") # temporary lookup data frame
  
  # Build a row for each check and comp star.
  CH_rows <- df[df$StarType %in% c("C","H"),]
  for (irow in 1:nrow(CH_rows)) {
    mag_entries <- trimws(unlist(strsplit(CH_rows$Mags[irow],  "|", fixed=TRUE)))
    for (mag in mag_entries) {
      mag_key_value <- trimws(unlist(strsplit(mag,"_",fixed=TRUE)))
      column_name   <- mag_xref[mag_key_value[1],]
      if (!is.na(column_name)) {
        CH_rows[irow,column_name] <- as.numeric(mag_key_value[2])
      }
    }
  }
  curated_df <- CH_rows                                     # Start with rows for Check and Comp stars.
  curated_df <- rbind(curated_df, df[df$StarType=="T",])    # Add rows for Target (unknown) stars.
  curated_df$StarType[curated_df$StarType=="C"] <- "Comp"   # Rename types...inelegant of course
  curated_df$StarType[curated_df$StarType=="H"] <- "Check"  #  "
  curated_df$StarType[curated_df$StarType=="T"] <- "Target" #  "
  curated_df <- curated_df[order(curated_df$StarType),]     # Sort rows by star type.
  curated_df$Mags <- NULL
  return(curated_df)
}

make_sequence_master_df <- function(sequences=c("CF Cas sparse"),
                                    folder="C:/Dev/Photometry/VPhot/") {
  df <- data.frame()
  for (sequence in sequences) {
    df <- rbind(df,get_one_AAVSO_sequence(sequence=sequence, folder=folder))
  }
  return(df)
}

# load entire cross reference between Catalog star names and AAVSO sequence (comps, etc) stars.
get_xref <- function (filepathname="C:/Astro/AAVSO/xref/xref.csv") {
  df <- read.table(filepathname, header=TRUE, sep=",", skip=0, fill=FALSE, strip.white=TRUE,
                   blank.lines.skip = TRUE, stringsAsFactors = FALSE)
  }

# Find best match(es) in PI catalog for each AAVSO sequence star.
find_xref_stars <- function (PI_csv="C:/Dev/Photometry/PixInsight/csv/NGC 7790-S001-R001-C002-V.csv",
                             sequence="CF Cas") {
  PI_df  <- read_one_PI_csv(PI_csv)
  seq_df <- make_sequence_master_df(sequence)
  df <- data.frame()
  for (irow in 1:nrow(seq_df)) {
    dRA <- 3600 * cos(PI_df$ImageDec*pi/180) * (seq_df$degRA[irow] - PI_df$ImageRA)
    dDec <- 3600 * (seq_df$degDec[irow] - PI_df$ImageDec)
    d <- sqrt(dRA^2 + dDec^2)
    minRow <- which.min(d)
    df_row <- data.frame(Sequence=sequence,
                         AAVSO=seq_df$StarID[irow],
                         CatalogID=(PI_df$CatalogStar)[minRow],
                         arcsec=d[minRow],
                         X=PI_df$X[minRow],
                         Y=PI_df$Y[minRow])
   cat(seq_df$StarID[irow], minRow, PI_df$CatalogStar[minRow], d[minRow], "\n")
   df <- rbind(df, df_row)
  }
  df
}



make_safe_path <- function (folder, filename, extension="") {
  # Pastes correctly regardless of duplicated "/". Disregard extension if it's in filename.
  gsub("/+","/",paste(trimws(folder),"/",trimws(filename),trimws(extension),sep=""),fixed=FALSE)
}