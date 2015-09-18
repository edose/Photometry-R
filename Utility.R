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
  lines <- readLines(filepathname)
  key <- "Aperture;"
  ap_line <- lines[substring(lines,1,nchar(key))==key][1]
  ap <- trimws(unlist(strsplit(ap_line, ";", fixed=TRUE)))[2]
  FluxName <- paste("FLUX", trimws(as.character(ap)), sep="")
  SNRName  <- paste("SNR",  trimws(as.character(ap)), sep="")
  t <- read.table(filepathname, header=TRUE, sep=";", skip=4, strip.white=TRUE, stringsAsFactors = FALSE)
  df <- data.frame(JD=t$DATE_OBS, CatalogStar=t$NAME, Filter=t$FILTER, 
                   ImageRA=t$IMGRA, ImageDec=t$IMGDEC, X=t$IMGX, Y=t$IMGY,
                   Bkgd=t$BKGROUND, Flux=t[,FluxName], SNR=t[,SNRName], 
                   Flag=as.integer(as.hexmode(t$FLAG)), Aperture=as.integer(ap),
                   stringsAsFactors = FALSE)
  keep_rows <- (df$SNR >= minSNR) & (df$Flag == 0)
  return(df[keep_rows,])
}

get_PI_csv_aperture <- function(){
  line <- lines[substring(lines,1,nchar(key))==key][1]
  trimws(substring(line,nchar(key)+1))
}

get_one_AAVSO_sequence <- function(filepathname=
                         "C:/Dev/Photometry/VPhot/CF Cas sparse.txt"){
  df <- read.table(filepathname, header=FALSE, sep="\t", skip=0, fill=TRUE, strip.white=TRUE, 
                   stringsAsFactors = FALSE)
  df <- df[,1:5]
  colnames(df) <- c("StarID", "degRA", "degDec", "Mags", "StarType")
  df$MagU <- NA
  df$MagB <- NA
  df$MagV <- NA
  df$MagR <- NA
  df$MagI <- NA
  mag_xref <- data.frame(passband=c("MagU","MagB","MagV","MagR","MagI"), stringsAsFactors = FALSE)
  rownames(mag_xref) <- c("1024","1","2","4","8")
  
  # Build a row in curated_df for each check and comp star.
  CH_rows <- df[df$StarType %in% c("C","H"),]
  for (irow in 1:nrow(CH_rows)) {
    mag_entries <- trimws(unlist(strsplit(CH_rows$Mags[irow],  "|", fixed=TRUE)))
    for (mag in mag_entries) {
      mag_key_value <- trimws(unlist(strsplit(mag,"_",fixed=TRUE)))
      column_name                <- mag_xref[mag_key_value[1],]
      if (!is.na(column_name)) {
        CH_rows[irow,column_name] <- as.numeric(mag_key_value[2])
      }
    }
  }
  curated_df <- CH_rows                                    # Start with rows for check and comp stars.
  curated_df <- rbind(curated_df, df[df$StarType=="T",])    # Add rows for target (unknown) stars.
  curated_df$StarType[curated_df$StarType=="C"] <- "Comp"   # ...inelegant of course
  curated_df$StarType[curated_df$StarType=="H"] <- "Check"
  curated_df$StarType[curated_df$StarType=="T"] <- "Target"
  curated_df <- curated_df[order(curated_df$StarType),]
  curated_df$Mags <- NULL
  
  
  return(curated_df)
}

