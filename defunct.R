##### defunct.R: Photometry code fragments lost by the wayside. Started early October 2015.

# Reads a folder of PixInsight photometry output (.csv) files, 
#    aggregates to one master data frame and curates.
# This version depends on: 
#    csv files from PixInxight (AperturePhotometryEVD.js script), one per image,
#    sequence files exported from VPHOT, one per RA-Dec field, and
#    cross-reference between catalog & AAVSO star names (file manually generated).
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
  return (df)
}