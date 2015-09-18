##### Utility.R, Input/output files for Photometry
##### Eric Dose, Bois d'Arc Observatory, Kansas, USA -- begun September 18 2015.

##### This version depends on: 
#        csv files from PixInxight (AperturePhotometryEVD.js script), one per image,
#        sequence files exported from VPHOT, one per RA-Dec field, and
#        cross-reference between catalog & AAVSO star names (file manually generated).

read_one_PI_csv <- function (filepathname=
                            "C:/Dev/Photometry/PixInsight/csv/NGC 7790-S001-R001-C002-V.csv") {
  lines <- readLines(filepathname)
  key <- "Aperture;"
  ap_line <- lines[substring(lines,1,nchar(key))==key][1]
  ap <- trimws(unlist(strsplit(ap_line, ";", fixed=TRUE)))[2]
  FluxName <- paste("FLUX", trimws(as.character(ap)), sep="")
  SNRName  <- paste("SNR",  trimws(as.character(ap)), sep="")
  t <- read.table(filepathname, header=TRUE, sep=";", skip=4, strip.white=TRUE, stringsAsFactors = FALSE)
  t2 <- data.frame(JD=t$DATE_OBS, CatalogStar=t$NAME, Filter=t$FILTER, 
                   ImageRA=t$IMGRA, ImageDec=t$IMGDEC, X=t$IMGX, Y=t$IMGY,
                   Bkgd=t$BKGROUND, Flux=t[,FluxName], SNR=t[,SNRName], Flag=t$FLAG)
  t2
}

get_PI_csv_aperture <- function(){
  line <- lines[substring(lines,1,nchar(key))==key][1]
  trimws(substring(line,nchar(key)+1))
}