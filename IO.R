##### IO.R, Input/output files for Photometry
##### Eric Dose, Bois d'Arc Observatory, Kansas, USA -- begun August 2015.

##### read_one_VPHOT_photometry_report(): Reads one tab-delimited file from VPHOT Photometry Report,
#####    returns one R dataframe holding all data.
get_one_VPHOT_photometry_report <- function (filepathname=
                      "C:/Dev/Photometry/NGC 7790 I 1.txt"){
  lines        <- readLines(filepathname)
  target       <- parse_VPHOT_header_line(lines,"Primary target:")
  exposure     <- parse_VPHOT_header_line(lines,"Exposure time:")
  filter       <- parse_VPHOT_header_line(lines,"Filter:")
  obs_datetime <- parse_VPHOT_header_line(lines,"Observation date/time:")
  JD           <- parse_VPHOT_header_line(lines,"JD:")
  decimal_date <- parse_VPHOT_header_line(lines,"Decimal date:")
  RA           <- parse_VPHOT_header_line(lines,"R.A.:")
  Dec          <- parse_VPHOT_header_line(lines,"Dec.:")
  airmass      <- parse_VPHOT_header_line(lines,"Airmass:")
  calibration  <- parse_VPHOT_header_line(lines,"Calibration:")
  ap_radius    <- parse_VPHOT_header_line(lines,"Apeture radius:") # yes, misspelled in AAVSO report.
  VPHOT_file   <- parse_VPHOT_header_line(lines,"File name:")

  table_start_key <- "Star\tIM\tSNR\t"
  table_start_line <- which(substring(lines,1,nchar(table_start_key))==table_start_key)[1]
  table <- read.table(filepathname,skip=table_start_line-1,header=TRUE,sep="\t")
  cat("   ", filepathname,"table has",nrow(table),"rows.\n")

  # Get index of column with catalog magnitudes by finding ".mag" at end of column's name.
  colnamelengths <- nchar(colnames(table))
  mag_column <- which(substring(colnames(table),colnamelengths-3,colnamelengths)==".mag")[1]

  # X and Y have spaces for thousands separators--remove them.
  table$X <- gsub(" ", "", table$X, fixed = TRUE)
  table$Y <- gsub(" ", "", table$Y, fixed = TRUE)
  
  exposure  <- trimws(gsub(" s","",exposure,fixed=TRUE))      # Remove " s" from end of exposure.
  ap_radius <- trimws(gsub(" pixels","",ap_radius,fixed=TRUE)) # Remove " pixels" from end of ap radius.
  
  #Construct the data frame.
  df <- data.frame(
    target=target,
    star=as.character(table$Star),
    
    InstMag=table$IM,
    CatMag=table[,mag_column],
    X=as.numeric(table$X),
    Y=as.numeric(table$Y),
    SNR=table$SNR,
    Sky=table$Sky,
    
    VPHOT_file=VPHOT_file,
    Filter=filter,
    Exposure=as.numeric(exposure),
    Airmass=as.numeric(airmass),
    JD=as.numeric(JD),
    Obs_datetime=obs_datetime,
    Dec_date=decimal_date,
    RA=RA,
    Dec=Dec,
    Cal=calibration,
    Ap_radius=as.numeric(ap_radius),
    BV_color=as.numeric(table$B.V),
    
    
    stringsAsFactors=FALSE)
}

parse_VPHOT_header_line <- function(lines, key){
  line <- lines[substring(lines,1,nchar(key))==key][1]
  trimws(substring(line,nchar(key)+1))
}