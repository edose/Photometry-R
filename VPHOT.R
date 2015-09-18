##### VPhot.R, Input & data frames from AAVSO's VPhot photometry file format
##### Eric Dose, Bois d'Arc Observatory, Kansas, USA -- begun August 2015.

##### utility: reads a folder of VPhot photometry files, aggregates to one master data frame.
make_VPhot_master_df <- function (VPhotFolder="C:\\") {
  ##### Argument "folder" must be a folder in which every .txt file is a transform VPHOT file.
  filenames <- trimws(list.files(VPhotFolder, pattern=".txt$", full.names=TRUE, 
                                 recursive=FALSE, ignore.case=TRUE))
  df <- data.frame()
  for (filename in filenames){
    df <- rbind(df, get_one_VPhot_photometry_report(filename)) # get next raw data frame and append it.
  }
  return(df)
}

##### utility: returns number of images (VPHOT files) in VPHOT master data frame, 
#####    or -1 if df appears invalid.
count_images_in_VPhot_master_df<- function(VPhot_master_df){
  n_filename <- length(unique(VPhot_master_df$VPhot_file))
  n_subset <- nrow(unique(data.frame(VPhot_master_df$VPhot_file, VPhot_master_df$JD, 
                                     VPhot_master_df$target)))
  if (n_filename != n_subset) { return(-1) } 
  else                        { return (n_subset)  }
}

##### utility: reads one tab-delimited file from VPhot Photometry Report,
#####    returns one R dataframe holding all data.
get_one_VPhot_photometry_report <- function (filepathname=
                                               "C:/Dev/Photometry/NGC 7790 I 1.txt"){
  # Get header lines.
  lines        <- readLines(filepathname)
  target       <- parse_VPhot_header_line(lines,"Primary target:")
  exposure     <- parse_VPhot_header_line(lines,"Exposure time:")
  filter       <- parse_VPhot_header_line(lines,"Filter:")
  obs_datetime <- parse_VPhot_header_line(lines,"Observation date/time:")
  JD           <- parse_VPhot_header_line(lines,"JD:")
  decimal_date <- parse_VPhot_header_line(lines,"Decimal date:")
  RA           <- parse_VPhot_header_line(lines,"R.A.:")
  Dec          <- parse_VPhot_header_line(lines,"Dec.:")
  airmass      <- parse_VPhot_header_line(lines,"Airmass:")
  calibration  <- parse_VPhot_header_line(lines,"Calibration:")
  ap_radius    <- parse_VPhot_header_line(lines,"Apeture radius:") # yes, misspelled in AAVSO report.
  VPhot_file   <- parse_VPhot_header_line(lines,"File name:")
  
  #Get Target (unknown) lines, if any.
  table_start_key <- "Star\tIM\tSNR\t"
  table_start_line <- which(substring(lines,1,nchar(table_start_key))==table_start_key)[1]
  table <- read.table(filepathname,skip=table_start_line-1,header=TRUE,sep="\t")
  cat("   ", filepathname,"table has",nrow(table),"rows.\n")
  
  # Get index of column with catalog magnitudes by finding ".mag" at end of column's name.
  colnamelengths <- nchar(colnames(table))
  mag_column <- which(substring(colnames(table),colnamelengths-3,colnamelengths)==".mag")[1]
  
  # Some fields as read have units at strings' end--remove them.
  exposure  <- trimws(gsub(" s",     "",exposure, fixed=TRUE)) # Remove " s" from end of exposure.
  ap_radius <- trimws(gsub(" pixels","",ap_radius,fixed=TRUE)) # Remove " pixels" from end of ap radius.

  # Some fields have spaces for thousands separators (inside the strings)--remove them.
  exposure        <- gsub(" ", "", exposure,  fixed = TRUE)
  JD              <- gsub(" ", "", JD,        fixed = TRUE)
  airmass         <- gsub(" ", "", airmass,   fixed = TRUE)
  ap_radius       <- gsub(" ", "", ap_radius, fixed = TRUE)  
  table$X         <- gsub(" ", "", table$X,   fixed = TRUE)
  table$Y         <- gsub(" ", "", table$Y,   fixed = TRUE)
  table$Sky       <- gsub(" ", "", table$Sky, fixed = TRUE)
  table$SNR       <- gsub(" ", "", table$SNR, fixed = TRUE)
  
  # Construct the data frame.
  df <- data.frame(
    target=target,
    star=as.character(table$Star),
    
    InstMag=table$IM,
    CatMag=table[,mag_column],
    X=as.numeric(table$X),
    Y=as.numeric(table$Y),
    SNR=as.numeric(table$SNR),
    Sky=as.numeric(table$Sky),
    
    VPhot_file=VPhot_file,
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

parse_VPhot_header_line <- function(lines, key){
  line <- lines[substring(lines,1,nchar(key))==key][1]
  trimws(substring(line,nchar(key)+1))
}
