##### VPhot.R, Input & data frames from AAVSO's VPhot photometry file format
##### Eric Dose, Bois d'Arc Observatory, Kansas, USA -- begun August 2015.

##### User functions in this file: [none]

# make_VPhot_master_df(): [support, not called by user; called by Transform.R::make_transform_df()].
#    Reads a folder of VPhot photometry-report files, aggregates to one master data frame.
#    Argument "VPhotFolder" must be a folder in which every .txt file is a VPhot photometry report file
#    intended for use in determining filter transforms.
make_VPhot_master_df <- function (VPhotFolder="C:\\") {

  filenames <- trimws(list.files(VPhotFolder, pattern=".txt$", full.names=TRUE, 
                                 recursive=FALSE, ignore.case=TRUE))
  df <- data.frame()
  for (filename in filenames){   # TODO: redo this as list-build and then do.call() on rbind().
    df <- rbind(df, get_one_VPhot_photometry_report(filename)) # get next raw data frame and append it.
  }
  return(df)
}

# count_images_in_VPhot_master_df(): [support, not called by user].
#    Returns number of images (VPhot photometry reports) represented in VPhot master data frame, 
#    or -1 if data frame appears invalid.
count_images_in_VPhot_master_df<- function(VPhot_master_df){
  n_filename <- length(unique(VPhot_master_df$VPhot_file))
  n_subset <- nrow(unique(data.frame(VPhot_master_df$VPhot_file, VPhot_master_df$JD, 
                                     VPhot_master_df$target)))
  if (n_filename != n_subset) { return(-1) } 
  else                        { return (n_subset)  }
}

# get_one_VPhot_photometry_report(): [support, not called by user].
#    Reads one tab-delimited file from VPhot photometry report, returns one R dataframe holding file's data.
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

# parse_VPhot_header_line(): [support, not called by user].
#    Locates one header line (of VPhot photometry report) that matches a key string, returns value string.
parse_VPhot_header_line <- function(lines, key){
  line <- lines[substring(lines,1,nchar(key))==key][1]
  trimws(substring(line,nchar(key)+1))
}

# get_one_VPhot_sequence(): [support, not typically called by user].
#    Reads one sequence (star list with RA/dec & magnitudes) as constructed in VPhot 
#    and saved locally as a .txt file. Returns a small data frame.
get_one_VPhot_sequence <- function(sequence=c("CF Cas"),
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

# make_VPhot_sequence_master_df(): Reads a user-supplied string vector of VPhot sequence names, and
#    returns a data frame of all data for all those sequences.
make_VPhot_sequence_master_df <- function(sequences=c("CF Cas"),
                                    folder="C:/Dev/Photometry/VPhot/") {
  df <- data.frame()
  for (sequence in sequences) {
    df <- rbind(df,get_one_VPhot_sequence(sequence=sequence, folder=folder))
  }
  return(df)
}