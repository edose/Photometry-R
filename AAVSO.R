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