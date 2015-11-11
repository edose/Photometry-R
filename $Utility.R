##### $Utility.R, Support for VS Photometry
##### Eric Dose, Bois d'Arc Observatory, Kansas, USA -- begun September 18 2015.

##### SUPPORT FILES ONLY in this file. ##########################################

make_safe_path <- function (folder, filename, extension="") {
# Pastes correctly regardless of duplicated "/". Disregard extension if it's in filename.
  gsub("/+","/",paste(trimws(folder),"/",trimws(filename),trimws(extension),sep=""),fixed=FALSE)
}

read_FOV_file <- function (FOV_name) {
  require(stringi)
  require(dplyr)
  FOV_folder <- "C:/Dev/Photometry/FOV"
  FOV_path   <- make_safe_path(FOV_folder,trimws(FOV_name),".txt")
  lines <- readLines(FOV_path)
  for (line in lines) {
    line <- trimws(strsplit(line,";",fixed=TRUE)[1])  # remove comments, then trim white space
  }
  
  # Parse DIRECTIVE LINES -> FOV_data (a list)
  directiveLines <- lines[stri_detect_regex(lines,'^#')]
  # Nested function:
  directive_value <- function(key) {
    line  <- directiveLines[substring(directiveLines,1,nchar(key))==key][1]
    directive <- (unlist(strsplit(line,"[ \t]")))[1]
    value <- trimws(substring(line,nchar(directive)+1))  # all but the directive
  }
  FOV_data <- list()
  FOV_data$Sequence <- directive_value("#SEQUENCE")
  center <- directive_value("#CENTER")
  FOV_data$RA_center  <- get_RA_deg(center[1])
  FOV_data$Dec_center <- get_Dec_deg(center[2])
  # TODO: parse all the other FOV-file directives.
  
  # Parse STAR LINES (embedded lines from VPhot sequence):
  df_star <- read.table(FOV_path,header=FALSE, sep="\t", skip=0, fill=TRUE, strip.white=TRUE, 
                        comment.char="#", stringsAsFactors = FALSE, 
                        col.names=c("StarID", "degRA", "degDec", "Mags", "StarType",6:10))
  df_star <- df_star[-stri_startswith_fixed(df_star[,1],";"),1:5] # remove comments, keep only 1st 5 cols
  df_star$MagU <- NA
  df_star$MagB <- NA
  df_star$MagV <- NA
  df_star$MagR <- NA
  df_star$MagI <- NA
  df_star$Sequence <- FOV_data$Sequence
  
  # temporary lookup/cross-reference data frame
  mag_xref <- data.frame(passband=c("MagU","MagB","MagV","MagR","MagI"), stringsAsFactors = FALSE)
  rownames(mag_xref) <- c("1024","1","2","4","8") 
  
  # Make CH_rows (a row for each check and comp star) & Target (unknown) rows.
  CH_rows <- df_star[df_star$StarType %in% c("C","H"),]
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
  
  # Add target (unknown) rows to make new df_star data frame.
  df_star <- rbind(CH_rows, df_star[df_star$StarType=="T",])  # Add rows for Target (unknown) stars.
  df_star$StarType[df_star$StarType=="C"] <- "Comp"           # Rename types...inelegant of course
  df_star$StarType[df_star$StarType=="H"] <- "Check"          #  "
  df_star$StarType[df_star$StarType=="T"] <- "Target"         #  "
  df_star <- df_star[order(df_star$StarType),]                # Sort rows by star type.
  df_star$Mags <- NULL                                        # Remove no-longer-needed Mags column.
  return(list(FOV_data=FOV_data, star_data=df_star))
}

get_RA_deg <- function (RA_string) {
  # argument must be either (1) one string of format "12:34:56.232323" (full hex string), or
  #    (2) already a numeric string between 0 and 360 (degrees string).
  str <- unlist(strsplit(RA_string,":",fixed=TRUE))
  if (length(str)==1) {
    RA_deg = as.double(RA_string)
  } else {
    str <- c(str,rep("0",3))
    RA_deg <- 15 * sum(as.numeric(str[1:3]) * c(1, 1/60, 1/3600))
  }
  if (RA_deg < 0 | RA_deg > 360 ) RA_deg = NA
  return(RA_deg)
}

get_Dec_deg <- function (Dec_string) {
  # argument must be either (1) one string of format "+12:34:56.232323" or "-12:34:56.232323" or "0:0:0"
  #   (full hex string), or (2) already a numeric string between -90 and +90 (degrees string).  str <- unlist(strsplit(RA_string,":",fixed=TRUE))
  str <- unlist(strsplit(Dec_string,":",fixed=TRUE))
  if (length(str)==1) {
    Dec_deg = as.double(Dec_string)
  } else {
    str <- c(str,rep("0",3))
    print (str)
    Dec_sign <- sign(as.numeric(paste(str[1],"1",sep=""))) # append digit to prevent "-0" problem.
    print(Dec_sign)
    Dec_deg <- Dec_sign * sum(abs(as.numeric(str[1:3])) * c(1, 1/60, 1/3600))
  }
  if (Dec_deg < -90 | Dec_deg > +90 ) Dec_deg = NA
  return(Dec_deg)
}
