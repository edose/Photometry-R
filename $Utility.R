##### $Utility.R, Support for VS Photometry
##### Eric Dose, Bois d'Arc Observatory, Kansas, USA -- begun September 18 2015.

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
