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
  if (!file.exists(FOV_path)) {
    return (NA)
  }
  lines <- readLines(FOV_path, warn=FALSE)   # read last line even without EOL character(s).
  for (iLine in 1:length(lines)) {
    lines[iLine] <- lines[iLine] %>% 
      strsplit(";",fixed=TRUE) %>% unlist() %>% first() %>% trimws()  # remove comments
  }

  # Parse DIRECTIVE LINES -> FOV_data (a list)
  directiveLines <- lines[stri_detect_regex(lines,'^#')] # detect and collect directive text lines.
  # Nested function:
  directive_value <- function(key) {
    line  <- directiveLines[substring(directiveLines,1,nchar(key))==key][1]
    directive <- line %>% strsplit("[ \t]",fixed=FALSE) %>% unlist() %>% first()
    value <- line %>% substring(nchar(directive)+1) %>% trimws()  # all but the directive
  }
  FOV_data <- list()
  
  # REQUIRED directives here.
  FOV_data$Sequence <- directive_value("#SEQUENCE")
  center <- directive_value("#CENTER") %>% strsplit("[ \t]+",fixed=FALSE) %>% unlist() %>% trimws()
  FOV_data$RA_center  <- get_RA_deg(center[1])
  FOV_data$Dec_center <- get_Dec_deg(center[2])
  FOV_data$Chart <- directive_value("#CHART")
  FOV_data$Date <- directive_value("#DATE")
  
  # OPTIONAL directives below here.
  exps <- directive_value("#EXP")
  if (is.na(exps[1])) {
    #df_exps <- NA
    FOV_data$Exposures <- NA
  } else {
    exps <- exps %>% strsplit("[ \t]+",fixed=FALSE) %>% unlist() %>% trimws()
    df_exps <- data.frame(Filter="", Seconds="", stringsAsFactors = FALSE) # primer row; removed below.
    for (iExp in 1:length(exps)) {
      strs <- exps[iExp] %>% strsplit("=",fixed=TRUE) %>% unlist() %>% trimws()
      filter <- strs[1]
      exps_filter <- strs[2] %>% strsplit(",",fixed=TRUE) %>% unlist() %>% trimws()
      for (iExpFilter in 1:length(exps_filter)) {
        thisExp <- c(Filter=filter,Seconds=exps_filter[iExpFilter])
        df_exps <- rbind(df_exps,thisExp)
      }
    }
    FOV_data$Exposures <- df_exps %>% filter(Filter!="") %>% filter(Seconds!="") # remove primer row.
  }
  cad_strs <- directive_value("#CADENCE") %>% strsplit("[ \t]+",fixed=FALSE) %>% unlist() %>% trimws()
  if(length(cad_strs) >= 2) {
    cad_num  <- cad_strs[1] %>% as.numeric()
    cad_unit <- cad_strs[2] %>% substr(1,1) %>% tolower()
    if ((cad_num <= 0) | (!cad_unit %in% c("s","m","h","d"))) {
      FOV_data$Cadence     <- NA
      FOV_data$CadenceUnit <- NA
    } else {
      FOV_data$Cadence     <- cad_num
      FOV_data$CadenceUnit <- cad_unit
    }
  } else {
    FOV_data$Cadence     <- NA
    FOV_data$CadenceUnit <- NA
  }
  stare_strs <- directive_value("#STAREFOR") %>% strsplit("[ \t]+",fixed=FALSE) %>% unlist() %>% trimws()
  if (length(stare_strs) >= 2) {
    stare_num  <- stare_strs[1] %>% as.numeric()
    stare_unit <- stare_strs[2] %>% substr(1,1) %>% tolower()
    if ((stare_num <= 0) | (!stare_unit %in% c("s","m","h","d"))) {
      FOV_data$Stare     <- NA
      FOV_data$StareUnit <- NA
    } else {
      FOV_data$Stare     <- stare_num
      FOV_data$StareUnit <- stare_unit
    }
  } else {
    FOV_data$Stare     <- NA
    FOV_data$StareUnit <- NA
  }
  ACP_strs <- directive_value("#ACP") %>% strsplit("|",fixed=TRUE) %>% unlist() %>% trimws()
  if(length(ACP_strs) < 4) {
    FOV_data$ACP <- NA
  } else {
    ACP_strs[length(ACP_strs)] <- paste0(FOV_data$Sequence,"\t", 
                                         get_RA_hours(FOV_data$RA_center),"\t", 
                                         format(FOV_data$Dec_center,nsmall=6)," ; ", 
                                         ACP_strs[length(ACP_strs)])  
    FOV_data$ACP <- paste0( paste(ACP_strs,collapse="\n"), "\n")
  }
  
  # Parse STAR LINES (embedded lines from VPhot sequence):
  starsDirectiveAt <- charmatch("#STARS", lines)
  df_star <- read.table(FOV_path,header=FALSE, sep="\t", skip=starsDirectiveAt, fill=TRUE, strip.white=TRUE, 
                        comment.char=";", stringsAsFactors = FALSE, 
                        col.names=c("StarID", "degRA", "degDec", "Mags", "StarType",6:10))
  df_star <- df_star[!stri_startswith_fixed(df_star[,1],";"),1:5] # remove comments, keep only 1st 5 cols
  if (nrow(df_star)<=0) {
    df_star <- NA  # no star data, prob if VPhot sequence has not yet been developed, e.g.,first-time FOV).
  } else {         # normal case where star data was read properly.
    df_star$MagU <- NA
    df_star$MagB <- NA
    df_star$MagV <- NA
    df_star$MagR <- NA
    df_star$MagI <- NA

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
          mag_value <- as.numeric(mag_key_value[2])
          CH_rows[irow,column_name] <- ifelse(mag_value==0, NA, mag_value) # VPhot writes zero to mean NA.
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
  
    # Diagnostic checks & messages before returning results.
    if (sum(df_star$StarType=="Check")<=0) {
      cat(">>>>> Warning: FOV file ",FOV_name," has NO CHECK STAR.\n")
    }
  }
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
  #   (full hex string), or (2) already a numeric string between -90 and +90 (degrees string).
  str <- unlist(strsplit(Dec_string,":",fixed=TRUE))
  if (length(str)==1) {
    Dec_deg = as.double(Dec_string)
  } else {
    str <- c(str,rep("0",3))
    Dec_sign <- sign(as.numeric(paste(str[1],"1",sep=""))) # append digit to prevent "-0" problem.
    Dec_deg <- Dec_sign * sum(abs(as.numeric(str[1:3])) * c(1, 1/60, 1/3600))
  }
  if (Dec_deg < -90 | Dec_deg > +90 ) Dec_deg = NA
  return(Dec_deg)
}

get_RA_hours <- function (RA_deg) {
  deg <- RA_deg %% 360
  hours <- floor(deg/15.0)
  minutes <- floor(4 * (deg-15*hours))
  seconds <- (deg-15*hours-minutes/4) * 240
  return (paste0(stri_pad_left(hours,width=2,pad="0"), ":",
                 stri_pad_left(minutes,width=2,pad="0"), ":",
                 stri_pad_left(format(seconds,nsmall=1),width=4,pad="0")))
}

read_AAVSO_Chart <- function(chartID) {
  ##### Tests OK 20151220.
  require(rvest)
  require(dplyr)
  url <- paste("https://www.aavso.org/apps/vsp/photometry/?chartid=", chartID, 
               "&Ic=on&B=on&Rc=on&maglimit=16", sep="")
  df_url <- url %>% read_html() %>% html_table() %>% nth(1)
  if (ncol(df_url)!=10) {
    stop("Table is corrupt, does not have 10 columns.")
  }
  # Convert text in table to more useful quantities.
  df <- df_url %>% 
    mutate(ChartID=chartID) %>%
    mutate(RA_deg=NA, Dec_deg=NA, Bmag=NA, Berr=NA, Vmag=NA, Verr=NA, Rmag=NA, Rerr=NA, Imag=NA, Ierr=NA)
  for (iRow in 1:nrow(df)) {
    # Unfortunately, strsplit() & unlist() don't work with %>%, thus must use loop.
    df$RA_deg[iRow]  <- strsplit(df$RA[iRow],"[",fixed=TRUE) %>% unlist() %>% first() %>% get_RA_deg()
    df$Dec_deg[iRow] <- strsplit(df$Dec[iRow],"[",fixed=TRUE) %>% unlist() %>% first() %>% get_Dec_deg()
    splits <- strsplit(df$B[iRow],"[()]") %>% unlist()
    df$Bmag[iRow] <- splits[1] %>% getNumericOrNA()
    df$Berr[iRow] <- splits[2] %>% getNumericOrNA()
    splits <- strsplit(df$V[iRow],"[()]") %>% unlist()
    df$Vmag[iRow] <- splits[1] %>% getNumericOrNA()
    df$Verr[iRow] <- splits[2] %>% getNumericOrNA()
    splits <- strsplit(df$Rc[iRow],"[()]") %>% unlist()
    df$Rmag[iRow] <- splits[1] %>% getNumericOrNA()
    df$Rerr[iRow] <- splits[2] %>% getNumericOrNA()
    splits <- strsplit(df$Ic[iRow],"[()]") %>% unlist()
    df$Imag[iRow] <- splits[1] %>% getNumericOrNA()
    df$Ierr[iRow] <- splits[2] %>% getNumericOrNA()
  }
  return(df %>% select(-RA, -Dec, -B, -V, -Rc, -Ic, -matches("B-V")))
}

getNumericOrNA <- function(string) {
  suppressWarnings(as.numeric(string)) # if not numeric, just return NA silently.
}
