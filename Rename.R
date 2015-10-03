##### Rename.R, renames FITS files from ACP (set-repeat-count naming) to sequential naming.
##### Eric Dose, Bois d'Arc Observatory, Kansas, USA -- begun October 3 2015.

renameACP <- function (folder="J:/Astro/Images/C11/2015/20150825/Renaming/") {
  df <- data.frame(ACP_name=list.files(path=folder,pattern="[.]*.f[[:alpha:]]*t"), 
                   target=NA, stringsAsFactors=FALSE)

    # Parse ACP filenames for all FITS files in this folder.
    pattern <- paste("(.+)-S[[:digit:]]{3}-R[[:digit:]]{3}-C[[:digit:]]{3}-",
                         "([[:alpha:]][[:alnum:]]*)[_]*(.*)f[[:alpha:]]+", sep="")
    for (iRow in 1:nrow(df)) {
      filename <- df$ACP_name[iRow]
      substrings <- regmatches(filename, regexec(pattern, filename))[[1]]
      df$target[iRow] <- substrings[2]
    }
    
  # Extract multiple metadata from FITS files' headers.  
  require(FITSio)
  get_header_value <- function(header,key) {  # nested function.
    value <- header[which(header==key)+1]
    if (length(value)==0) value <- NA
    value
    }
  df <- cbind(df, Object=NA, JD_start=NA, JD_mid=NA, Filter=NA, Airmass=NA, FWHM=NA, Exp_secs=NA)
  for (iRow in 1:nrow(df)) {  
    fullfilename <- make_safe_path(folder,df$ACP_name[iRow])
    fileHandle <- file(description=fullfilename, open="rb")
    header <- parseHdr(readFITSheader(fileHandle))
    close(fileHandle)
    df$Object[iRow]   <- get_header_value(header, "OBJECT")
    df$JD_start[iRow] <- get_header_value(header, "JD")
    df$Filter[iRow]   <- get_header_value(header, "FILTER")
    df$Airmass[iRow]  <- get_header_value(header, "AIRMASS")
    df$FWHM[iRow]     <- get_header_value(header, "FWHM")
    df$Exp_secs[iRow] <- get_header_value(header, "EXPTIME")
  }
  df$JD_mid <- as.character(as.numeric(df$JD_start) + as.numeric(df$Exp_secs) / (24*3600))
  cat ("Number of target-Object mismatches (s/b zero) =",sum(df$target != df$Object), "\n")
  
  # Construct time-based index for files *within* each target
  df2 <- df %>% 
    group_by(Object) %>% 
    arrange(JD_start) %>% 
    mutate(one=1) %>% 
    mutate(index=trimws(as.character(cumsum(one)))) %>%
    mutate(index=substring(paste("000",index,sep=""),nchar(index)+1,nchar(index)+3)) %>%
    select(-target,-one)

  # Construct new file names.
  df2$newName <- paste(df2$Object,"-",df2$index,"-",df2$Filter,sep="")
  
  return(df2)
}