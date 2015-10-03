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
    
  # Get exposure-start Julian Date from FITS files.  
  require(FITSio)
  get_header_value <- function(header,key) {  # nested function.
    value <- header[which(header==key)+1]
    if (length(value)==0) value <- NA
    value
    }
  df <- cbind(df, JD_start=NA, Object=NA, Filter=NA, Airmass=NA, FWHM=NA, Exp_secs=NA)
  for (iRow in 1:nrow(df)) {  
    fullfilename <- make_safe_path(folder,df$ACP_name[iRow])
    cat("Opening >",fullfilename,"<\n",sep="")
    fileHandle <- file(description=fullfilename, open="rb")
    header <- parseHdr(readFITSheader(fileHandle))
    close(fileHandle)
    df$JD_start[iRow] <- get_header_value(header, "JD")
    df$Object[iRow]   <- get_header_value(header, "OBJECT")
    df$Filter[iRow]   <- get_header_value(header, "FILTER")
    df$Airmass[iRow]  <- get_header_value(header, "AIRMASS")
    df$FWHM[iRow]     <- get_header_value(header, "FWHM")
    df$Exp_secs[iRow] <- get_header_value(header, "EXPTIME")
  }  
  return(df)
}