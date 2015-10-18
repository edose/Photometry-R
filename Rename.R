##### Rename.R: Rename or renumber FITS files in the user-specified folder.
#####    Eric Dose, Bois d'Arc Observatory, Kansas, USA -- begun October 3 2015.


##### renumberACP(): Renames FITS files from ACP naming ("T Cep-S001-R001-C001-I.fts")
#####    to sequential naming ("T Cep-001-I.fts").
#####    Also returns data frame of FITS files' metadata, including new file names.
renumberACP <- function (folder="J:/Astro/Images/C11/2015/20150825/Renaming/") {
  df <- data.frame(ACP_name=list.files(path=folder,pattern="[.]*.f[[:alpha:]]*t"), 
                   target=NA, stringsAsFactors=FALSE)
  
  # Parse ACP filenames for all FITS files in this folder.
  pattern <- paste("(.+)-S[[:digit:]]{3}-R[[:digit:]]{3}-C[[:digit:]]{3}-",
                   "([[:alpha:]][[:alnum:]]*)[_]*(.*)f[[:alpha:]]+", sep="")
  filesDetected <- nrow(df)
  for (iRow in 1:filesDetected) {
    filename <- df$ACP_name[iRow]
    substrings <- regmatches(filename, regexec(pattern, filename))[[1]]
    df$target[iRow] <- substrings[2]
  }
  df <- df[!is.na(df$target),]  # keep only files with ACP-style names.
  filesToRename <- nrow(df)
  cat(filesDetected, "FITS files detected, of which", filesToRename, "shall be renamed.\n")
  if (filesToRename <= 0) {
    stop(cat("STOPPING: no FITS files to rename in folder", folder, 
             "(have you already renamed them?)"))
  }
  
  # Extract multiple metadata from FITS files' headers.  
  source('C:/Dev/Photometry/Utility.R')
  require(FITSio)
  require(dplyr)
  get_header_value <- function(header, key) {  # nested function.
    value <- header[which(header==key)+1]
    if (length(value)==0) value <- NA
    value
    }
  df <- cbind(df, Object=NA, JD_start=NA, JD_mid=NA, Filter=NA, Airmass=NA, FWHM=NA, Exp_secs=NA)
  for (iRow in 1:nrow(df)) {  
    fullfilename <- make_safe_path(folder, df$ACP_name[iRow])
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
  JD_half_duration <- (as.numeric(df$Exp_secs)/2) / (24*3600) # convert seconds to JD days.
  df$JD_mid <- as.character(as.numeric(df$JD_start) + JD_half_duration)
  
  # Report any mismatches between ACP name and FITS-header Object.
  mismatches <- df %>%
    select(ACP_name,target,Object) %>%
    filter(target!=Object)
  cat ("Object-target mismatches (s/b zero) =",nrow(mismatches), "\n")
  if (nrow(mismatches) > 0) {
    write.table(mismatches,file="",quote=FALSE,row.names=FALSE)
  }

  # Construct time-based index for files *within* each target, then construct new names.
  df2 <- df %>% 
    group_by(Object) %>% 
    arrange(JD_start) %>% 
    mutate(one=1) %>% 
    mutate(index=trimws(as.character(cumsum(one)))) %>%
    mutate(index=substring(paste("000",index,sep=""),nchar(index)+1,nchar(index)+3)) %>%
    select(-target,-one)
  df2$newName <- paste(df2$Object,"-",df2$index,"-",df2$Filter,".fts",sep="") # (dplyr mutate didn't work)
  
  # Rename all FITS files (after confirmation).
  userConfirmedRename <- "Y"==toupper(trimws(readline(
    cat("Do you really want to rename all",nrow(df2), "FITS files? (y/n)"))))
  if (userConfirmedRename) {
    df_rename <- df2 %>%
      select(ACP_name,newName) %>%
      mutate(oldName=make_safe_path(folder,ACP_name)) %>%
      mutate(newName=make_safe_path(folder,newName))
    outcome <- file.rename(df_rename$oldName, df_rename$newName)
    cat(sum(outcome),"files have been renamed...")
    write.csv(df2,file=make_safe_path(folder,"$catalog.csv"),row.names=FALSE)
    cat("and file $catalog.csv has been written.\n")
  }
  return(df2)
}
