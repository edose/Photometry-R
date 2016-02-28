##### $Archive.R: code no longer needed but possibly useful for reference.

make_VPhot_sequence_master_df <- function(sequences=c("CF Cas"), folder="C:/Dev/Photometry/VPhot/") {
  # [Archived 11/10/2015]
  # Reads a user-supplied string vector of VPhot sequence names, and
  # returns a data frame of all data for all those sequences.
  df <- data.frame()
  for (sequence in sequences) {
    df <- rbind(df,get_one_VPhot_sequence(sequence=sequence, folder=folder))
  }
  return(df)
}

renumberACP <- function (folder="J:/Astro/Images/C11/2015/20150825/Renaming/") {
  # Renames FITS files from ACP naming ("T Cep-S001-R001-C001-I.fts")
  # to sequential naming ("T Cep-001-I.fts").
  # Also returns data frame of FITS files' metadata, including new file names,
  # and writes this data to a file $catalog.csv stored with the FITS files.
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
  source('C:/Dev/Photometry/$Utility.R')
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

apertureXXX <- function(image, X0, Y0, discR=9, innerR=12, outerR=16) {
  # [was aperture(), but renamed to avoid later name conflicts]
  require(dplyr)
  Xsize <- dim(image)[1]
  Ysize <- dim(image)[2]
  testRadius  <- outerR + 1.5
  Xlow  <- max(1, floor(X0-testRadius))
  Xhigh <- min(Xsize, ceiling(X0+testRadius))
  Ylow  <- max(1, floor(Y0-testRadius))
  Yhigh <- min(Ysize, ceiling(Y0+testRadius))
  subImage <- image[Xlow:Xhigh, Ylow:Yhigh]
  discR2 <- discR^2
  innerR2 <- innerR^2
  outerR2 <- outerR^2
  
  X <- matrix((0:(ncol(subImage)-1))+Xlow, nrow=nrow(subImage), ncol=ncol(subImage), byrow=TRUE)
  Y <- matrix((0:(nrow(subImage)-1))+Ylow, nrow=nrow(subImage), ncol=ncol(subImage), byrow=FALSE)
  Xoffset <- 0  # will be a loop, later.
  Yoffset <- 0  #   ""
  dX <- X - (X0 + Xoffset)
  dY <- Y - (Y0 + Yoffset)
  dist2 <- ((dX*dX) + (dY*dY))
  discMask <- dist2 <= discR2
  annulusMask <- (dist2>=innerR2) & (dist2<=outerR2)
  blotMask <- blot(X, Y)
  discMask <- discMask & blotMask
  annulusMask <- annulusMask & blotMask
  discNA <- ifelse(discMask, subImage, NA)        # either pixel value or NA
  annulusNA <- ifelse(annulusMask, subImage, NA)  #    ""
  annulusMedian <- median(annulusNA, na.rm=TRUE)  # this could be more sophisticated later
  annulusMean   <- mean(annulusNA, na.rm=TRUE)
  discNet <- discNA - annulusMedian
  totalNetFlux <- sum(discNet, na.rm=TRUE)
  if(totalNetFlux <= 0) {stop(paste("totalNetFlux=", totalNetFlux))}
  Xcentroid <- sum(discNet * X, na.rm=TRUE) / totalNetFlux
  Ycentroid <- sum(discNet * Y, na.rm=TRUE) / totalNetFlux
  # Test here for convergence. Loop back if change > TOL, drop through if <= TOL.
  sigma <- sqrt(sum(discNet * dist2, na.rm=TRUE) / totalNetFlux) / 2
  FWHM  <- sigma * (2 * sqrt(2 * log(2)))
  maxPixel <- max(discNA, na.rm=TRUE)
  return(list(Xcentroid=Xcentroid, Ycentroid=Ycentroid, Flux=totalNetFlux,
              SkyMedian=annulusMedian, SkyMean=annulusMean,
              sigma=sigma, FWHM=FWHM, maxPixel=maxPixel))
}
