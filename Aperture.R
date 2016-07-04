##### Aperture.R  Supporting Aperture and sky-background funcions, no longer using APT at all.
#####    TESTS OK March 12 2016.
#####    NOTE: In FITS image matrices: X is *FIRST* index, Y is SECOND index.

addPunch <- function (FOV_name=NULL, starID, RA, Dec, punchRA, punchDec, RADecAsStrings=FALSE,
                      userMustConfirm=TRUE) {
  # If RADecAsStrings==TRUE, all RA and Dec are character scalars in form "+nn:nn:nn.nn".
  require(stringi, quietly=TRUE)
  if (is.null(FOV_name)) {stop(">>>>> FOV file name is NULL.")}
  FOV_folder <- "C:/Dev/Photometry/FOV"
  FOV_path   <- make_safe_path(FOV_folder,trimws(FOV_name),".txt")
  if (!file.exists(FOV_path)) {stop(">>>>> Cannot open FOV file.")}
  lines <- readLines(FOV_path, warn=FALSE)   # read last line even without EOL character(s).
  starsDirectiveAt <- charmatch("#STARS", lines)
  if(starsDirectiveAt <= 0) {stop(">>>>> Cannot find #STARS directive in FOV file.")}
  starID <- trimws(starID)
  
  # Check that given starID is in the FOV's list of stars.
  df_FOV_stars <- (read_FOV_file(FOV_name))$star_data # data frame of stars in FOV file
  if (!starID %in% df_FOV_stars$StarID) {stop(paste0(">>>>> StarID '", starID, "' not in FOV star list."))}

  # Convert to degrees if necessary.
  if (RADecAsStrings) {
    RA       <- get_RA_deg(RA)
    Dec      <- get_Dec_deg(Dec)
    punchRA  <- get_RA_deg(punchRA)
    punchDec <- get_Dec_deg(punchDec)
  }
  
  # Make line to insert into FOV file.
  dNorth <- 3600 * ((punchDec) - (Dec))                       # in arcsec
  dEast  <- 3600 * ((punchRA)   - (RA)) * cos((pi/180)*(Dec)) # in arcsec
  punchLine <- paste0("#PUNCH ", starID, " : ",
                    (round(dNorth,2) %>% format(nsmall=2) %>% stri_pad_left(width=8)),
                    (round(dEast,2)  %>% format(nsmall=2) %>% stri_pad_left(width=8)),
                    "   ;   dNorth dEast of punch center vs target star, in arcsec")

  # Confirm if required, then write into FOV file.
  if (userMustConfirm) {
    cat(paste0("FOV '", FOV_name, "':\n     ", punchLine, "\n"))
    OKtoWrite <- "Y" == (cat("Proceed? (y/n)") %>% readline() %>% trimws() %>% toupper())
  } else {
    OKtoWrite <- TRUE
  }
  if (OKtoWrite) {
    cat(c(lines[1:starsDirectiveAt-1], punchLine, lines[starsDirectiveAt:length(lines)]), 
        file=FOV_path, sep="\n")
    if (userMustConfirm) {cat("Line inserted.\n")}
  }
}

addPunchesFromText <- function(textPath="C:/Dev/Photometry/Punches.txt", delim="\t") {
  # the file at textPath is a delimited (tab-delimited default) file of Punches,
  #    most easily generated in Excel and re-saved as tab-delimited text.
  # Header line and data lines must contain: FOV, Target, RA, Dec. 
  # blank line between target groups are helpful but not strictly required.
  
  require(dplyr, quietly = TRUE)
  df_in <- read.table(textPath, header=TRUE, sep="\t", stringsAsFactors=FALSE, strip.white=TRUE,
                      comment.char=";") %>%
    select(FOV, Target, RA, Dec)  
  df_punches <- data.frame()
  targetRA  <- NA
  targetDec <- NA
  for (i in 1:nrow(df_in)) {
    if (df_in$RA[i]=="") { # skip blank lines
      next
    }
    lineFOV    <- df_in$FOV[i]
    lineTarget <- as.character(df_in$Target[i])  
    if (is.na(lineTarget)) {lineTarget <- ""}  # because R may have tried to read Target as an integer
    lineRA     <- df_in$RA[i]
    lineDec    <- df_in$Dec[i]
    if (nchar(lineFOV) >= 1) {             # update target data on (new) target line
      targetFOV  <- lineFOV
      target     <- lineTarget
      targetRA   <- get_RA_deg(lineRA)
      targetDec  <- get_Dec_deg(lineDec)
      next
    }
    if (lineFOV=="" & lineTarget=="" & nchar(lineRA)>=1 & nchar(lineDec)>=1) { # punch line
      punchRA  <- get_RA_deg(lineRA)    # in degrees
      punchDec <- get_Dec_deg(lineDec)  # in degrees
      distance <- 3600 * distanceRADec(punchRA, punchDec, targetRA, targetDec)     # in arcsec
      dNorth <- round(3600 * (punchDec - targetDec),2)                             # in arcsec
      dEast  <- round((3600 * (punchRA  - targetRA) * cos((pi/180)*targetDec)) ,2) # in arcsec
      # distance <- round(sqrt(dNorth^2 + dEast^2),2) # in arcsec
      df_punches <- df_punches %>% 
        rbind(data.frame(i=i, FOV=targetFOV, Target=target, 
                         TargetRA=targetRA, TargetDec=targetDec,        # in degrees
                         PunchRA=punchRA, PunchDec=punchDec,            # in degrees
                         DNorth=dNorth, DEast=dEast, Distance=distance, # in arcsec
                         stringsAsFactors=FALSE))
      next  
    }  
    stop(paste0(">>>>> addPunchesFromText() can't handle input line ", i, ": ",
                lineFOV, ", ", lineTarget, ", ", lineRA, ", ", lineDec))
  }
  
  # --> BEFORE writing into FOV files, verify here that in master data frame df_punches: 
  #    (1) each FOV entry actually has a FOV file, 
  #    (2) each punch target is present in FOV file, and
  #    (3) each punch target RA & Dec is close to that in FOV file.
  areAllOK <- TRUE  # default value, to be falsified by problem.
  FOVs <- df_punches %>% select(FOV) %>% unique() %>% unlist()
  for (FOV_name in FOVs) {
    FOV_list <- read_FOV_file(FOV_name)
    # Verify FOV file exists for targets with this FOV name.
    if (is.na(FOV_list[[1]][[1]])) {
      cat(paste0(">>>>> FOV file '", FOV_name, "' does not exist.\n"))
      areAllOK = FALSE
      next
    }
    FOV_stars <- FOV_list$star_data$StarID
    
    df_thisFOV <- df_punches %>% filter(FOV==FOV_name)
    
    # Verify each target for this FOV is present and proper.
    for (i_target in 1:nrow(df_thisFOV)) {
      target    <- df_thisFOV$Target[i_target]
      targetRA  <- df_thisFOV$TargetRA[i_target]
      targetDec <- df_thisFOV$TargetDec[i_target]
      
      # Verify target's RA,Dec position in punch file is very close to that in FOV file.
      i_FOV <- match(target, FOV_stars)
      if (is.na(i_FOV)) {
        cat(paste0(">>>>> Target ", target, " is MISSING from FOV file '", FOV_name), "\n")
        areAllOK = FALSE
        next
      }

      distance <- distanceRADec(targetRA, targetDec, 
                                FOV_list$star_data$degRA[i_FOV], FOV_list$star_data$degDec[i_FOV])
      if (distance > 5/3600) { # 5 arcseconds (even this may be way too generous).
        cat(paste0(">>>>> Target ", target, " in FOV '", FOV_name, 
                   "' has punch list RA,Dec=(", targetRA, ",", targetDec, ") which does not match ",
                   "FOV-file RA,Dec=(", FOV_list$star_data$degRA[i], ",", 
                   FOV_list$star_data$degDec[i], ")"), "\n")
        areAllOK = FALSE
      }
    }
  }

  # Give any warnings, confirm that we're really going to modify FOV files.
  minPunchDistance <- 2  # in arcsec
  numTooClose <- sum((df_punches$Distance) < minPunchDistance)
  if (numTooClose >= 1) {
    cat(paste(">>>>> WARNING:", numTooClose, "Target-Punch pairs are too close to be real.\n"))
    df_punches %>% 
      filter(Distance < minPunchDistance) %>%
      select(FOV, Target, Distance) %>% 
      print()
  }
  cat(paste("Largest distance (Target-Punch) = ", max(df_punches$Distance),"arcsec\n"))
  recommendYes <- (areAllOK) & (numTooClose == 0) & (max(df_punches$Distance) <= 20)
  
  answerYES <- "Y" == (cat("Proceed? Recommend", ifelse(recommendYes, "Y", "NO!!!"), "(y/n)") %>% 
                         readline() %>% trimws() %>% toupper())
  if (answerYES) {
    # Write all the #PUNCH lines into the FOV files.
    cat(paste("[Write all", nrow(df_punches), "#PUNCH directives into", length(FOVs), "FOV files.]\n"))
    for (i in 1:nrow(df_punches)) {
      cat(paste("#PUNCH", df_punches$FOV[i], df_punches$Target[i], 
                round(df_punches$TargetRA[i],6), round(df_punches$TargetDec[i],6),
                round(df_punches$PunchRA[i],6), round(df_punches$PunchDec[i],6) ),"\n")
      addPunch(df_punches$FOV[i], df_punches$Target[i], 
               round(df_punches$TargetRA[i],6), round(df_punches$TargetDec[i],6),
               round(df_punches$PunchRA[i],6), round(df_punches$PunchDec[i],6),
               userMustConfirm = FALSE)
      }
    cat(paste("[All", nrow(df_punches), "#PUNCH directives written into", length(FOVs), "FOV files.]\n"))
  } else {
    {cat("[No change to FOV files.]\n")}
  }
  return(df_punches)  # return for testing...
}


################################################################################
##### Support-only functions not normally called by user. ######################

makeRawAperture <- function (image, Xcenter, Ycenter, Rdisc=10, Rinner=15, Router=20) {
  require(dplyr, quietly=TRUE)
  Xsize <- dim(image)[1]
  Ysize <- dim(image)[2]
  testRadius  <- Router + 1.5
  Xlow  <- max(1, floor(Xcenter-testRadius))
  Xhigh <- min(Xsize, ceiling(Xcenter+testRadius))
  Ylow  <- max(1, floor(Ycenter-testRadius))
  Yhigh <- min(Ysize, ceiling(Ycenter+testRadius))
  R2disc <- Rdisc^2
  R2inner <- Rinner^2
  R2outer <- Router^2
  # cat(paste("XYsize=", Xsize, Ysize, "\n"))
  # cat(paste("   X,Y=", Xcenter, Ycenter, "Xlims=", Xlow, Xhigh, "Ylims=", Ylow, Yhigh, "\n"))
  subImage <- image[Xlow:Xhigh, Ylow:Yhigh] # FITS image matrices: X is *FIRST* index, Y is SECOND index.
  X <- matrix((0:(ncol(subImage)-1))+Xlow, nrow=nrow(subImage), ncol=ncol(subImage), byrow=TRUE)
  Y <- matrix((0:(nrow(subImage)-1))+Ylow, nrow=nrow(subImage), ncol=ncol(subImage), byrow=FALSE)
  dX <- X - Xcenter
  dY <- Y - Ycenter
  dist2 <- ((dX*dX) + (dY*dY))
  discMask <- ifelse(dist2 <= R2disc, 1, 0)                         # 1=include pixel, else 0
  skyMask  <- ifelse((dist2 >= R2inner) & (dist2 <= R2outer), 1, 0) # 1=include pixel, else 0
  return (list(subImage=subImage, Xlow=Xlow, Ylow=Ylow, Xcenter=Xcenter, Ycenter=Ycenter,
               discMask=discMask, skyMask=skyMask,
               discPixels=sum(discMask), skyPixels=sum(skyMask)))
}

punch <- function (aperture, df_punch=NULL, Rpunch=7) {
  # df_punch: only columns X & Y are needed.
  if (is.null(df_punch)) return (aperture)
  if (nrow(df_punch)<=0 | Rpunch<=0) return (aperture)
  require(dplyr, quietly=TRUE)
  nCol <- ncol(aperture$subImage)
  nRow <- nrow(aperture$subImage)
  Xlow     <- aperture$Xlow
  Ylow     <- aperture$Ylow
  XsubImageCenter <- Xlow + (nCol-1)/2
  YsubImageCenter <- Ylow + (nRow-1)/2
  Xcenter  <- aperture$Xcenter   # centroid of main star
  Ycenter  <- aperture$Ycenter   # centroid of main star
  skyMask  <- aperture$skyMask
  R2punch <- Rpunch^2
  RmaxPotentialPunch <-sqrt((Xlow-XsubImageCenter)^2 + (Ylow-YsubImageCenter)^2) + Rpunch + 2
  X <- matrix((0:(nCol-1))+Xlow, nrow=nRow, ncol=nCol, byrow=FALSE)
  Y <- matrix((0:(nRow-1))+Ylow, nrow=nRow, ncol=nCol, byrow=TRUE)
  for (iPunch in 1:nrow(df_punch)) {
    R2punchToCenter <- df_punch$dX[iPunch]^2 + df_punch$dY[iPunch]^2
    if (R2punchToCenter <= RmaxPotentialPunch^2) {
      dX <- X - (Xcenter + df_punch$dX[iPunch])
      dY <- Y - (Ycenter + df_punch$dY[iPunch])
      dist2 <- dX*dX + dY*dY
      punchMask <- ifelse(dist2 <= R2punch, 0, 1)  # boolean, 0 to exclude from skyMask, else 1
      skyMask <- skyMask * punchMask               # do the punch, i.e., exclude pixels
    }
  }
  # Modify aperture's sky mask and its # pixels ONLY.
  aperture$skyMask   <- skyMask    
  aperture$skyPixels <- sum(aperture$skyMask)
  return (aperture)
}

evalAperture <- function (aperture, evalSkyFunction=evalSky005) {
  require(dplyr, quietly=TRUE)
  maxADU  <- max(aperture$subImage[aperture$discMask>0])
  discADU <- aperture$subImage * aperture$discMask
  skyADU <- evalSkyFunction(aperture)
  netADU <- ifelse(aperture$discMask>0, discADU - skyADU, NA)
  netFlux <- sum(netADU, na.rm=TRUE)
  if(netFlux <= 0) {netFlux = NA}  # return NA if star undetected
  gain <- 1.57 # e-/ADU, this value is for current scope BOREA
  skySigma <- sd( ifelse(aperture$skyMask==0,NA,aperture$subImage), na.rm=TRUE)
  # netFluxSigma equation after APT paper, PASP 124, 737 (2012), but pi/2 in 3rd term set to 1 as
  #   pi/2 seems hard to justify, and as 1 gives S/N closer to VPhot's values.
  netFluxSigma <- sqrt( 
    (netFlux/gain) + 
    (aperture$discPixels*(skySigma^2)) + 
    ( (1) * ((aperture$discPixels*skySigma)^2) /aperture$skyPixels ) # the 1 here was pi/2 in paper.
  )
  nCol <- ncol(aperture$subImage)
  nRow <- nrow(aperture$subImage)
  Xlow     <- aperture$Xlow
  Ylow     <- aperture$Ylow
  X <- matrix((0:(nCol-1))+Xlow, nrow=nRow, ncol=nCol, byrow=TRUE)
  Y <- matrix((0:(nRow-1))+Ylow, nrow=nRow, ncol=nCol, byrow=FALSE)
  Xcentroid <- sum(X*netADU, na.rm=TRUE) / netFlux
  Ycentroid <- sum(Y*netADU, na.rm=TRUE) / netFlux
  dX <- X - Xcentroid
  dY <- Y - Ycentroid
  dist2 <- (dX*dX + dY*dY)
  wtSumDist2 <- sum(dist2*netADU, na.rm=TRUE) / netFlux
  sigma <- sqrt(max(0,wtSumDist2)/2)
  FWHM  <- sigma * (2 * sqrt(2 * log(2)))
  return(list(Xcentroid=Xcentroid, Ycentroid=Ycentroid, maxADU=maxADU,
              netFlux=netFlux, netFluxSigma=netFluxSigma, skyADU=skyADU, skySigma=skySigma,
              FWHM=FWHM, discPixels=aperture$discPixels, skyPixels=aperture$skyPixels))
}

getImageMatrix <- function (FITS_path) {
  #####    NOTE: In FITS image matrices: X is *FIRST* index, Y is SECOND index.
  require(FITSio, quietly=TRUE)
  zz <- file(description=FITS_path, open="rb")
  header <- readFITSheader(zz)
  D <- readFITSarray(zz, header)
  close(zz)
  image <- D$imDat
  return (image)
}

getXYfromWCS <- function (df_RADec, FITS_path=NULL) {
  require(dplyr,  quietly=TRUE)
  require(FITSio, quietly=TRUE)
  
  hdr <- getFITSheaderValues(FITS_path, c("CD1_1", "CD1_2", "CD2_1", "CD2_2",
                                          "CRVAL1", "CRVAL2", "CRPIX1", "CRPIX2"))
  cd11   <- as.numeric(hdr$CD1_1)
  cd12   <- as.numeric(hdr$CD1_2)
  cd21   <- as.numeric(hdr$CD2_1)
  cd22   <- as.numeric(hdr$CD2_2)
  crval1 <- as.numeric(hdr$CRVAL1)
  crval2 <- as.numeric(hdr$CRVAL2)
  crpix1 <- as.numeric(hdr$CRPIX1)
  crpix2 <- as.numeric(hdr$CRPIX2)
  df_XY <- df_RADec %>% select(Number, StarID)
  
  dRA  <- df_RADec$degRA  - crval1
  dDec <- df_RADec$degDec - crval2
  degEW <- dRA * cos((pi/180)*df_RADec$degDec)
  degNS <- dDec
  a <- cd22 / cd12
  dX <- (degNS-degEW*a)  / (cd21-cd11*a)
  dY <- (degEW-cd11*dX)  / cd12
  df_XY$Xcentroid <- crpix1 + dX
  df_XY$Ycentroid <- crpix2 + dY
  return(df_XY)
}

evalSky005 <- function (aperture) {
  # trimmed-mean of 12 slice trimmed-means
  # This is the algorithm that won the test tournament of Feb 21 2016.
  sliceADUs <- make_sky_slices(aperture, nSlices=12, method="trimmedMean")
  skyADU    <- mean(sliceADUs, trim=0.3)
  return (skyADU) 
}

make_sky_slices <- function (aperture, nSlices=8, method="trimmedMean") {
  require(dplyr, quietly=TRUE)
  minPixelsPerSlice <- (aperture$skyPixels / nSlices) / 2
  skyMask <- aperture$skyMask
  nCol <- ncol(aperture$subImage)
  nRow <- nrow(aperture$subImage)
  Xlow <- aperture$Xlow
  Ylow <- aperture$Ylow
  X <- matrix((0:(nCol-1))+Xlow, nrow=nRow, ncol=nCol, byrow=TRUE)
  Y <- matrix((0:(nRow-1))+Ylow, nrow=nRow, ncol=nCol, byrow=FALSE)
  Xcenter <- Xlow + (nCol-1)/2
  Ycenter <- Ylow + (nRow-1)/2
  dX <- X - Xcenter
  dY <- Y - Ycenter
  arcTan <- atan2(dY,dX)
  arcTanLow <- ((2 * pi) * (0:(nSlices)) / nSlices) - pi #atan2() ranges -pi to pi
  sliceList <- list()
  for (i in 1:nSlices) {
    sliceMask <- skyMask * ifelse(between(arcTan,arcTanLow[i],arcTanLow[i+1]),1,0)
    slicePixels <- sum(sliceMask)
    if (slicePixels > minPixelsPerSlice) {
      sliceADUs <- ifelse(sliceMask==0,NA,aperture$subImage)
      if (method=="trimmedMean") {
        sliceList[[i]] <- mean(sliceADUs, trim=0.4, na.rm=TRUE)
      } else {
        sliceList[[i]] <- median(sliceADUs, na.rm=TRUE)
      }
      i <- 1
      # plotImage(sliceADUs, Xlow, Ylow)
    }
  }
  return(sliceList %>% unlist())
}


################################################################################
##### Test frameworks only. ####################################################

plotImage <- function (image, Xlow, Ylow, title="") {
  # image needs to be a matrix of values.
  # Note: In FITS image matrices: X is *FIRST* index, Y is SECOND index.
  require(ggplot2, quietly=TRUE)
  require(dplyr, quietly=TRUE)
  df <- expand.grid(X=Xlow:(Xlow+ncol(image)-1), Y=Ylow:(Ylow+nrow(image)-1)) %>%
    mutate(Z=as.vector(image))
  p <- ggplot(df, aes(X,Y)) + geom_raster(aes(fill=Z)) + 
    scale_fill_gradientn(colours=c("#111111", "#EEEEEE")) +
    scale_y_reverse() + ggtitle(title)
  print (p)
}

plotAperture <- function (aperture, title="") {
  # aperture is an "aperture" object.
  require(ggplot2, quietly=TRUE)
  require(dplyr, quietly=TRUE)
  apMod <- aperture
  minADU <- min(apMod$subImage)
  apMod$subImage <- (apMod$subImage - (minADU-1)) ^ 0.2
  plotImage (apMod$subImage, apMod$Xlow, apMod$Ylow, title)
  plotImage (aperture$discMask+0.8*aperture$skyMask, aperture$Xlow, aperture$Ylow, title)
}

#----------------------------------------------------------------------------------------------------------
plotAperturePresentation <- function (aperture, title="") {
  # aperture is an "aperture" object.
  require(ggplot2, quietly=TRUE)
  require(dplyr, quietly=TRUE)
  apMod <- aperture
  minADU <- min(apMod$subImage)
  apMod$subImage <- (apMod$subImage - (minADU-1)) ^ 0.2
  plotImagePresentation (apMod$subImage, apMod$Xlow, apMod$Ylow, title)
  plotImagePresentation (aperture$discMask+0.8*aperture$skyMask, aperture$Xlow, aperture$Ylow, title)
  plotSlicesPresentation(aperture$skyMask, apMod$Xlow, apMod$Ylow, title)
}

# I never figured out how to tile the slices in different colors.
plotSlicesPresentation <- function (image, Xlow, Ylow, title="") {
  # image needs to be a matrix of values, typically passed in from plotAperturePresentation().
  # Note: In FITS image matrices: X is *FIRST* index, Y is SECOND index.
  require(ggplot2, quietly=TRUE)
  require(dplyr, quietly=TRUE)
  nCol = ncol(image)
  nRow = nrow(image)
  df <- expand.grid(X=Xlow:(Xlow+ncol(image)-1), Y=Ylow:(Ylow+nrow(image)-1)) %>%
    mutate(Z=as.vector(image))
  Xcenter <- Xlow + (ncol(image)-1)/2
  Ycenter <- Ylow + (nrow(image)-1)/2
  X <- matrix((0:(nCol-1))+Xlow, nrow=nRow, ncol=nCol, byrow=TRUE)
  Y <- matrix((0:(nRow-1))+Ylow, nrow=nRow, ncol=nCol, byrow=FALSE)
  dX <- X - Xcenter
  dY <- Y - Ycenter
  mask <- ifelse(image==0, NA, image)
  arcTan <- atan2(dY,dX) * mask # to mask for the skyMask only.
  nSlices = 12
  arcTanPerSlice <- (2 * pi) / nSlices  # as atan2() ranges -pi to pi
  arcTanLowest <- -pi
  df$Islice = 1 + floor((as.vector(arcTan) - arcTanLowest) / arcTanPerSlice)
  df$Islice = pmax(1, df$Islice)
  df$Islice = pmin(nSlices, df$Islice)
  df$Islice = as.factor(df$Islice)
  p <- ggplot(df, aes(X,Y)) + geom_tile(aes(fill=Islice)) +
    # scale_fill_gradientn(colours=these_colors) +
    scale_fill_brewer(palette="Paired") +
    scale_y_reverse() + ggtitle(title)
  p <- p + SAS2016.bw.theme()
  p <- p + theme(panel.background = element_rect(fill = "gray10"),
                 panel.grid.major = element_line(colour="grey30"),
                 panel.grid.minor = element_line(colour="grey20"))
  print (p)
}

plotImagePresentation <- function (image, Xlow, Ylow, title="", 
                                   these_colors=c("#111111", "#EEEEEE")) {
  # image needs to be a matrix of values, typically passed in from plotAperturePresentation().
  # Note: In FITS image matrices: X is *FIRST* index, Y is SECOND index.
  require(ggplot2, quietly=TRUE)
  require(dplyr, quietly=TRUE)
  df <- expand.grid(X=Xlow:(Xlow+ncol(image)-1), Y=Ylow:(Ylow+nrow(image)-1)) %>%
    mutate(Z=as.vector(image))
  p <- ggplot(df, aes(X,Y)) + geom_raster(aes(fill=Z)) + 
    scale_fill_gradientn(colours=these_colors) +
    scale_y_reverse() + ggtitle(title)
  p <- p + SAS2016.bw.theme()
  print (p)
}

#----------------------------------------------------------------------------------------------------------

makeMockAperture <- function (Xcentroid=303.112, Ycentroid=201.8333, FWHM=5) {
  image <- matrix(data=123, nrow=2*1024, ncol=3*1024)
  sigma <- FWHM / (2 * sqrt(2 * log(2)))
  amplitude <- 328
  for (iRow in 150:250) {
    for (iCol in 250:350) {
      d2 <- ((iRow-Ycentroid)^2 / (2*sigma^2)) + ((iCol-Xcentroid)^2 / (2*sigma^2))
      image[iRow,iCol] <- image[iRow,iCol] + amplitude * exp(-d2)
    }
  }
  return(makeRawAperture(image, Xcentroid, Ycentroid))
}

evalSky000 <- function (aperture) {
  # mean
  maskedImage <- ifelse(aperture$skyMask>0, aperture$subImage, NA)
  skyADU <- mean(maskedImage, na.rm=TRUE)
  return (skyADU) 
}

evalSky001 <- function (aperture) {
  # median
  maskedImage <- ifelse(aperture$skyMask>0, aperture$subImage, NA)
  skyADU <- median(maskedImage, na.rm=TRUE)
  return (skyADU) 
}

evalSky002 <- function (aperture) {
  # trimmed mean (nearly a median)
  maskedImage <- ifelse(aperture$skyMask>0, aperture$subImage, NA)
  skyADU <- mean(maskedImage, trim=0.45, na.rm=TRUE)
  return (skyADU) 
}

evalSky003 <- function (aperture) {
  # trim away >3 sigma above, then trimmed mean
  maskedImage <- ifelse(aperture$skyMask>0, aperture$subImage, NA)
  thisMean <- mean(maskedImage, na.rm=TRUE)
  thisSD   <- sd  (maskedImage, na.rm=TRUE)
  maskedImage <- ifelse(maskedImage <= thisMean+3*thisSD, maskedImage, NA)
  skyADU <- mean(maskedImage, trim=0.45, na.rm=TRUE)
  return (skyADU) 
}

evalSky004 <- function (aperture) {
  # median of 12 slice medians
  sliceADUs <- make_sky_slices(aperture, nSlices=12, method="median")
  skyADU    <- median(sliceADUs)
  return (skyADU) 
}

evalSky006 <- function (aperture) {
  # median of 12 outlier-trimmed slice trimmed-means
  sliceADUs <- make_sky_slices(aperture, nSlices=12, method="trimmedMean")
  upperLimitADU <- mean(sliceADUs) + 3*sd(sliceADUs)
  sliceADUs <- sliceADUs[sliceADUs <= upperLimitADU]
  skyADU    <- mean(sliceADUs, trim=0.25)
  return (skyADU) 
}

evalSky007 <- function (aperture) {
  # trimmed-mean of 8 (not 12) slice trimmed-means
  sliceADUs <- make_sky_slices(aperture, nSlices=8, method="trimmedMean")
  skyADU    <- mean(sliceADUs, trim=0.3)
  return (skyADU) 
}

evalSky008 <- function (aperture) {
  # trimmed-mean of 16 (not 12) slice trimmed-means
  sliceADUs <- make_sky_slices(aperture, nSlices=16, method="trimmedMean")
  skyADU    <- mean(sliceADUs, trim=0.3)
  return (skyADU) 
}

testEvalAperture <- function(evalSkyFunction=NULL, Rdisc=8, Rinner=11, Router=16) {
  require(dplyr)
  image     <- getImageMatrix("J:/Astro/Images/C14/20151016/Calibrated/V1804 Cyg-0007-V.fts")
  df_aptest <- read.table("C:/Dev/Photometry/df_aptest.txt", header=TRUE) %>% select(-N)
  df_punch  <- df_aptest %>% select(X=Xpunch, Y=Ypunch)
  df_output <- matrix(0, nrow=nrow(df_aptest), ncol=3) %>% as.data.frame()
  colnames(df_output) <- c("netFluxClear","netFluxInterf","netFluxInterfPunched")
  
  for (i in 1:nrow(df_aptest)) {
    df_output$netFluxClear[i] <- 
      (makeRawAperture(image, df_aptest$Xclear[i], df_aptest$Yclear[i], Rdisc=8, Rinner=11, Router=16) %>% 
         evalAperture(evalSkyFunction))$netFlux
    df_output$netFluxInterf[i] <- 
      (makeRawAperture(image, df_aptest$Xinterf[i], df_aptest$Yinterf[i], Rdisc=8, Rinner=11, Router=16) %>% 
         evalAperture(evalSkyFunction))$netFlux
    df_output$netFluxInterfPunched[i] <- 
      (makeRawAperture(image, df_aptest$Xinterf[i], df_aptest$Yinterf[i], Rdisc=8, Rinner=11, Router=16) %>% 
         punch(df_punch) %>%
         evalAperture(evalSkyFunction))$netFlux
  }
  cat("mean: ", paste0(df_output %>% select(netFluxClear, netFluxInterf, netFluxInterfPunched) %>% 
                         apply(2,mean) %>% round(digits=0), collapse="   "),"\n")
  cat("  sd: ", paste0(df_output %>% select(netFluxClear, netFluxInterf, netFluxInterfPunched) %>% 
                         apply(2,sd) %>% round(digits=0), collapse="   "),"\n")
  return (cbind(df_aptest, df_output))
}

