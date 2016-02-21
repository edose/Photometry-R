##### Aperture.R  Supporting Aperture and sky-background funcions, perhaps to replace APT.
#####    In development February 2016.
#####    NOTE: In FITS image matrices: X is *FIRST* index, Y is SECOND index.

makeRawAperture <- function (image, Xcenter, Ycenter, Rdisc=8, Rinner=11, Router=16) {
  require(dplyr)
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
  subImage <- image[Xlow:Xhigh, Ylow:Yhigh] # FITS image matrices: X is *FIRST* index, Y is SECOND index.
  X <- matrix((0:(ncol(subImage)-1))+Xlow, nrow=nrow(subImage), ncol=ncol(subImage), byrow=TRUE)
  Y <- matrix((0:(nrow(subImage)-1))+Ylow, nrow=nrow(subImage), ncol=ncol(subImage), byrow=FALSE)
  dX <- X - Xcenter
  dY <- Y - Ycenter
  dist2 <- ((dX*dX) + (dY*dY))
  discMask <- ifelse(dist2 <= R2disc, 1, 0)                         # 1=include pixel, else 0
  skyMask  <- ifelse((dist2 >= R2inner) & (dist2 <= R2outer), 1, 0) # 1=include pixel, else 0
  return (list(subImage=subImage, Xlow=Xlow, Ylow=Ylow, discMask=discMask, skyMask=skyMask,
               discPixels=sum(discMask), skyPixels=sum(skyMask)))
}

evalAperture <- function (aperture, evalSkyFunction=NULL) {
  discADU <- aperture$subImage * aperture$discMask
  if (is.null(evalSkyFunction)) {
    skyADU <- weighted.mean(aperture$subImage, w=aperture$skyMask) # mean of sky ADUs
  } else {
    skyADU <- evalSkyFunction(aperture)
  }
  netADU <- ifelse(aperture$discMask>0, discADU - skyADU, NA)
  netFlux <- sum(netADU, na.rm=TRUE)
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
  maxADU <- ifelse(aperture$discMask>0, aperture$subImage, NA) %>% max(na.rm=TRUE)
  return(list(Xcentroid=Xcentroid, Ycentroid=Ycentroid, netFlux=netFlux, skyADU=skyADU,
              sigma=sigma, FWHM=FWHM, maxADU=maxADU, 
              discPixels=aperture$discPixels, skyPixels=aperture$skyPixels))
}


################################################################################
##### Support-only functions not normally called by user. ######################

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

evalSky005 <- function (aperture) {
  # trimmed-mean of 12 slice trimmed-means
  # This is the algorithm that won the test tournament of Feb 21 2016.
  sliceADUs <- make_sky_slices(aperture, nSlices=12, method="trimmedMean")
  skyADU    <- mean(sliceADUs, trim=0.3)
  return (skyADU) 
}

make_sky_slices <- function (aperture, nSlices=8, method="trimmedMean") {
  require(dplyr)
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
      # plotImage(sliceADUs, Xlow, Ylow)
    }
  }
  return(sliceList %>% unlist())
}

punch <- function (aperture, df_punch=NULL, Rpunch=7) {
  # df_punch: columns X & Y only.
  if (is.null(df_punch)) return (aperture)
  if (nrow(df_punch)<=0 | Rpunch<=0) return (aperture)
  require(dplyr)
  nCol <- ncol(aperture$subImage)
  nRow <- nrow(aperture$subImage)
  Xlow     <- aperture$Xlow
  Ylow     <- aperture$Ylow
  Xcenter  <- aperture$Xlow + (1+ncol(aperture$subImage))/2
  Ycenter  <- aperture$Ylow + (1+nrow(aperture$subImage))/2
  skyMask  <- aperture$skyMask
  R2punch <- Rpunch^2
  RmaxPotentialPunch <-sqrt((Xlow-Xcenter)^2 + (Ylow-Ycenter)^2) + Rpunch + 2
  X <- matrix((0:(nCol-1))+Xlow, nrow=nRow, ncol=nCol, byrow=FALSE)
  Y <- matrix((0:(nRow-1))+Ylow, nrow=nRow, ncol=nCol, byrow=TRUE)
  for (iPunch in 1:nrow(df_punch)) {
    R2punchToCenter <- (df_punch$X[iPunch]-Xcenter)^2 + (df_punch$Y[iPunch]-Ycenter)^2
    if (R2punchToCenter <= RmaxPotentialPunch^2) {
      dX <- X - df_punch$X[iPunch]
      dY <- Y - df_punch$Y[iPunch]
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


################################################################################
##### Test frameworks only. ####################################################

plotImage <- function (image, Xlow, Ylow) {
  # image needs to be a matrix of values.
  # Note: In FITS image matrices: X is *FIRST* index, Y is SECOND index.
  df <- expand.grid(X=Xlow:(Xlow+ncol(image)-1), Y=Ylow:(Ylow+nrow(image)-1)) %>%
    mutate(Z=as.vector(image))
  p <- ggplot(df, aes(X,Y)) + geom_raster(aes(fill=Z)) + scale_y_reverse()
  print (p)
}

plotAperture <- function (aperture) {
  # aperture is an "aperture" object.
  plotImage (aperture$subImage, aperture$Xlow, aperture$Ylow)
  plotImage (aperture$discMask+0.8*aperture$skyMask, aperture$Xlow, aperture$Ylow)
}

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

