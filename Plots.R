##### Plots.R
##### All diagnostic and presentation plots for project Photometry.
##### Initiated Jan 17 2016.
##### Eric Dose, Bois d'Arc Observatory, Kansas

eclipserPlot <- function(df=df_predictions, starID=NULL, filter=NULL) {
  # Reversed mag scale to put minima at bottom as customary 20160813.
  require(dplyr, quietly=TRUE)
  if (is.null(df) | is.null(starID)){ stop("Please give non-NULL parms.")}
  if (!is.null(filter)) {
    df <- df %>% filter(Filter==filter)
  }
  df_eclipser <- (df %>% filter(StarID==starID) %>% arrange(JD_num))
  x <- df_eclipser$JD_num
  y <- df_eclipser$TransformedMag
  JD_mid_floor <-df_eclipser$JD_mid %>% as.numeric() %>% min() %>% floor()
  y_top <- min(y) - 0.08 * (max(y)-min(y))
  y_base <- max(y) + 0.08 * (max(y)-min(y))
  plot(x,y, ylim=c(y_base, y_top),
       main=starID,
       xlab=paste0("JD(mid) - ", JD_mid_floor),
       ylab="best Mag")
}

modelPlots <- function(modelList) {
  ##### Reads modelList from modelOneFilter() run, then performs diagnostics.
  ##### Testing needed.
  ##### Typical usage: modelPlots(listV) where listV is output of modelOneFilter() for one filter (e.g. V).
  require(dplyr, quietly=TRUE)
  # require(lme4, quietly=TRUE)
  require(ggplot2, quietly=TRUE)
  
  model      <- modelList$model
  obs        <- modelList$obs
  AN         <- modelList$AN
  image      <- modelList$image
  star       <- modelList$star
  filter     <- modelList$filter
  transform  <- modelList$transform
  extinction <- modelList$extinction
  
  plotList <- list()
  sigmaResidual <- sd(obs$Residual)
  
  ##### Plots looking for outlier individual observations.
  
  # Residuals vs Fitted plot.
  x <- obs$Fitted
  xRange <- max(x) - min(x)
  y <- obs$Residual
  ptLabels <- ifelse(abs(y)>=2*sigmaResidual, obs$Serial, "")
  df_plot <- data.frame(x=x, y=y, ptLabels=ptLabels)
  p <- ggplot(df_plot, aes(x,y*1000)) +
    geom_point() +
    geom_text(aes(label=ptLabels), hjust=-0.1) +
    xlim(min(x)-0.03*xRange, max(x)+0.1*xRange) +
    ggtitle(paste0("Residuals vs Fitted: ", AN, "     ", filter, " filter ")) +
    xlab(paste0("Fitted (Mag)")) +
    ylab("Residual\n(mMag)") +
    Photometry.bw.theme()
  plotList$ResidualsFitted <- p
  
  # Residuals vs Instrument Magnitude plot.
  x <- obs$InstMag
  xRange <- max(x) - min(x)
  y <- obs$Residual
  ptLabels <- ifelse(abs(y)>=2*sigmaResidual, obs$Serial, "")
  df_plot <- data.frame(x=x, y=y, ptLabels=ptLabels)
  p <- ggplot(df_plot, aes(x,y*1000)) +
    geom_point() +
    geom_text(aes(label=ptLabels), hjust=-0.1) +
    xlim(min(x)-0.03*xRange, max(x)+0.1*xRange) +
    ggtitle(paste0("Residuals vs Instrument Magnitude: ", AN, "     ", filter, " filter ")) +
    xlab(paste0("Instrument Magnitude (1-second basis)")) +
    ylab("Residual\n(mMag)") +
    Photometry.bw.theme()
  plotList$ResidualsInstMag <- p
  
  # Q-Q plot.
  y=obs$Residual
  ptLabels <- obs$Serial # get this now, before sorted on Residuals
  df_plot <- data.frame(y=y, ptLabels=ptLabels) %>% arrange(y)
  prob <- ((1:length(y)-0.5)/length(y))
  x <- qnorm(prob)
  df_plot <- df_plot %>%
    mutate(x=x) %>%
    mutate(ptLabels=ifelse((abs(y)>=2*sd(y))|(abs(x)>=2), paste0("  ",ptLabels), ""))
  xRange <- max(x) - min(x)
  p <- ggplot(data = df_plot, aes(x=x, y=y*1000)) + 
    geom_point() + 
    geom_text(aes(hjust=0.1,label=ptLabels,angle=-25)) +
    xlim(min(x)-0.03*xRange, max(x)+0.05*xRange) +
    ggtitle(paste0("Q-Q Plot of Residuals: ", AN, "     ", filter, " filter ")) +
    geom_abline(slope=sd(y*1000)) +
    xlab(paste0("t  (sigma.residuals = ", format(sd(1000*y),nsmall=1,digits=3), " mMag)")) +
    ylab("Residual\n(mMag)") +
    annotate("text", x=mean(range(df_plot$x)), y=+Inf, label=paste0(length(x)," observations in model."), 
             hjust=+0.5, vjust=+2) +
    Photometry.bw.theme()
  plotList$QQ <- p
  
  # Residuals vs Sky Background plot.
  x <- obs$SkyADU
  xRange <- max(x) - min(x)
  y <- obs$Residual
  ptLabels <- ifelse(abs(y)>=2*sigmaResidual, obs$Serial, "")
  df_plot <- data.frame(x=x, y=y*1000, ptLabels=ptLabels)
  p <- ggplot(df_plot, aes(x,y)) +
    geom_point() +
    geom_text(aes(label=ptLabels), hjust=-0.1) +
    xlim(min(x)-0.03*xRange, max(x)+0.1*xRange) +
    ggtitle(paste0("Residuals vs Sky Background: ", AN, "     ", filter, " filter ")) +
    xlab(paste0("Sky Background (median, ADUs)")) +
    ylab("Residual\n(mMag)") +
    Photometry.bw.theme()
  plotList$ResidualsSky <- p
  

  # Residuals vs JD (time) plot.
  JD_mid_num <- as.numeric(obs$JD_mid) # not image$JD_mid...this is a plot by observations.
  JD_mid_floor <- floor(min(JD_mid_num))
  x <- JD_mid_num - JD_mid_floor
  xRange <- max(x) - min(x)
  y <- obs$Residual
  ptLabels <- ifelse(abs(y)>=2*sigmaResidual, obs$Serial, "")
  df_plot <- data.frame(x=x, y=y, ptLabels=ptLabels)
  p <- ggplot(df_plot, aes(x,y*1000)) +
    geom_point() +
    geom_text(aes(label=ptLabels), hjust=-0.1) +
    xlim(min(x)-0.03*xRange, max(x)+0.1*xRange) +
    ggtitle(paste0("Residuals vs JD(mid): ", AN, "     ", filter, " filter ")) +
    xlab(paste0(paste0("JD(mid) - ", JD_mid_floor))) +
    ylab("Residual\n(mMag)") +
    Photometry.bw.theme()
  plotList$ResidualsJD <- p
  
  # Sky Background vs JD (time) plot.
  JD_mid_num <- as.numeric(obs$JD_mid) # not image$JD_mid...this is a plot by observations.
  JD_mid_floor <- floor(min(JD_mid_num))
  x <- JD_mid_num - JD_mid_floor
  xRange <- max(x) - min(x)
  y <- obs$SkyADU
  sigmaSky <- (y-mean(y))/sd(y)
  ptLabels <- ifelse(abs(sigmaSky)>=2, obs$Serial, "")
  df_plot <- data.frame(x=x, y=y, ptLabels=ptLabels)
  p <- ggplot(df_plot, aes(x,y)) +
    geom_point() +
    geom_text(aes(label=ptLabels), hjust=-0.25) +
    xlim(min(x)-0.03*xRange, max(x)+0.1*xRange) +
    ggtitle(paste0("Sky Background vs JD(mid): ", AN, "     ", filter, " filter ")) +
    xlab(paste0(paste0("JD(mid) - ", JD_mid_floor))) +
    ylab("Sky Background\n(median ADUs)") +
    Photometry.bw.theme()
  plotList$SkyJD <- p
  
  # Residuals vs X-pixel plot.
  x <- obs$X1024*1024
  xRange <- max(x) - min(x)
  y <- obs$Residual
  ptLabels <- ifelse(abs(y)>=2*sigmaResidual, obs$Serial, "")
  df_plot <- data.frame(x=x, y=y, ptLabels=ptLabels)
  p <- ggplot(df_plot, aes(x,y*1000)) +
    geom_point() +
    geom_text(aes(label=ptLabels), hjust=-0.1) +
    xlim(min(x)-0.03*xRange, max(x)+0.1*xRange) +
    ggtitle(paste0("Residuals vs X: ", AN, "     ", filter, " filter ")) +
    xlab(paste0("X (pixels from image center)")) +
    ylab("Residual\n(mMag)") +
    Photometry.bw.theme()
  plotList$ResidualsX <- p
  
  # Residuals vs Y-pixel plot.
  x <- obs$Y1024*1024
  xRange <- max(x) - min(x)
  y <- obs$Residual
  ptLabels <- ifelse(abs(y)>=2*sigmaResidual, obs$Serial, "")
  df_plot <- data.frame(x=x, y=y, ptLabels=ptLabels)
  p <- ggplot(df_plot, aes(x,y*1000)) +
    geom_point() +
    geom_text(aes(label=ptLabels), hjust=-0.1) +
    xlim(min(x)-0.03*xRange, max(x)+0.1*xRange) +
    ggtitle(paste0("Residuals vs Y: ", AN, "     ", filter, " filter ")) +
    xlab(paste0("Y (pixels from image center)")) +
    ylab("Residual\n(mMag)") +
    Photometry.bw.theme()
  plotList$ResidualsY <- p
  
  # Residuals vs vignette plot.
  x <- sqrt(obs$Vignette)*1024
  xRange <- max(x) - min(x)
  y <- obs$Residual
  ptLabels <- ifelse((x>0.8*1024)|(abs(y)>=2*sigmaResidual), obs$Serial, "")
  df_plot <- data.frame(x=x, y=y, ptLabels=ptLabels)
  p <- ggplot(df_plot, aes(x,y*1000)) +
    geom_point() +
    geom_text(aes(label=ptLabels), hjust=-0.1) +
    xlim(min(x)-0.03*xRange, max(x)+0.1*xRange) +
    ggtitle(paste0("Residuals vs Vignette: ", AN, "     ", filter, " filter ")) +
    xlab(paste0("Vignette (pixels from image center)")) +
    ylab("Residual\n(mMag)") +
    Photometry.bw.theme()
  plotList$ResidualsVignette <- p
  
  # Residuals vs CI (color index) plot.
  x <- obs$CI
  xRange <- max(x) - min(x)
  y <- obs$Residual
  ptLabels <- ifelse(abs(y)>=2*sigmaResidual, obs$Serial, "")
  df_plot <- data.frame(x=x, y=y, ptLabels=ptLabels)
  p <- ggplot(df_plot, aes(x,y*1000)) +
    geom_point() +
    geom_text(aes(label=ptLabels), hjust=-0.1) +
    xlim(min(x)-0.03*xRange, max(x)+0.1*xRange) +
    ggtitle(paste0("Residuals vs Color Index: ", AN, "     ", filter, " filter ")) +
    xlab(paste0("Color Index (V-I)")) +
    ylab("Residual\n(mMag)") +
    Photometry.bw.theme()
  plotList$ResidualsCI <- p
  
  # Residuals vs Airmass plot.
  x <- obs$Airmass
  xRange <- max(x) - min(x)
  y <- obs$Residual
  ptLabels <- ifelse(abs(y)>=2*sigmaResidual, obs$Serial, "")
  df_plot <- data.frame(x=x, y=y, ptLabels=ptLabels)
  p <- ggplot(df_plot, aes(x,y*1000)) +
    geom_point() +
    geom_text(aes(label=ptLabels), hjust=-0.1) +
    xlim(min(x)-0.03*xRange, max(x)+0.1*xRange) +
    ggtitle(paste0("Residuals vs Airmass: ", AN, "     ", filter, " filter ")) +
    xlab(paste0("Airmass")) +
    ylab("Residual\n(mMag)") +
    Photometry.bw.theme()
  plotList$ResidualsAirmass <- p
  
  # Residuals vs Exposure time plot.
  x <- obs$Exposure
  xRange <- max(x) - min(x)
  y <- obs$Residual
  ptLabels <- ifelse(abs(y)>=2*sigmaResidual, obs$Serial, "")
  df_plot <- data.frame(x=x, y=y, ptLabels=ptLabels)
  p <- ggplot(df_plot, aes(x,y*1000)) +
    geom_point() +
    geom_text(aes(label=ptLabels), hjust=-0.1) +
    xlim(min(x)-0.03*xRange, max(x)+0.1*xRange) +
    ggtitle(paste0("Residuals vs Exposure Time: ", AN, "     ", filter, " filter ")) +
    xlab(paste0("Exposure Time (seconds)")) +
    ylab("Residual\n(mMag)") +
    Photometry.bw.theme()
  plotList$ResidualsExposure <- p
  
  # Residuals vs Max ADU (uncalibrated) plot.
  library(scales)
  x <- obs$MaxADU
  # x <- obs$MaxADU_Ur
  xRange <- max(x) - min(x)
  y <- obs$Residual
  ptLabels <- ifelse(abs(y)>=2*sigmaResidual, obs$Serial, "")
  df_plot <- data.frame(x=x, y=y, ptLabels=ptLabels)
  p <- ggplot(df_plot, aes(x,y*1000)) +
    geom_point() +
    geom_text(aes(label=ptLabels), hjust=-0.1) +
    #xlim(min(x)-0.03*xRange, max(x)+0.1*xRange) +
    scale_x_log10(breaks=c(100,200,500,1000,2000,5000,10000,20000,50000,100000)) +
    ggtitle(paste0("Residuals vs Max ADU: ", AN, "     ", filter, " filter ")) +
    xlab(paste0("Maximum ADU in raw (uncalibrated) image [log scale]")) +
    ylab("Residual\n(mMag)") +
    Photometry.bw.theme()
  plotList$ResidualsMaxADU_Ur <- p
  
  ##### Plots looking for outlier images etc.
  
  # The cirrus plot.
  JD_mid_num <- as.numeric(image$JD_mid)
  JD_mid_floor <- floor(min(JD_mid_num))
  x <- JD_mid_num - JD_mid_floor
  xRange <- max(x) - min(x)
  y <- image$CirrusEffect
  ptLabels <- ifelse(abs(y)>=2*sd(image$CirrusEffect), image$FITSfile, "")
  df_plot <- data.frame(x=x, y=y, ptLabels=ptLabels)
  p <- ggplot(df_plot, aes(x,y*1000)) +
    geom_point() +
    geom_text(aes(label=ptLabels), hjust=-0.05) +
    xlim(min(x)-0.03*xRange, max(x)+0.25*xRange) +
    ggtitle(paste0("Image Cirrus Plot: ", AN, "     ", filter, " filter ")) +
    xlab(paste0("JD(mid) - ", JD_mid_floor, "   of images")) +
    ylab("Cirrus Effect\n(mMag)") +
    annotate("text", x=mean(range(df_plot$x)), y=+Inf, label=paste0(length(x)," images in model."), 
             hjust=-0, vjust=+2) +
    Photometry.bw.theme()
  plotList$Cirrus <- p
  
  ##### Now print, in desired order, the plot objects.
  # Find outlier observations & model bias.
  print(plotList$ResidualsMaxADU_Ur)
  print(plotList$ResidualsFitted)
  print(plotList$ResidualsInstMag)
  print(plotList$ResidualsExposure)
  
  # Find additional sources of problems.
  print(plotList$ResidualsSky)
  print(plotList$ResidualsJD)
  print(plotList$SkyJD)

  # Check for non-linearity of included model parameters.
  print(plotList$ResidualsX)
  print(plotList$ResidualsY)
  print(plotList$ResidualsVignette)
  print(plotList$ResidualsCI)
  print(plotList$ResidualsAirmass)
  
  # Find outlier images.
  print(plotList$QQ)
  print(plotList$Cirrus)
}


################################################################################################
##### Below are test or support-only functions, rarely or not typically called by user. ########

Photometry.bw.theme <- function(){ # use later for publication. (but first: copy elements from SAS grey theme)
  th <- theme_bw()
  text_color = "gray20"
  th <- th + theme(plot.title=element_text(size=rel(1.6),face="bold",color=text_color))
  th <- th + theme(axis.text=element_text(size=rel(1.4),color=text_color))  # both axes.
  th <- th + theme(axis.title.x=element_text(angle=0,size=rel(1.4),face="italic",color=text_color))
  th <- th + theme(axis.title.y=element_text(angle=0,size=rel(1.4),face="italic",color=text_color))
  th <- th + theme(legend.position="none")
}