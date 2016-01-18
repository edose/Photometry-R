##### Plots.R
##### All diagnostic and presentation plots for project Photometry.
##### Initiated Jan 17 2016.
##### Eric Dose, Bois d'Arc Observatory, Kansas


diagnostics <- function(modelList) {
  ##### Reads modelList from modelOneFilter() run, then performs diagnostics.
  ##### Testing needed.
  ##### Typical usage: diagnostics(listV) where listV is output of modelOneFilter() for one filter (e.g. V).
  require(dplyr, quietly=TRUE)
  require(ggplot2, quietly=TRUE)
  
  model      <- modelList$model
  obs        <- modelList$obs
  AN         <- modelList$AN
  image      <- modelList$image
  star       <- modelList$star
  filter     <- modelList$filter
  transform  <- modelList$transform
  extinction <- modelList$extinction

  ##### Plots looking for outlier individual observations.
  
  # Residuals vs Fitted plot.
  x <- obs$Fitted
  xRange <- max(x) - min(x)
  y <- obs$Residual
  ptLabels <- ifelse(abs(y)>=0.05, obs$Serial, "")
  df_plot <- data.frame(x=x, y=y, ptLabels=ptLabels)
  p <- ggplot(df_plot, aes(x,y)) +
    geom_point() +
    geom_text(aes(label=ptLabels), hjust=-0.1) +
    xlim(min(x)-0.03*xRange, max(x)+0.1*xRange) +
    ggtitle(paste0("Residuals vs Fitted: ", AN, "     filter=", filter)) +
    xlab(paste0("Fitted (Mag)")) +
    ylab("Residual\n(Mag)") +
    Photometry.bw.theme()
  print(p)
  
  # Residuals vs vignette plot.
  x <- obs$Vignette
  xRange <- max(x) - min(x)
  y <- obs$Residual
  ptLabels <- ifelse(abs(y)>=0.05, obs$Serial, "")
  df_plot <- data.frame(x=x, y=y, ptLabels=ptLabels)
  p <- ggplot(df_plot, aes(x,y)) +
    geom_point() +
    geom_text(aes(label=ptLabels), hjust=-0.1) +
    xlim(min(x)-0.03*xRange, max(x)+0.1*xRange) +
    ggtitle(paste0("Residuals vs Vignette: ", AN, "     filter=", filter)) +
    xlab(paste0("Vignette (squared fraction of distance from center to corner)")) +
    ylab("Residual\n(Mag)") +
    Photometry.bw.theme()
  print(p)
  
  # Residuals vs Max ADU (uncalibrated) plot.
  library(scales)
  x <- obs$MaxADU
  xRange <- max(x) - min(x)
  y <- obs$Residual
  ptLabels <- ifelse(abs(y)>=0.05, obs$Serial, "")
  df_plot <- data.frame(x=x, y=y, ptLabels=ptLabels)
  p <- ggplot(df_plot, aes(x,y)) +
    geom_point() +
    geom_text(aes(label=ptLabels), hjust=-0.1) +
    #xlim(min(x)-0.03*xRange, max(x)+0.1*xRange) +
    scale_x_log10(breaks=c(100,200,500,1000,2000,5000,10000,20000,50000,100000)) +
    ggtitle(paste0("Residuals vs Max ADU: ", AN, "     filter=", filter)) +
    xlab(paste0("Maximum ADU in raw (uncalibrated) image [log scale]")) +
    ylab("Residual\n(Mag)") +
    Photometry.bw.theme()
  print(p)
  
  ##### Plots looking for outlier images etc.
  
  # The cirrus plot.
  JD_mid_num <- as.numeric(image$JD_mid)
  JD_mid_floor <- floor(JD_mid_num)
  x <- JD_mid_num - JD_mid_floor
  xRange <- max(x) - min(x)
  y <- image$CirrusEffect
  ptLabels <- ifelse(abs(y)>=0.02, image$FITSfile, "")
  df_plot <- data.frame(x=x, y=y, ptLabels=ptLabels)
  p <- ggplot(df_plot, aes(x,y)) +
    geom_point() +
    geom_text(aes(label=ptLabels), hjust=-0.05) +
    xlim(min(x)-0.03*xRange, max(x)+0.25*xRange) +
    ggtitle(paste0("Cirrus Plot: ", AN, "     filter=", filter)) +
    xlab(paste0("JD(mid) - ", JD_mid_floor)) +
    ylab("Cirrus Effect\n(Mag)") +
    Photometry.bw.theme()
  print(p)
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