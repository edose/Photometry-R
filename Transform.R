##### Transform.R, Filter transform routine(s) for Photometry
##### Eric Dose, Bois d'Arc Observatory, Kansas, USA -- begun August 2015.

##### User functions in this file:
#####   make_transform_df(): first R step in transform process, run once per process.
#####   transform(): second R step in transform process, run once per filter.
#####   plot_t(): plots residuals from transform() so that user can omit outlier point in repeat
#####      run of transform() until satisfactory.

### The transform process (probably run only twice per year, per telescope+camera rig): 
###    1. Image in every filter: NGC 7790 or other suitable (dense) star field having a VPhot sequence.
###    2. Copy images' FITS files into a folder with only those files.
###    3. Run once: dft <- make_transform_df("C:/...yourFolderOfFITS").
###    4. For *each* filter, run (e.g. for filter R): summary(transform(dft, filter="R")).
###    5. After each transform() run, plot with plot_t(), omit outlier points, rerun transform().
###    6. The transform value for each filter (given the color index definition, e.g.: R given V-I) 
###          is the "CI" (Color Index) coefficient. 

##### Model improved Oct 21 2015 to fix CatMag coefficient to 1.

# make_transform_df(): from VPhot photometry files, construct data frame ready for subsetting, then 
#    apply either linear model lm() or mixed-model lmer().
#    Argument "folder": a folder in which every .txt file is a transform VPhot file.
make_transform_df <- function (VPhotFolder="C:\\") {
  source("C:/Dev/Photometry/VPhot.R")
  dft <- make_VPhot_master_df(VPhotFolder)
  cat(nrow(df),"rows.\n")

  ##### Add V-I color field (B-V colors are already given by VPhot), or NA if V or I not available.
  V_rows  <- df[df$Filter=="V",]
  V_stars <- V_rows$star
  I_rows  <- df[df$Filter=="I",]
  I_stars <- I_rows$star
  V_CatMags <- V_rows[match(dft$star, V_rows$star),]$CatMag
  I_CatMags <- I_rows[match(dft$star, I_rows$star),]$CatMag
  dft$VI_color <- V_CatMags - I_CatMags
  return(dft)
}

##### from transform_df (made by make_transform_df()), run transform with options.
transform <- function (dft, filter="V", color="V-I", minSNR=30, omitStars=c("")) {
  df_fit <- dft
  # color MUST be either "V-I" or "B-V". Put proper color index in field "CI".
  if (color=="B-V") { df_fit$CI <- df_fit$BV_color } 
               else { df_fit$CI <- df_fit$VI_color } 
  
  # Keep only rows with valid and appropriate data for this fit.
  df_fit <- df_fit[toupper(df_fit$Filter)==toupper(filter), ]  # select rows for current filter.
  df_fit <- df_fit[!is.na(df_fit$CI), ]             # keep only rows having color (required for transform).
  df_fit <- df_fit[df_fit$SNR>=minSNR, ]            # keep only stars with sufficient S/N ratio.
  df_fit <- df_fit[!(df_fit$star %in% omitStars), ] # remove stars per user.

  # Apply mixed-model lmer() if multiple files for this filter; if only one file, apply regular lm().
  numFiles = length(unique(df_fit$VPhot_file)) # number of files represented in df_fit.
  #cat(">", df_fit$Vphot_file, "<", "\n")
  #cat(">", df_fit$Vphot_file[1], "<", "\n")
  #cat(">", unique(df_fit$Vphot_file), "<", "\n")
  #cat(numF, "\n")
  cat("transform() using",nrow(df_fit), "rows in", numFiles, "files of filter", filter, "\n")
  cat("df_fit has",sum(is.na(df_fit)),"values=NA\n")
  if (numFiles == 1) {
    model <- lm (InstMag ~ offset(CatMag) + CI, data=df_fit)  # regular linear model.
    model$star <- df_fit$star                         # add star names for plot_t().
  } else {
    model <- lm (InstMag ~ offset(CatMag) + CI + as.factor(JD), data=df_fit)
    model$star <- df_fit$star                         # add star names for plot_t().
  }
  return(model)
}

##### Custom transform plot to label extreme points by star name.
plot_t <- function (model_transform) {
  plot(model_transform, labels.id=model_transform$star)
}