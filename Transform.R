##### Transform.R, Filter transform routine(s) for Photometry
##### Eric Dose, Bois d'Arc Observatory, Kansas, USA -- begun August 2015.

##### from VPHOT photometry files, construct data frame ready for subsetting, then 
#####    applying either linear model lm() or mixed-model lmer().
#####    
make_transform_df <- function (VPHOTfolder="C:\\") {
  ##### Argument "folder" must be a folder in which every .txt file is a transform VPHOT file.
  dft <- make_master_df (VPHOTfolder)
  cat(nrow(df),"rows.\n")

  ##### Add V-I color field (B-V colors are already given by VPHOT), or NA if V or I not available.
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
  # color MUST be either "V-I" or "B-V". Put proper color in field "CI".
  if (color=="B-V") { df_fit$CI <- df_fit$BV_color } 
               else { df_fit$CI <- df_fit$VI_color } 
  
  # Keep only rows with valid and appropriate data for this fit.
  df_fit <- df_fit[df_fit$Filter==filter, ]         # select rows for current filter.
  df_fit <- df_fit[!is.na(df_fit$CI), ]             # select rows having color (required for transform).
  df_fit <- df_fit[df_fit$SNR>=minSNR, ]            # only bright stars in this filter.
  df_fit <- df_fit[!(df_fit$star %in% omitStars), ] # remove stars per user.

  # Apply mixed-model lmer() if > one file for this filter, else regular lm() if only one.
  numFiles = length(unique(df_fit$VPHOT_file))  
  if (numFiles == 1) {
    model <- lm (InstMag ~ CatMag + CI, data=df_fit)  # regular linear model.
  } else {
    require(lme4)
    model <- lmer(InstMag ~ CatMag + CI + (1 | VPHOT_file), data=df_fit) # mixed-model.
  }
  model$star <- df_fit$star
  return(model)
}

plot_t <- function (model_transform) {
  plot(model_transform, labels.id=model_transform$star)
}