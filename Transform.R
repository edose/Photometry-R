##### Transform.R, Filter transform routine(s) for Photometry
##### Eric Dose, Bois d'Arc Observatory, Kansas, USA -- begun August 2015.

##### from VPHOT photometry files, construct data frame ready for mixed-model lmer().
#####    
make_transform_df <- function (VPHOTfolder="C:\\") {
  ##### Argument "folder" must be a folder in which every .txt file is a transform VPHOT file.
  df <- make_master_df (VPHOTfolder)

  ##### Get V-I colors, make data frame including only stars with catalog V & I mags for color.
  V_rows  <- df[df$Filter=="V",]
  V_stars <- V_rows$star
  I_rows  <- df[df$Filter=="I",]
  I_stars <- I_rows$star
  color_stars <- unique(intersect(V_stars,I_stars))  # ID of all stars having catalog mags for both V & I.
  
  df_with_color <- df[df$star %in% color_stars,]      # keep only rows that can be used for color.
  V_CatMags <- V_rows[match(df_with_color$star, V_rows$star),]$CatMag
  I_CatMags <- I_rows[match(df_with_color$star, I_rows$star),]$CatMag
  df_with_color$VI_color <- V_CatMags - I_CatMags
  row.names(df_with_color) <- paste(df_with_color$star,df_with_color$Filter,sep="_")
  return(df_with_color)
}

##### from transform_df (made by make_transform_df()), run transform with options.
transform <- function (dft, filter="V", color="V-I", minSNR=30, omitStars="") {
  if (color=="B-V") { dft$CI <- dft$BV_color } 
               else { dft$CI <- dft$VI_color }

  dft <- dft[dft$Filter==filter & dft$SNR>=minSNR,] # dft holds only bright stars in target filter.
  dft <- dft[!(dft$star %in% omitStars),]        # remove stars per user.
  
  numFiles = length(unique(dft$VPHOT_file))  # mixed-model lmer() if > one file, regular lm() if only one.
  if (numFiles == 1) {
    model <- lm (InstMag ~ CatMag + CI, data=dft)  # regular linear model.
  } else {
    require(lme4)
    model <- lmer(InstMag ~ CatMag + CI + (1 | VPHOT_file), data=dft) # mixed-model.
  }
  return(model)
}