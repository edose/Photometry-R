##### Transform.R, Filter transform routine(s) for Photometry
##### Eric Dose, Bois d'Arc Observatory, Kansas, USA -- begun August 2015.

##### read VPHOT files, construct data frame, invoke mixed-model lmer().
make_transform_df <- function (files="^NGC 7790") {
  ##### Argument files can be a vector of file names or a regex to find files.
  if (!is.character(files)) {
    cat("ERROR: Argument must be a string (for regex) or a vector of 2 or more filenames.\n")
    stop()
  }
  
  ##### Get the vector of file names.
  if (length(files)==1) {
    cat("Regex = '", files, "'\n")
    filelist <- trimws(list.files(".", pattern=files, recursive=FALSE))
  } else {
    cat("List of",length(files),"files.")
    filelist <- trimws(files)
  }
  cat(length(filelist), "files found and will be included.\n")
  if (length(filelist)<2) cat("ERROR: at least 2 files of the same field must be included.\n")
  
  ##### Read VPHOT files and construct a raw (superset, unedited) data frame.
  df <- data.frame()
  for (filename in filelist){
    df <- rbind(df, get_one_VPHOT_photometry_report(filename)) # get next raw data frame and append it.
  }
  
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

transform <- function (df, filter="V", color="V-I", minSNR=30, omitStars="") {
  dft <- df[df$Filter==filter & df$SNR>=minSNR,] # dft holds only bright stars in target filter.
  if (color=="B-V") { dft$CI <- dft$BV_color } 
               else { dft$CI <- dft$VI_color }
  
  dft <- dft[!(dft$star %in% omitStars),]
  
  numFiles = length(unique(dft$VPHOT_file))  # we use mixed-model if more than one file, lm() if only one.
  if (numFiles == 1) {
    model <- lm (InstMag ~ CatMag + CI, data=dft)  # regular linear model.
  } else {
    require(lme4)
    model <- lmer(InstMag ~ CatMag + CI + (1 | VPHOT_file), data=dft) # mixed-model.
  }
  return(model)
}