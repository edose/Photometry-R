##### Transform.R, Filter transform routine(s) for Photometry
##### Eric Dose, Bois d'Arc Observatory, Kansas, USA -- begun August 2015.

##### Model improved Oct 21 2015 to fix CatMag coefficient to 1.
##### Support functions (for VPhot Photometric Report files) moved here Nov 10 2015 from old VPhot.R.

### The transform process (probably run only twice per year, per telescope+camera rig): 
###    1. Image in every filter: NGC 7790 or other suitable (dense) star field having a VPhot sequence.
###    2. Calibrate images, copy their FITS files into a folder with *only* those files.
###    3. Run once: dft <- make_transform_df("C:/...yourFolderOfFITS").
###    4. For *each* filter, run (e.g. for filter R): summary(transform(dft, filter="R")).
###    5. After each transform() run, plot with plot_t(), omit outlier points, resume at step 4 if needed.
###    6. The transform value for each filter (given the color index definition, e.g.: R given V-I) 
###          is the "CI" (Color Index) coefficient. 

make_transform_df <- function (VPhotFolder="C:\\") {
# From VPhot photometry files, construct data frame ready for subsetting, then 
# apply either linear model lm() or mixed-model lmer().
# Argument "folder": a folder in which every .txt file is a transform VPhot file.
  source("C:/Dev/Photometry/VPhot.R")
  dft <- make_VPhot_master_df(VPhotFolder)
  cat(nrow(dft),"rows.\n")

  ##### Add V-I color field (B-V colors are already given by VPhot), or NA if V or I not available.
  V_rows  <- dft[dft$Filter=="V",]
  V_stars <- V_rows$star
  I_rows  <- dft[dft$Filter=="I",]
  I_stars <- I_rows$star
  V_CatMags <- V_rows[match(dft$star, V_rows$star),]$CatMag
  I_CatMags <- I_rows[match(dft$star, I_rows$star),]$CatMag
  dft$VI_color <- V_CatMags - I_CatMags
  return(dft)
}

transform <- function (dft, filter="V", color="V-I", minSNR=30, omitStars=c("")) {
# Compute transform in one filter. 
# Input is from transform_df (made by make_transform_df()).
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

plot_t <- function (model_transform) {
# Plot tranform model (one filter), to find then remove extreme points by star name.
  plot(model_transform, labels.id=model_transform$star)
}

##### SUPPORT FILES ONLY below this line. ##########################################

make_VPhot_master_df <- function (VPhotFolder="C:\\") {
  # [support, not called by user; called by make_transform_df()].
  #    Reads a folder of VPhot photometry-report files, aggregates to one master data frame.
  #    Argument "VPhotFolder" must be a folder in which every .txt file is a VPhot photometry report file
  #    intended for use in determining filter transforms.
  filenames <- trimws(list.files(VPhotFolder, pattern=".txt$", full.names=TRUE, 
                                 recursive=FALSE, ignore.case=TRUE))
  df <- data.frame()
  for (filename in filenames){
    df <- rbind(df, get_one_VPhot_photometry_report(filename)) # get next raw data frame and append it.
  }
  return(df)
}

count_images_in_VPhot_master_df<- function(VPhot_master_df){
  # [support, not called by user].
  # Returns number of images (VPhot photometry reports) represented in VPhot master data frame, 
  # or -1 if data frame appears invalid.
  n_filename <- length(unique(VPhot_master_df$VPhot_file))
  n_subset <- nrow(unique(data.frame(VPhot_master_df$VPhot_file, VPhot_master_df$JD, 
                                     VPhot_master_df$target)))
  if (n_filename != n_subset) { return(-1) } 
  else                        { return (n_subset)  }
}

get_one_VPhot_photometry_report <- function (filepathname="C:/Dev/Photometry/NGC 7790 I 1.txt"){
  # [support, not called by user].
  # Reads one tab-delimited file from VPhot photometry report, returns one R dataframe holding file's data.
  
  # Get header lines.
  lines        <- readLines(filepathname)
  target       <- parse_VPhot_header_line(lines,"Primary target:")
  exposure     <- parse_VPhot_header_line(lines,"Exposure time:")
  filter       <- parse_VPhot_header_line(lines,"Filter:")
  obs_datetime <- parse_VPhot_header_line(lines,"Observation date/time:")
  JD           <- parse_VPhot_header_line(lines,"JD:")
  decimal_date <- parse_VPhot_header_line(lines,"Decimal date:")
  RA           <- parse_VPhot_header_line(lines,"R.A.:")
  Dec          <- parse_VPhot_header_line(lines,"Dec.:")
  airmass      <- parse_VPhot_header_line(lines,"Airmass:")
  calibration  <- parse_VPhot_header_line(lines,"Calibration:")
  ap_radius    <- parse_VPhot_header_line(lines,"Apeture radius:") # yes, misspelled in AAVSO report.
  VPhot_file   <- parse_VPhot_header_line(lines,"File name:")
  
  #Get Target (unknown) lines, if any.
  table_start_key <- "Star\tIM\tSNR\t"
  table_start_line <- which(substring(lines,1,nchar(table_start_key))==table_start_key)[1]
  table <- read.table(filepathname,skip=table_start_line-1,header=TRUE,sep="\t")
  cat("   ", filepathname,"table has",nrow(table),"rows.\n")
  
  # Get index of column with catalog magnitudes by finding ".mag" at end of column's name.
  colnamelengths <- nchar(colnames(table))
  mag_column <- which(substring(colnames(table),colnamelengths-3,colnamelengths)==".mag")[1]
  
  # Some fields as read have units at strings' end--remove them.
  exposure  <- trimws(gsub(" s",     "",exposure, fixed=TRUE)) # Remove " s" from end of exposure.
  ap_radius <- trimws(gsub(" pixels","",ap_radius,fixed=TRUE)) # Remove " pixels" from end of ap radius.
  
  # Some fields have spaces for thousands separators (inside the strings)--remove them.
  exposure        <- gsub(" ", "", exposure,  fixed = TRUE)
  JD              <- gsub(" ", "", JD,        fixed = TRUE)
  airmass         <- gsub(" ", "", airmass,   fixed = TRUE)
  ap_radius       <- gsub(" ", "", ap_radius, fixed = TRUE)  
  table$X         <- gsub(" ", "", table$X,   fixed = TRUE)
  table$Y         <- gsub(" ", "", table$Y,   fixed = TRUE)
  table$Sky       <- gsub(" ", "", table$Sky, fixed = TRUE)
  table$SNR       <- gsub(" ", "", table$SNR, fixed = TRUE)
  
  # Construct the data frame.
  df <- data.frame(
    target=target,
    star=as.character(table$Star),
    
    InstMag=table$IM,
    CatMag=table[,mag_column],
    X=as.numeric(table$X),
    Y=as.numeric(table$Y),
    SNR=as.numeric(table$SNR),
    Sky=as.numeric(table$Sky),
    
    VPhot_file=VPhot_file,
    Filter=filter,
    Exposure=as.numeric(exposure),
    Airmass=as.numeric(airmass),
    JD=as.numeric(JD),
    Obs_datetime=obs_datetime,
    Dec_date=decimal_date,
    RA=RA,
    Dec=Dec,
    Cal=calibration,
    Ap_radius=as.numeric(ap_radius),
    BV_color=as.numeric(table$B.V),
    
    stringsAsFactors=FALSE)
}

parse_VPhot_header_line <- function(lines, key){
  # [support, not called by user].
  # Locates one header line (of VPhot photometry report) that matches a key string, returns value string.
  line <- lines[substring(lines,1,nchar(key))==key][1]
  trimws(substring(line,nchar(key)+1))
}

get_one_VPhot_sequence <- function(sequence=c("CF Cas"), folder="C:/Dev/Photometry/VPhot/") {
  # [support, not typically called by user].
  # Reads one sequence (star list with RA/dec & magnitudes) as constructed in VPhot 
  # and saved locally as a .txt file. Returns a small data frame.
  filepathname <- make_safe_path(folder, sequence, ".txt")
  df <- read.table(filepathname, header=FALSE, sep="\t", skip=0, fill=TRUE, strip.white=TRUE, 
                   stringsAsFactors = FALSE)
  df <- df[,1:5]
  colnames(df) <- c("StarID", "degRA", "degDec", "Mags", "StarType")
  df$MagU <- NA
  df$MagB <- NA
  df$MagV <- NA
  df$MagR <- NA
  df$MagI <- NA
  df$Sequence <- trimws(sequence)
  mag_xref <- data.frame(passband=c("MagU","MagB","MagV","MagR","MagI"), stringsAsFactors = FALSE)
  rownames(mag_xref) <- c("1024","1","2","4","8") # temporary lookup data frame
  
  # Build a row for each check and comp star.
  CH_rows <- df[df$StarType %in% c("C","H"),]
  for (irow in 1:nrow(CH_rows)) {
    mag_entries <- trimws(unlist(strsplit(CH_rows$Mags[irow],  "|", fixed=TRUE)))
    for (mag in mag_entries) {
      mag_key_value <- trimws(unlist(strsplit(mag,"_",fixed=TRUE)))
      column_name   <- mag_xref[mag_key_value[1],]
      if (!is.na(column_name)) {
        CH_rows[irow,column_name] <- as.numeric(mag_key_value[2])
      }
    }
  }
  curated_df <- CH_rows                                     # Start with rows for Check and Comp stars.
  curated_df <- rbind(curated_df, df[df$StarType=="T",])    # Add rows for Target (unknown) stars.
  curated_df$StarType[curated_df$StarType=="C"] <- "Comp"   # Rename types...inelegant of course
  curated_df$StarType[curated_df$StarType=="H"] <- "Check"  #  "
  curated_df$StarType[curated_df$StarType=="T"] <- "Target" #  "
  curated_df <- curated_df[order(curated_df$StarType),]     # Sort rows by star type.
  curated_df$Mags <- NULL
  return(curated_df)
}

