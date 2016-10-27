listEclipsers <- function(JD=2457633.75, hours_tol=4) {
  require(dplyr, quietly=TRUE)
  require(stringi)
  source("C:/Dev/Photometry/$Utility.R")
  days_tol <- hours_tol / 24.0
  FOV_files <- list_FOV_files()
  df <- data.frame(hrs_diff=NULL, FOV_name=NULL, RA_hr=NULL, Dec=NULL, Mag_V=NULL, Priority=NULL, Period=NULL,
                   Comments=NULL, stringsAsFactors = FALSE)
  
  ##### Nested function:
  add_row <- function (df, JD_diff, fov_name, fov) {
    df <- bind_rows(
      df, list(hrs_diff=round(24*JD_diff,2), FOV_name=fov_name, 
               RA=round(fov$FOV_data$RA_center/15,1), Dec=round(fov$FOV_data$Dec_center,1), 
               Mag_V=fov$FOV_data$Mag_V_bright, Priority=fov$FOV_data$Priority, 
               Period=round(fov$FOV_data$Period,3), Comments=fov$FOV_data$ACP_comments )
    )
  }
  
  
   for (filename in FOV_files) {
    fov_name <- filename %>% strsplit(".txt") %>% unlist()
    fov <- read_FOV_file(fov_name)
    # Treat this fov.
    type <- fov$FOV_data$Target_type
    if (type == "Eclipser") {
      n_exact <- (JD-fov$FOV_data$JD_faint) / fov$FOV_data$Period
      n_closest <- round(n_exact)
      if (fov$FOV_data$Period>0) {
        JD_closest <- fov$FOV_data$JD_faint + n_closest * fov$FOV_data$Period
        JD_diff <- JD_closest-JD
         if (abs(JD_diff) <= days_tol) {
          df <- add_row(df, JD_diff, fov_name, fov)
          for (i in 1:10) {
            JD_i = JD_closest - i * fov$FOV_data$Period
            JD_diff <- JD_i - JD
            if (abs(JD_diff) > days_tol) {
              break
            }
            df <- add_row(df, JD_diff, fov_name, fov)
          }
          for (i in 1:10) {
            JD_i = JD_closest + i * fov$FOV_data$Period
            JD_diff <- JD_i - JD
            if (abs(JD_diff) > days_tol) {
              break
            }
            df <- add_row(df, JD_diff, fov_name, fov)
          }
        }
      }
    }
    df <- df %>% arrange(FOV_name, hrs_diff)
   }
  #cat(">>>>> Now, 'print(df,right=FALSE)'\n")
  print(df, right=FALSE)
  df
}

list_FOV_files <- function()  {
  # Make list of filenames (as strings) of all FOVs now available. Omit those starting with $.
  FOV_input_folder = "C:/Dev/Photometry/FOV"
  require(dplyr)
  require(stringi)
  all_files <- list.files(FOV_input_folder, full.names=FALSE, recursive=FALSE, include.dirs=FALSE) %>% 
    setdiff(list.dirs("C:/Dev/Photometry/FOV", full.names=FALSE))
  df <- data.frame(FOV_files=all_files, stringsAsFactors = FALSE) %>% 
    filter(!stri_startswith_fixed(FOV_files,"$")) %>% arrange(FOV_files)
  df$FOV_files
}

list_FOV_objects <- function () {
  # Make list of all FOVs (objects) now available.
  require(dplyr)
  require(stringi)
  FOV_files <- list_FOV_files()
  objects <- c()
  for (filename in FOV_files) {
    fov_name <- filename %>% strsplit(".txt") %>% unlist()
    fov <- read_FOV_file(fov_name)
    objects <- c(objects,fov)
  }
  objects
}

extract_FOV_data <- function () {
  # Valid for FOV version 1.2. Modify as needed for desired data field.
  require(dplyr, quietly=TRUE)
  require(stringi)
  source("C:/Dev/Photometry/$Utility.R")
  FOV_files <- list_FOV_files()
  types <- c()
  for (filename in FOV_files) {
    fov_name <- filename %>% strsplit(".txt") %>% unlist()
    fov <- read_FOV_file(fov_name)
    # This section does the extraction...alter as desired.
    class <- fov$FOV_data$Target_class
    type <- fov$FOV_data$Target_type
    types <- c(types, paste0(class,"::",type))
    if (type=="Delta Scuti") {
      cat(fov_name, "::", class, type)
    }
  }
  for (type in unique(types)) {
    cat(type,"\n")
  }
}
