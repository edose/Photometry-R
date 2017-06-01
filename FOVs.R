# Thoroughly tested with FOV 1.5 on April 27 2017.

listEclipsers2 <- function(JD=2457633.75, hours_tol=4, secondaries='integrated', min_priority=0.5) {
  # secondaries expected to be one of 'integrated', 'separate', or 'none'.
  require(dplyr, quietly=TRUE)
  require(stringi)
  source("C:/Dev/Photometry/$Utility.R")
  days_tol <- hours_tol / 24.0
  FOV_files <- list_FOV_files()
  df <- data.frame(hrs_diff=NULL, FOV_name=NULL, Min_type=NULL, RA_hr=NULL, Dec=NULL, Mag_V=NULL, 
                   Priority=NULL, Period=NULL, Motive=NULL, stringsAsFactors = FALSE)
  
  ##### Nested function:
  add_row <- function (df, JD_diff, fov_name, fov, min_type) {
    df <- bind_rows(
      df, list(hrs_diff=round(24*JD_diff,2), FOV_name=fov_name, Min_type=min_type,
               RA=round(fov$FOV_data$RA_center/15,1), Dec=round(fov$FOV_data$Dec_center,1), 
               Mag_V=fov$FOV_data$Mag_V_bright, Priority=fov$FOV_data$Priority, 
               Period=round(fov$FOV_data$Period,3), Motive=fov$FOV_data$Motive )
    )
    df
  }
  ##### (End nested function.)
  
  ##### Nested function:
  add_all_rows_this_min <- function (df, JD, days_tol, fov, JD_min, min_type) {
    period = fov$FOV_data$Period
    n_exact <- (JD - JD_min) / period
    n_closest <- round(n_exact)
    JD_closest <- JD_min + n_closest * period
    JD_diff <- JD_closest - JD
    if (abs(JD_diff) <= days_tol) {
      df <- add_row(df, JD_diff, fov_name, fov, min_type)  # row for min closest to given JD (if in tol).
      for (i in 1:10) {
        JD_i = JD_closest - i * period
        JD_diff <- JD_i - JD
        if (abs(JD_diff) > days_tol) {
          break
        }
        df <- add_row(df, JD_diff, fov_name, fov, min_type)  # rows for qualifying mins before closest JD.
      }
      for (i in 1:10) {
        JD_i = JD_closest + i * period
        JD_diff <- JD_i - JD
        if (abs(JD_diff) > days_tol) {
          break
        } # if abs
        df <- add_row(df, JD_diff, fov_name, fov, min_type)  # rows for qualifying mins after closest JD.  
      } # for i
    } # if abs
    df
  } # (fn)
  ##### (End nested function.)
  
  # Build dataframe of all minima within time window (tol).
  for (filename in FOV_files) {
    fov_name <- filename %>% strsplit(".txt") %>% unlist()
    fov <- read_FOV_file(fov_name)
    # Treat this fov.
    type <- fov$FOV_data$Target_type
    if (type == "Eclipser") {
      if (fov$FOV_data$Period > 0) {
        if (fov$FOV_data$Priority >= min_priority) {
          
          # Treat PRIMARY minimum ("JD_faint"):
          df <- add_all_rows_this_min(df, JD, days_tol, fov, fov$FOV_data$JD_faint, 1)
          
          #Treat SECONDARY minimum ("JD_second") IF available in FOV and requested by user.
          if (secondaries %in% c('integrated', 'separate')) {
            if (!is.na(fov$FOV_data$JD_second)) {
              df <- add_all_rows_this_min(df, JD, days_tol, fov, fov$FOV_data$JD_second, 2)
            }
          } # if secondaries
        } # if Priority
      } # if Period
    } # if type
  } # for filename
  
  if (secondaries == 'separate') {
    df <- df %>% arrange(Min_type, FOV_name, hrs_diff)
  } else {
    df <- df %>% arrange(FOV_name, hrs_diff)  # works for both 'integrated' and 'none'.
  }
  print(df, right=FALSE)
  df
}

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
               # Period=round(fov$FOV_data$Period,3), Comments=fov$FOV_data$ACP_comments )
               Period=round(fov$FOV_data$Period,3), Motive=fov$FOV_data$Motive )
    )
  }
  ##### (End nested function.)
  
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
  # Tests OK as is Jan 9 2017 for FOV version 1.4. Modify as needed for desired data field.
  require(dplyr, quietly=TRUE)
  require(stringi)
  source("C:/Dev/Photometry/$Utility.R")
  FOV_files <- list_FOV_files()
  types <- c()
  for (filename in FOV_files) {
    fov_name <- filename %>% strsplit(".txt") %>% unlist()
    fov <- read_FOV_file(fov_name)
    
    # The section below does the extraction...alter as needed.
    type <- fov$FOV_data$Target_type
    types <- c(types, paste0(type))
    if (type=="Delta Scuti") {
      cat(fov_name, " :: ", type, "\n")
    }
  }
  cat('All target_types found:\n')
  for (type in unique(types)) {
    cat("  ", type,"\n")
  }
}

compDiag <- function(AN_top_folder="J:/Astro/Images/C14", AN_rel_folder, filename) {
  # Usage, e.g.: I
  if (is.null(AN_rel_folder)) {stop(">>>>> You must provide a AN_rel_folder, ",
                                    "e.g., AN_rel_folder='20151216'.")}
  require(dplyr, quietly=TRUE)
  require(stringi, quietly=TRUE)
  AN_folder   <- make_safe_path(AN_top_folder, AN_rel_folder)
  photometry_folder <- make_safe_path(AN_folder, "Photometry")
  path_df_master <- make_safe_path(photometry_folder, "df_master.Rdata")
  load(path_df_master)
  df <- df_master %>%
    filter(stri_startswith_fixed(FITSfile, filename)) %>%
    filter(StarType=='Comp') %>%
    filter(MaxADU_Ur < 57000) %>%
    select(FITSfile, FOV, Filter, Exposure, StarID, CatMag, InstMag, Sigma=InstMagSigma, 
           FWHM, MaxADU_Ur, Serial) %>%
    arrange(StarID) %>%
    mutate(diff = InstMag - CatMag)
  avg = mean(df$diff,na.rm=TRUE)
  stdev = sd(df$diff,na.rm=TRUE)
  df <- df %>% mutate(Z_value=round((diff-avg)/stdev,2)) %>% select(-FITSfile)
  print(df)
  cat('--- Std Dev ', 1000*stdev, ' mMag.\n\n')
}

diag_1_5 <- function(df_new=df_master, df_old=df_master_old) {
  cols_to_test <- colnames(df_master_new)
  for (this_col in cols_to_test) {
    new_col_data <- unlist(df_new[this_col])
    old_col_data <- unlist(df_old[this_col])
    # df_new[2, colname] <- NA   # for testing only; remove for production use
    # df_new[3, colname] <- 999  # for testing only; remove for production use
    # df_new[4, colname] <- NA   # for testing only; remove for production use
    # df_old[4, colname] <- NA   # for testing only; remove for production use
    same_value <- new_col_data == old_col_data  # but yields NA if either is NA or both are NA
    bothNA <- is.na(new_col_data) & is.na(old_col_data)
    cols_identical = all((same_value | bothNA) & (!(xor(is.na(new_col_data), is.na(old_col_data)))))
    # cols_identical <- all(length(new_col_data)==length(old_col_data)) && all(new_col_data==old_col_data)
    cat(paste0(this_col, " --> ", cols_identical, "\n"))
  }  
}

view_1_5 <- function(df_new=df_master, df_old=df_master_old, colname) {
  require(dplyr, quietly=TRUE)
  # df_new[2, colname] <- NA   # for testing only; remove for production use
  # df_new[3, colname] <- 999  # for testing only; remove for production use
  # df_new[4, colname] <- NA   # for testing only; remove for production use
  # df_old[4, colname] <- NA   # for testing only; remove for production use
  df <- data.frame(new=unlist(df_new[colname]), old=unlist(df_old[colname]), stringsAsFactors=FALSE)
  df$same_value <- df$new == df$old  # but yields NA if either is NA or both are NA
  df$bothNA <- is.na(df$new) & is.na(df$old)
  df$identical <- (df$same_value | df$bothNA) & (!(xor(is.na(df$new), is.na(df$old))))
  df$ModelStarID <- df_new$ModelStarID
  df$Filter <- df_new$Filter
  cat(all(unlist(df$identical)), " ", sum(unlist(!df$identical)))
  View(df)
}