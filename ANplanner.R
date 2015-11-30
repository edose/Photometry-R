##### ANplanner.R, Astronight planning tools.
#####    Eric Dose, Bois d'Arc Observatory, Kansas, USA -- begun October 3 2015.

##### USER FUNCTIONS in this file:
#####    obsPlanner: get VS data from AAVSO Observation Planner; to rationally select a night's targets.
#####    eveningMiras: wrapper for obsPlanner() for evening Miras.

obsPlanner <- function (VStype="%", faintMagLimit=15, localStdTime=22, maxHoursEW=2,
                        decLimitS=0, decLimitN=60, selectBest=75) {
# Gets VS data from AAVSO Observation Planner, selected by type, timespan, etc.
# Returns a data frame sorted by RA. 
# Parameter "selectBest" gives the maximum number of targets, the best currently assumed to be 
#    the maximum number of previous AAVSO observations; if selectBest <= 0, return all targets found.
  require(rvest)
  require(dplyr)
  session <- html_session("https://www.aavso.org/observation-planner-tool")
  form_in <- (session %>% read_html() %>% html_form()) [[2]]
  form_out <- set_values(form_in, vartype=VStype, faintlim=trimws(as.character(faintMagLimit)), 
                         loctime=trimws(as.character(localStdTime)), maxh=trimws(as.character(maxHoursEW)), 
                         mindec=trimws(as.character(decLimitS)), maxdec=trimws(as.character(decLimitN)),
                         docsv=trimws(as.character("0")))
  result <- submit_form(session, form_out)
  table <- (result %>% read_html() %>% html_table()) [[2]]
  if (ncol(table) != 12) {
    stop("Table is corrupt, does not have 12 columns.")
  }
  table <- table[-1,]  # remove first row, which is only the old column names.
  colnames(table) <- c("Name", "RA", "Dec", "Type", "max", "min", "band", "period", 
                       "N_total", "N_30", "obs_30", "Days")
  if (nrow(table)<=0) {
    stop("ERROR: No data retrieved. Returning table with zero rows.\n")
  }
  
  # Parse RA and Dec (h:m:s) to degrees, add columns to table.
  source ("C:/Dev/Photometry/$Utility.R")
  RA_Deg  <- rep(NA, nrow(table))
  Dec_Deg <- rep(NA, nrow(table))
  for (iRow in 1:nrow(table)) {
    RA_deg[iRow]  <- get_RA_deg(table$RA[iRow])
    Dec_deg[iRow] <- get_Dec_deg(table$Dec[iRow])
  }
  table <- table %>% mutate(RA_deg=RA_deg, Dec_deg=Dec_deg)  
  
  # Convert strings representing numbers to numbers.
  table <- table %>% mutate(Days=ifelse(Days=="30+",NA,Days))  # N_30==0 already implies that Days>30.
  for (name in c("max", "min", "period", "N_total", "N_30", "obs_30", "Days")) {
    table[name] <- as.numeric(unlist(table[name]))
  }
  if (selectBest <= 0) {
    cat("Returning (all)",nrow(table),"VS targets of",length(unique(table$Type)),"distinct types.\n")
  } else {
    nrowOld <- nrow(table)
    typesOld <- length(unique(table$Type))
    table <- table %>% 
      filter(N_total>0) %>% 
      arrange(desc(N_total)) %>% 
      head(selectBest)
    cat("Returning", nrow(table), "VS targets (from", nrowOld, ") of", 
        length(unique(table$Type)), "(from", typesOld, ") distinct types.\n")
  }
  table <- table %>% arrange(RA_deg)
  return(table)
}

eveningMiras <- function (localStdTime=22, maxHoursEW=3, decLimitS=0, decLimitN=85, selectBest=150) {
  require(dplyr)
  obsPlanner(VStype="M%", localStdTime=localStdTime, maxHoursEW=maxHoursEW, 
             decLimitS=decLimitS, decLimitN=decLimitN, selectBest=selectBest) %>%
    select(-RA_deg, -Dec_deg)
}