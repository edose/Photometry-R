##### Planning.R, gets important data from a few astro web sites, including sites with forms.
#####    Eric Dose, Bois d'Arc Observatory, Kansas, USA -- begun October 3 2015.


##### obsPlanner(): gets VS data from AAVSO Observation Planner, selected by type, timespan, etc.
#####    Returns a data frame sorted by RA.
obsPlanner <- function (VStype="%", faintMagLimit=15, localStdTime=22, maxHoursEW=2,
                         decLimitS=0, decLimitN=60, selectBest=75) {
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
  RA_deg <- 15 * (as.numeric(as.difftime(table$RA,units="hours")))
  Dec_deg <- rep(NA, nrow(table))
  pattern <- "([+-]*)([[:digit:]]+)[:]{1}([[:digit:]]+)[:]{1}([[:digit:].]+)"
  for (iRow in 1:nrow(table)) {
    substrings <- regmatches(table$Dec[iRow], regexec(pattern, table$Dec[iRow]))[[1]]  
    Dec_deg[iRow] <- sum(as.numeric(substrings[3:5]) * c(1, 1/60, 1/3600)) *
      ifelse(trimws(substrings[2])=="-", -1, 1)    
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

eveningMiras <- function (localStdTime=22, maxHoursEW=3, decLimitS=0, decLimitN=85, selectBest=80) {
  obsPlanner(VStype="M%", localStdTime=localStdTime, maxHoursEW=maxHoursEW, 
             decLimitS=decLimitS, decLimitN=decLimitN, selectBest=selectBest)
}