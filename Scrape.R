##### Scrape.R, gets important data from a few astro web sites, including sites with forms.
#####    Eric Dose, Bois d'Arc Observatory, Kansas, USA -- begun October 3 2015.

ObsPlanner <- function (VStype="EB", faintMagLimit=15, localStdTime=22, maxHoursEW=2,
                         decLimitS=0, decLimitN=60) {
  require(rvest)
  session <- html_session("https://www.aavso.org/observation-planner-tool")
  form_in <- (session %>% read_html() %>% html_form()) [[2]]
  form_out <- set_values(form_in, vartype=VStype, faintlim=trimws(as.character(faintMagLimit)), 
                       loctime=trimws(as.character(localStdTime)), maxh=trimws(as.character(maxHoursEW)), 
                         mindec=trimws(as.character(decLimitS)), maxdec=trimws(as.character(decLimitN)))
  result <- submit_form(session, form_out)
  table <- (result %>% read_html() %>% html_table()) [[2]]
  if (ncol(table) != 12) {
    stop("Table is corrupt, does not have 12 columns.")
  }
  table <- table[-1,]  # remove first row, which is only the old column names.
  colnames(table) <- c("Name", "RA", "Dec", "Type", "max", "min", "band", "period", 
                         "N_total", "N_30", "obs_30", "Days")
  
  # Parse RA and Dec (h:m:s) to degrees, add columns to table.
  RA_deg <- 15 * (as.numeric(as.difftime(table$RA,units="hours")))
  Dec_deg <- rep(NA, nrow(table))
  pattern <- "([+-]*)([[:digit:]]+)[:]{1}([[:digit:]]+)[:]{1}([[:digit:].]+)"
  for (iRow in 1:nrow(table)) {
    substrings <- regmatches(table$Dec[iRow], regexec(pattern, table$Dec[iRow]))[[1]]  
    Dec_deg[iRow] <- sum(as.numeric(substrings[3:5]) * c(1, 1/60, 1/3600)) *
      ifelse(trimws(substrings[2])=="-", -1, 1)    
  }
  table <- table %>% 
    mutate(RA_deg=RA_deg, Dec_deg=Dec_deg)
}