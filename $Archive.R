##### $Archive.R: code no longer needed but possibly useful for reference.

make_VPhot_sequence_master_df <- function(sequences=c("CF Cas"), folder="C:/Dev/Photometry/VPhot/") {
  # [Archived 11/10/2015]
  # Reads a user-supplied string vector of VPhot sequence names, and
  # returns a data frame of all data for all those sequences.
  df <- data.frame()
  for (sequence in sequences) {
    df <- rbind(df,get_one_VPhot_sequence(sequence=sequence, folder=folder))
  }
  return(df)
}