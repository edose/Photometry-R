#####  APT.R -- functions to control & extract data from APT (Cal Tech aperture photometry software).
#####      This is to replace previous use of PI (PixInsight), which has proven unreliable in its
#####         identification of source signals, and is in any case too tied to catalogs, 
#####         when all we want is fast, accurate AP of our targets.

##### APT_one_FITS(): process one FITS file through APT. 
#####    Requires calibrated FITS file. 
#####    Also requires APT source-list file (previously constructed from VPhot sequence file).
#####    Writes output file "APT-[FITS filename].txt" into FITS file's directory.
APT_one_FITS <- function (APTfolder="C:/Programs/APT/APT_v2.5.8/", 
                          FITSfolder=NULL, FITSfile=NULL) {
  
  # Construct arguments.
  outputFile <- make_safe_path(FITSfolder, paste("APT-",FITSfile,sep=""))
  APT_arguments <- c("-Duser.language=en",
                     "-Duser.region=US",
                     "-mx1024M",
                     "-jar", make_safe_path(APTfolder, "APT.jar"),
                     "-i", FITSfile,
                     "-s", sourceListFile,
                     "-o", outputFile)
  
  # Run APT to generate output file.
  errorCode <- system2("java", shQuote(APT_arguments), wait=TRUE)
  if (errorCode == 0) return(outputFile) else return("")
}
