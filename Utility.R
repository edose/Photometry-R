##### Utility.R, Support for VS Photometry
##### Eric Dose, Bois d'Arc Observatory, Kansas, USA -- begun September 18 2015.

##### User functions in this file: [none]

# make_safe_path(): Pastes correctly regardless of duplicated "/". Disregard extension if it's in filename.
make_safe_path <- function (folder, filename, extension="") {
  gsub("/+","/",paste(trimws(folder),"/",trimws(filename),trimws(extension),sep=""),fixed=FALSE)
}