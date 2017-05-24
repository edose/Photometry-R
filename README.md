# Photometry-R
Eric Dose's R scripts etc for telescope photometry, especially for processing FITS image files and field-of-view data, all the way to AAVSO reports for submission.

## Status of this repository
The code in this repo still serves for my photometric data reduction and reporting to the AAVSO. But all code is being ported to python. I expect to use python only and to mothball all R code by about July 2017.

Why? Especially when R (with dplyr extensions) is clearly superior to python as a language*per se*, and when the RStudio IDE is superb. Well, python has 3 advantages for complex code, which this repo is becoming:

 1. **Unit testing.** Python has it, R doesn't--simple as that. I thought my photometry workflow might stay simple enough to do without this. Ha. Now I can have real confidence in at least the low-level methods and functions.
 2. **Object-orientation / real classes**: Python does have some real stupidities in its design--lack of private variables being first among them. You might as well build your home in the middle of a freeway and let everyone drive through at full speed. But if one must live with python's baby-talk mentality, at least classes lend some separation and focus of attention.
 3. **Execution speed**: Not the most important, but nice to have. R is blindingly fast for what it's designed for, but sometimes pretty slow for general programming.

That out of the way...June 2017. I'll make public (here) an early draft of the python code, repo 'photrix'. 

## My workflow
Wow, could it be a year since I last updated this README? Yes, almost to the day. For an overview, you might want to look at the poster PDF in this repository. 

While the first, basic end-to-end workflow was complete Feb 7 2016, and tagged v 0.1, I found in April-May that I needed to rewrite the modeling data prep to accommodate the poor quality of very many comp stars in VPhot sequences (as evaluated both by large regression residuals and by widely ranging mag errors listed in the matching VSP charts). 

Current working version is **V 1.1.4** (May 2017), and it works. Bug and typo fixes aside, this code is frozen. 

The V 1.1.4 workflow is:

**Input.R (pre-calibration)_______________________________________________**

 1. **renameObject()** --- Rarely used. Renames FITS files *and* FITS headers' "OBJECT" value from old to new Object name, across all FITS files (and subdirectories) within AN top folder. For example, renames all "Landolt_001" to "Std_001". This was needed whenever object name and file name in images (probably due to inclusion in ACP nightly observing plan) had been poorly chosen. renameObject() must be run first if it is used at all.
 2. **precheck()** --- Run on every data set without exception, iteratively until all errors fixed. Checks that a field-of-view FOV file exists for every FITS file. A FOV file is required for every FITS file included in the model and AAVSO reporting. FOV files are obviously not required for FITS files acquired for other uses, notably those acquired as "Burn" images to prepare new FOV files for use in future photometric runs.
 1. **beforeCal()** --- Calls all 3 functions in correct order: copyToUr(), renameACP(), and prepareForCal(), which takes renames and arranges FITS files (light, flat, bias, and darks) for the next step, which is calibration (I use MaxIm DL). Its steps in order are: [copyToUr()]: Copies (backs up) all (possibly renamed) target FITS files to a new /Ur directory, for safe keeping. [renameACP()]: Renames all _target_ FITS files (not calibration files) from ACP-format names (like TargetName-S001-R001-C003-V-dupe1.fts) to my own strictly sequential file naming system(like TargetName-0003-V.fts). Each new sequence number is based on time (JD) from the FITS header, and it is validated and sorted by time across all files of the same TargetName. Renaming is not fooled by dupe or duplicate naming across subdirectories. TargetNames and Filter IDs are verified for equality between the original file name and the FITS header. All target FITS files are then collected in the top AN folder. A table of old vs new file names is written to a text file, and the data frame is stored as a .RData file for immediate reload into R. [prepareForCal()]: Sets up folders /Calibration (all flat, dark, bias images for the AN), /CalibrationMasters (calibration images duplicated here & ultimate home of MaxIm generated calibration master frames), /Ur (backup of all renamed target FITS files; may move this forward to renameACP() (or similar fn) so that the Ur files are really Ur), and /Photometry (catchall for metadata about this AN folder).

**Calibration of all target images in MaxIm DL (not in R)______________________**

  1. Ensure all needed raw flat, dark, and bias frames are collected in /CalibrationMasters.
  1. In MaxIm, "Set Calibration" to this /CalibrationMasters folder, and "Replace w/Masters".
  1. Set up "Batch Convert and Save" function (under Edit): Get names of all FITS from /Uncalibrated folder, check the Calibration box, set Destination to /Calibrated folder, click OK. Now your calibrated FITS files (unfortunately with MaxIm's alternate extension .fit) are saved and ready. Go back to R.

**Input.R (post-calibration)______________________________________________**

 1. **finishFITS()**  ---  Resets calibrated FITS file extensions to .fts. Deletes all non-master files from /CalibrationMasters, verifies via each FITS header that all calibrations happened OK, greatly cleans up directory structure.
 1. **make_df_master()**  ---  Prepares inputs for and invokes APT software to build R data frame suitable for later use in photometric model building and AAVSO. APT is a bit limited in its aperture math, and I may later write my own aperture math in R, but for now APT is very reliable in its star marking in images. This function draws heavily on previously prepared FOV (field-of-view), one text file per field of view, which roughly corresponds to a VPhot "sequence" (badly named by AAVSO) of all the relevant target, check, and standard comp stars; my FOV file includes some additional metadata that might as well be bundled with the VPhot data, right in the same files.
**df_master** is the first immutable data frame; it is saved in the AN folder's "Photometry" subfolder. Once set up it should not be changed. This is not hard, as it is automated and hard to make a mistake in the steps constructing it. If it is changed, *all* steps below (Model.R & Predict.R) simply have to be rerun.

**Model.R_____________________________________________________________**

 1. **modelOneFilter()**  ---  Run once for each filter with observations of both comp stars and target stars. Preforms a mixed-model regression for one filter, and writes an empty template omit.txt file. The model uses all eligible comp stars, especially those with minimum signal-to-noise ratio and that are not saturated *in the original FITS image*. This function returns a one-filter model list for the entire Astronight. This one-filter model list will be compounded into a masterModelList once all filters have been modeled. The masterModelList includes all data that later will be needed to predict check-star and target-star magnitudes, plus all the other data required for an AAVSO report. 
 2. **omitSerial()** --- Adds #OBS lines to omit.txt so that observations with the specified Serial numbers are no longer used in model building. Model building is frequently an iterative process--model, plot, omit, and repeat--done independently in each filter.

**Plots.R ___________________________________________________________**

1. **modelPlots()** -- Writes a large number of ggplot2 diagnostic plots of model into RStudio's plot window. I use these to identify anything I want to omit from the model input--single observations (one comp star in one image, by Serial number), individual images, or all instances of a comp star either in one filter or in all filters. The user can specify these omissions by editing omit.txt directly, or sometimes more conveniently when omitting individual observations, by running omitSerial().
1. **[Later, we will also need some diagnostic plots for selecting/combining Target results for building
our AAVSO reports.]**

**Predict.R ______________________________________________________________**

 1. **predictAll()**  ---  Takes the final masterModelList and (1) makes predictions of target and check-star. Applies color-index / transform corrections. Reads user control file pre-predict.txt to ensure that a single set of comp stars is applied to all points in a one-night ("Stare") light curve, particularly for eclipser light curves. Does a lot of bookkeeping. This is the final part of the true photometric computation chain--later steps merely aid user curation and report generation (essentially: bookkeeping).
 2. **markupReport()** -- Prints a table of target observations, along with: magnitude errors of various types, Julian date (time), check-star magnitudes, and a lot more. Aids in the user's last chance to curate target data before final report generation. Most changes coming out of this report are of two types: Omitting unsuitable target observations, and combining target observations, which requires careful handling of observation time and reported magnitude errors. The user uses this report to edit omit.txt, which is input to the next and last processing step.
 1. **writeAAVSO()**  ---  Takes user control file 'omit.txt' and target magnitude etc data from predictAll(), and generates a format-perfect AAVSO report, ready for upload. (In some cases, the user will want to annotate the report, most often to make comments in the report header about consequences of the data, e.g., recapture of a lost lightcurve phase, or detected outburst of an observed star.

**[end of README] ______________________________________________________________**


