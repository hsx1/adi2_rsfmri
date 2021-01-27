# adi2_rsfmri

## Contributors: Hannah Sophie Heinrichs, Frauke Beyer

## About

Project 2 of the ADIPOSITAS study is about potential changes in resting state functional connectivity induced by bariatric surgery. This repository contains all scripts necessary for data analysis (descriptive and inference statistics) and the manuscript draft.

`abs_path` - absolute path in server

`.Rproj` - R project; open Rmd file in this project to ensure relative paths are correct

### code
This folder contains all code for the MATLAB analysis with the SwE toolbox. Explanations of code below.

`run_or_display.m`: analysis script; conducts analysis as specified in the beginning of the script - for non-parametric estimation it will run or display all contrasts in a loop. To use script read [How to run analysis in MATLAB](#how-to-run-the-analysis-in-matlab) below.

`RunModelGroupTime.m`, `RunModelBMI.m` and `RunModelFD.m` as well as `AllTPEval.m` and `SingleTPEval.m` are functions called by `run_or_display.m`

`spm2csv.m` - saves results depending on activation and deactivation depicted by SPM/SwE into csv table in output folder as `results_act.csv` or `results_deact.csv`, and `info.csv` containing estimation specifics.

### report
This folder contains scripts called upon by `draft_rsfmri.R` in Order to construct data frames or construct figures and tables.

`create_sample_df.R` - creates data frame of sample to analze.

`tables_FC.R` - constructs result tables on all estimated models with statistics, p-values and coordinates, that were retrieved by
 `run_or_display.m` calling `spm2csv.m`, and anatomical labels that have been previously been extract with **SPM Anatomy Toolbox 2.2c** (saving txt files).

`FDFCcorrelations.csv` - cached data from (time-consuming) analysis of `...Preprocessing/qa/rs_qa/group_level_QA/QC_FC_correlations/calc_qc_fc_correlations.R`

### manuscript
Folder with files to create a manuscript. CAVE: introduction and rest are in separate files.

`draft_rsfmri` - document for the manuscript comprising method, result, and discussion section.

`intro_rsfmri` - introduction to manuscript

`df2SwEtxt.R` - R file with function `get_txt_for_swe(group, tp)` that creates final data frame

## Analysis pipeline

(1) [Run analysis](#how-to-run-the-analysis-in-matlab) with matlab.

(2) [Extract anatomical labels](#how-to-spm-anatomy-toolbox) with SPM Anatomy toolbox

(3) [Create manuscript](#how-to-create-the-script-in-rmarkdown)

 - run `...Preprocessing/qa/rs_qa/group_level_QA/QC_FC_correlations/calc_qc_fc_correlations.R` from within `draft_rsfmri.Rmd` and make sure report folder contains `FDFCcorrelations.csv`

 - (make sure `tables_FC.R` contains correct path specifications)

 - knit `draft_rsfmri.Rmd` to pdf

## Detailed How To's for analysis pipeline

### How to run the analysis in MATLAB

(1) General

1. Specify the analysis you want to perform in the `run_or_display.m` file by setting the respective parameters (therefore view description in script)
2. Run script
3. If estimations were previously performed, set `ONLY_DISPLAY = true`
4. For display only and parametric estimation, follow the prompts by SPM

  - Prompts will differ depending on whether results were also displayed earlier.
	   - If so, contrasts were already specified. Only select the desired contrast or - in case contrasts are incomplete or inappropriate - re-run the analysis by running the script with `ONLY_DISPLAY = false` and `OVERWRITE = true` and Step 1-3. This will delete previouly specified contrasts.
	   - If not previously specified, please specify according to the information that was displayed in the console before starting SPM

(2) Non-parametric estimation
1. All contrasts will be evaluated successively
2. Contrast Type: <kbd><samp>Activation</samp></kbd> / <kbd><samp>Deactivation</samp></kbd>
3. Title for comparison: accept given title
4. Inference type: <kbd><samp>clusterwise</samp></kbd>
5. P value adjustment to control: <kbd><samp>FWE</samp></kbd>
6. P value (FWE): `0.05`

- saves `Swe.mat`, nii-files and csv file with output in output folder.

### How to SPM Anatomy toolbox

To extract and process significant results.

1. open the **SPM Anatomy toolbox 2.2c** by entering "Anatomy" into Matlab Console (if version doen't match, echeck spm12/toolbox folder)
2. import Image from result folder, `swe_vox_zTstat_c02.nii` for activation, `swe_vox_zTstat_c02.nii` for deactivation.
3. keep premultiply `1`, enter thresholds, e.g. Z = 3.09 and k = 55 (retrieve exact thresholds from `info.csv` or spm output image)[alternatively to step 2+3, you can save the thresholded image from SPM and don't need to give thesholds in the Anatomy toolbox then]
4. save tables in txt table by pressing <kbd><samp>Tab</samp></kbd>
5. in txt file naming "result_resport.txt" in folder "Results", save "model_name" "path_to_txt" and "cluster_threshold"


### How to create the script in RMarkdown

* open `draft_rsfmri.Rproj` (ensures correct relative paths)
* open `draft_rsfmri.Rmd` within the project
* check if packages are installed and package versions are consistent
* Knit markdown document to PDF
