# adi2_rsfmri

## Contributors: Hannah Sophie Heinrichs, Frauke Beyer

## About

Project 2 of the ADIPOSITAS study is about potential changes in resting state functional connectivity induced by bariatric surgery. This repository contains all scripts necessary for data analysis (descriptive and inference statistics) and the manuscript draft.

`abs_path` - absolute path in server
`.Rproj` - R project; open Rmd file in this project to ensure relative paths are correct

### code
This folder containts all code for the MATLAB analysis with the SwE toolbox. Explanations of code below.

`run_or_display.m`: analysis script; conducts analysis as specified in the beginning of the script - for non-parametric estimation it will run or display all contrasts in a loop. To use script read [How to run analysis in MATLAB](#how-to-run-the-analysis-in-matlab) below.

`RunModelGroupTime.m`, `RunModelBMI.m` and `RunModelFD.m` as well as `AllTPEval.m` and `SingleTPEval.m` are functions called by `run_or_display.m`


### manuscript
Folder with files to create a manuscript. CAVE: introduction and rest are in separate files.

`draft_rsfmri` - document for the manuscript comprising method, result, and discussion section.

`intro_rsfmri` - introduction to manuscript

`get_subject_ID_group_tp_forSwE.R` - R file with function `txt_for_swe(group, tp)` that creates final data frame

create_df_and_txt_for_swe.R

### report
Folder contains figures and tables to include into the manuscript that were previously calculate or set up.

## How to run the analysis in MATLAB

(1) General

1. Specify the analysis you want to perform in the run_or_display.m file by setting the respective parameters
2. Run script
3. If estimations were previously performed, set `ONLY_DISPLAY = true`
4. For display only and parametric estimation, follow the prompts by SPM

  - Prompts will differ depending on whether results were also displayed earlier.
	   - If so, contrasts were already specified. Only select the desired contrast or - in case contrasts are incomplete or inappropriate - re-run the analysis by running the script with `ONLY_DISPLAY = false` and `OVERWRITE = true` and Step 1-3. This will delete previouly specified contrasts.
	   - If not previously specified, please specify according to the information that was displayed in the console before starting SPM

(2a) Parametric estimation
1. Select the desired contrast.
2. Apply masking: [none]
3. Title for comparison: accept given title
4. Inference type: [voxelwise]
5. P value adjustment to control: [FDR]
6. P value (FDR): [0.05]

(2b) Non-parametric estimation
1. All contrasts will be evaluated successively
2. Contrast Type: [Activation] / [Deactivation]
3. Title for comparison: accept given title
4. Inference type: [clusterwise]
5. P value adjustment to control: [FWE]
6. P value (FWE): [0.05]

## How to create the script in RMarkdown

* open `draft_rsfmri.Rproj` (ensures correct relative paths)
* open `draft_rsfmri.Rmd` within the project
* check if packages are installed and package versions are consistent
* Knit markdown document to PDF
