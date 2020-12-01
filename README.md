# adi2_rsfmri

## Contributors: Hannah Sophie Heinrichs, Frauke Beyer

## About

Project 2 of the ADI study is about potential changes in resting state functional connectivity induced by bariatric surgery. This repository contains all scripts necessary for data analysis (descriptive and inference statistics) and the draft.

`run_or_display.m`: analysis script; conducts analysis as specified in the beginning of the script - for non-parametric estimation it will run or display all contrasts in a loop. To use script read > How to run analysis < below.

`RunModelGroupTime.m`, `RunModelBMI.m` and `RunModelFD.m` as well as 'AllTPEval.m' and 'SingleTPEval' are functions called by `run_or_display.m`

`draft_rsfmri` - Rproject, Rmarkdown and PDF document for the draft.

`get_subject_ID_group_tp_forSwE.R` - R file for descriptive analysis


## How to run the analysis

(1) General

1. Specify the analysis you want to perform in the run_or_display.m file by setting the respective parameters
2. Run script
3. If estimations were previously performed, set ONLY_DISPLAY = true
4. For display only and parametric estimation, follow the prompts by SPM

  - Prompts will differ depending on whether results were also displayed earlier.
	   - If so, contrasts were already specified. Only select the desired contrast or - in case contrasts are incomplete or inappropriate - re-run the analysis by running the script with "ONLY_DISPLAY = false" and "OVERWRITE = true" and Step 1-3. This will delete previouly specified contrasts.
	   - If not previously specified, please specify according to the information that was displayed in the console before starting SPM

(2a) Parametric estimation
1. Select the desired contrast.
2. Apply masking: >> none <<
3. Title for comparison: accept given title
4. Inference type: >> voxelwise <<
5. P value adjustment to control: >> FDR <<
6. P value (FDR): >> 0.05 <<

(2b) Non-parametric estimation
1. All contrasts will be evaluated successively
2. Contrast Type: >> Activation <<
3. Title for comparison: accept given title
4. Inference type: >> voxelwise <<
5. P value adjustment to control: >> FDR <<
6. P value (FDR): >> 0.05 <<
