# adi2_rsfmri

Necessary files for analysis of Project 2 with the ADI cohort. This repository contains all scripts necessary for data analysis (descriptive and inference statistics) and the draft.

`run_or_display.m`: analysis script; conducts analysis as specified in the beginning of the script - for non-parametric estimation it will run or display all contrasts in a loop

`RunModelGroupTime.m`, `RunModelBMI.m` and `RunModelFD.m` are functions called
by `run_or_display.m`

`draft_rsfmri` - Rproject, Rmarkdown and PDF document for the draft.

`get_subject_ID_group_tp_forSwE.R` - R file for descriptive analysis
