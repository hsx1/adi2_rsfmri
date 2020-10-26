%% Define parameters
% _Input_: Relevant directories (INFO_DIR, OUT_DIR, MASK_DIR), model 
% definition (MODEL, COVARIATES), target region (ROI), estimation and 
% inference specifics (MASK, WILD_BOOT, INFERENCE) and action aguments 
% (DISPLAY_ONLY, OVERWRITE)
% _Output_: saves SwE.mat and contrast imagies (*.nii)
% 
% Use string for stating *ROI and preprocessing step* and z-transform in
% in the following format, e.g. ROI_PREP = 'PCC_min_z' or 'NAcc_aroma'.
% ROI_PREP must be in a cell array of multiple or a single cell!
% For roi enter either 'PCC' or 'Nacc'. For prep enter either 'min', 
% 'aroma', 'cc', or 'gsr'. If desired, append 'z'.
% To evaluate all specific ROI/ preprocessing combinations type
% ROI_PREP = readcell(fullfile(INFO_DIR,'ROIs.txt'), 'Delimiter',' ','Whitespace',"'");
%
% Define the *MODEL* you want to test
% 'grouptime': FC ~ group + time + group*time
% 'bmi': FC ~ avgBMI + BMIcgn
% 'bmiIG': FC ~ avgBMI + BMIcgn - only for the intercention group
% 'bmi2tp': FC ~ avgBMI + BMIcgn - only for BL and FU1
% 'fd': FC ~ (avgBMI + BMIcgn +) avgFD + FDcgn + age + sex
% 'fdIG': FC ~ (avgBMI + BMIcgn) + avgFD + FDcgn + age + sex - only for the intercention group
%
% Specify the model further with an integer for *COVARIATES definition*
% 11: group-timef-age-sex-meanFD
% 12: group-timef-age-sex
% 13: group-timec-age-sex-meanFD
% 14: group-timec-age
% 21: bmi-age-sex-meanFD
% 22: bmi-age-sex
% 31: bmi-fd-age-sex
% 32: fd-age-sex
%
% Use string to define *output directory*, e.g. OUT_DIR = 
% '/data/hu_heinrichs/Documents/Swe_results/';
% The SwE.mat and output files will be saved under the OUT_DIR in a folder 
% named after roi_prep with a subfolder for the model will be created.
%
% INFO_DIR is the *directory containing all txt files* needed for the 
% regressors that will be entered in model according to the model 
% specification, e.g. '/data/pt_02161/Analysis/Project2_resting_state/seed-based/Second_level /SwE_files/'
%
% WILD_BOOT if set on true will use non-parametri Wild Bootstrap instead 
% of parametric estimation. It will automatically run all relevant contrast 
% of the model (as defined in the script).
% NOTE: the contrast implemented in the script are optimized for 
% 'bmi-age-sex-meanFD'.
% 
% When using WILD_BOOT = true, you must specify which kind of inference
% is used for bootstrapping:
% 'voxel': voxelwise 
% 'cluster': clusterwise
% 'tfce': TFCE
% 
% *** IMPORTANT ***
% If DISPLAY_ONLY is set on false and the model has not yet been estimated,
% it will be. In case of parametric estimation if will estimate beta
% coefficients; in case of non-parametric bootstrapping, it will estimate
% contrast coefficients.
% If DISPLAY_ONLY is set on false and the model has been estimated, the
% results will only then change, if OVERWRITE is set on true.
% define relevant input and output directories
%% ========================================================================

% create a struct, with all important parameters
ABS_DIR = readcell("abs_path.txt");
ABS_DIR = ABS_DIR{1};
param.OUT_DIR = fullfile(ABS_DIR,'/Results/Project2_resting_state/connectivity/Analysis/'); %'/data/pt_02161/Results/Project2_resting_state/connectivity/Analysis/preliminary_analysis/'; 
param.INFO_DIR = fullfile(ABS_DIR,'/Analysis/Project2_resting_state/seed-based/Second_level /SwE_files/');
param.MASK_DIR = fullfile(ABS_DIR, '/Analysis/Project2_resting_state/seed-based/Brain_masks/');
% define ROI
roi_prep = readcell(fullfile(param.INFO_DIR,'ROIs.txt'), 'Delimiter',' ','Whitespace',"'");
param.ROI_PREP = {roi_prep{[4]}}; % {ROI_PREP{[4, 6, 12, 14]}} or {'Nacc_cc_z','Nacc_gsr_z','PCC_cc_z','PCC_gsr_cc'}

% Model definition
% All three models have unique options for covariate definition, the
% association to a model is indicated by the tens digit (GroupTime_: 1_; 
% BMI_ = 2_; FD_ = 3_) the specific covariate combination by the ones digit

param.MODEL = {'grouptime'}; % {'grouptime','grouptime2tp'} % {'bmi','bmiIG','bmi2tp'} % {'fd','fdIG'}
param.COVARIATES = [11]; % [11, 12];                    % [21, 22];                % [31, 32]; 

% define masking and type of inference
param.MASK = 'brain';                  % 'brain' or 'gm'
param.WILD_BOOT = false;             % false
param.INFERENCE_TYPE = {'voxel'};   % {'voxel','cluster','tfce'};
% analysis parameter (estimate or display?)
param.ONLY_DISPLAY = true;         % false
param.OVERWRITE = false;            % false
param.VIEWSEC = 1; % for ONLY_DISPLAY: seconds you want to view the results

% set path for spm and path with my functions
addpath(genpath('/data/pt_life/data_fbeyer/spm-fbeyer'))
addpath('/data/pt_02161/Analysis/Project2_resting_state/seed-based/Second_level /adi_analysis')

% use function to run or display models
if param.COVARIATES(1) < 20
    RunModelGroupTime(param)
elseif param.COVARIATES(1) < 30
    RunModelBMI(param)
elseif param.COVARIATES(1) < 40
    RunModelFD(param)
end
