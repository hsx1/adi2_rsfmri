%% Define parameters
% _Input_: Relevant directories (INFO_DIR, OUT_DIR, MASK_DIR), model 
% definition (MODEL, COVARIATES), target region (ROI), estimation and 
% inference specifics (MASK, WILD_BOOT, INFERENCE) and action aguments 
% (DISPLAY_ONLY, OVERWRITE); MASKS are further specified 
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
% 'alltp'
% 'singletp'
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
% 41: network-age-sex
% 42: network-meanFD-age-sex
%
% Specify directories relativ to project folder. The absolute path to the
% project folder should be specified in a file called "abs_path.csv".
%
% Use string to define *output directory*, e.g. OUT_DIR = 
% '/Documents/Swe_results/';
% The SwE.mat and output files will be saved under the OUT_DIR in a folder 
% named after roi_prep with a subfolder for the model will be created.
%
% INFO_DIR is the *directory containing all txt files* needed for the 
% regressors that will be entered in model according to the model 
% specification, e.g. '/Analysis/Project2_resting_state/seed-based/Second_level /SwE_files/'
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
% Especially for results of non-parametric estimation, you need to specify
% how long you want to look at the results in VIEWSEC. Alternatively you
% can also set breakpoints in the funtion scripts.

%% ========================================================================
% set path for spm and path with my functions
% swe version 2.2.1 download of development 
addpath('/data/u_heinrichs_software/MATLAB/spm12/')
addpath('/data/u_heinrichs_software/MATLAB/spm12/toolbox/SwE-toolbox-2.2.1-1292020')
% swe version 2.1.1
% addpath(genpath('/data/pt_life/data_fbeyer/spm-fbeyer'))


% import absolute path for project
% SET PATH TO code_and_manuscript FOLDER AS CURRENT DIRECTORY !!
cd("/data/pt_02161/Analysis/Project2_resting_state/seed-based/Second_level /code_and_manuscript/code")
ABS_DIR = readcell("../abs_path.csv");
ABS_DIR = ABS_DIR{1};
addpath(fullfile(ABS_DIR,'/Analysis/Project2_resting_state/seed-based/Second_level /code_and_manuscript')) % for functions
% create a struct, with all important parameters
param.OUT_DIR = fullfile(ABS_DIR,'/Results/Project2_resting_state/connectivity/Analysis/'); %'/data/pt_02161/Results/Project2_resting_state/connectivity/Analysis/preliminary_analysis/'; 
param.INFO_DIR = fullfile(ABS_DIR,'/Analysis/Project2_resting_state/seed-based/Second_level /SwE_files/');
param.MASK_DIR = fullfile(ABS_DIR, '/Analysis/Project2_resting_state/seed-based/Brain_masks/');
% file names of masks
param.MASK_GM = 'mni_icbm152_gm_tal_nlin_sym_09a_resampl_bin.nii,1';
param.MASK_B = 'MNI_resampled_brain_mask.nii,1';

% define ROI
roi_prep = readcell(fullfile(param.INFO_DIR,'ROIs.txt'), 'Delimiter',' ','Whitespace',"'");
param.ROI_PREP = {roi_prep{[4, 6, 12, 14]}}; % {roi_prep{[4, 6, 12, 14]}} or {'Nacc_cc_z','Nacc_gsr_z','PCC_cc_z','PCC_gsr_z'}

% Model definition
% All three models have unique options for covariate definition, the
% association to a model is indicated by the tens digit (GroupTime_: 1_; 
% BMI_ = 2_; FD_ = 3_) the specific covariate combination by the ones digit
param.MODEL = {'grouptime'}; % {'grouptime','grouptime2tp'} % {'bmi','bmiIG','bmi2tp'} % {'fd','fdIG'} % {'alltp'} % {'singletp} 
param.COVARIATES = [11];     % [11, 12];                    % [21, 22];                % [31, 32];  % [41, 42]    % []

% define masking and type of inference
param.MASK = 'brain';               % 'brain' or 'gm'
param.EXCLFD = false;                % false
param.WILD_BOOT = false;             % false
param.INFERENCE_TYPE = {'cluster'};   % {'voxel','cluster','tfce'};
% analysis parameter (estimate or display?)
param.ONLY_DISPLAY = false;         % false
param.OVERWRITE = true;            % false
param.VIEWSEC = 1; % for ONLY_DISPLAY: seconds you want to view the results
% param.ACTION = 'estimate','display','overwrite'

% use function to run or display models
if strcmp(param.MODEL,'singletp')
    SingleTPEvalAge(param)
elseif strcmp(param.MODEL,'alltp')
    AllTPEval(param)
else
    % for all other models
    if param.COVARIATES(1) < 20
        RunModelGroupTime(param);
    elseif param.COVARIATES(1) < 30
        RunModelBMI(param);
    elseif param.COVARIATES(1) < 40
        RunModelFD(param);
    end
end
