%% ========================================================================
%% Define parameters
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
%
% Specify the model further with an integer for *covariate definition*
% 11: Model1_group_timef_cage_sex_meanFD
% 12: Model2_group_timef_cage_sex
% 13: Model3_group_timec_cage_sex_meanFD
% 14: Model4_group_timec_cage
% 21: Model1_bmi_cage_sex_meanFD
% 22: Model2_bmi_cage_sex
%
% Use string to define *output directory*, e.g. OUT_DIR = 
% '/data/hu_heinrichs/Documents/Swe_results/';
% The SwE.mat and output files will be saved under the OUT_DIR in a folder 
% named after one_roi_prep with a subfolder for the model will be created.
%
% INFO_DIR is the *directory containing all txt files* needed for the 
% regressors that will be entered in model according to the model 
% specification, e.g. '/data/pt_02161/Analysis/Project2_resting_state/seed-based/Second_level /SwE_files/'
%
% WILD_BOOT if set on true will use non-parametri Wild Bootstrap instead 
% of parametric estimation. It will automatically run all relevant contrast 
% of the model (as defined in the script).
% NOTE: the contrast implemented in the script are optimized for 
% 'Model1_group_timef_cage_sex_meanFD' and 'Model1_bmi_cage_sex_meanFD'.
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

%% ========================================================================
% INFORMATION: bootstrapping is changed to 999 (time consuming)


%% ========================================================================
% please define parameters (default as comment (%) behind each variable)
INFO_DIR = '/data/pt_02161/Analysis/Project2_resting_state/seed-based/Second_level /SwE_files/';
OUT_DIR = '/data/pt_02161/Results/Project2_resting_state/connectivity/Analysis/'; %'/data/pt_02161/Results/Project2_resting_state/connectivity/Analysis/preliminary_analysis/'; 
MODEL = {'grouptime'}; % {'grouptime','bmi'};
COVARIATES = [11];  % [11,12,21,22]; 
ROI_PREP = readcell(fullfile(INFO_DIR,'ROIs.txt'), 'Delimiter',' ','Whitespace',"'");
ROI_PREP = {ROI_PREP{[4]}}; % {ROI_PREP{[4, 6, 12, 14]}} or  {'Nacc_cc_z','Nacc_gsr_z','PCC_cc_z','PCC_gsr_cc'}

WILD_BOOT = true; %true
INFERENCE_TYPE = {'cluster','tfce'}; %{'voxel','cluster','tfce'};

ONLY_DISPLAY = false; %false
OVERWRITE = true; %false


%% ========================================================================
addpath(genpath('/data/pt_life/data_fbeyer/spm-fbeyer'))
%addpath(genpath('/data/hu_heinrichs/Documents/MATLAB/spm12'))

% loop over all models, only all covariates, all inference types for
% wild_boot = true
for i = 1:length(MODEL)
    if strcmp(MODEL{i},'bmi')
        covariates = COVARIATES(COVARIATES > 15);
    elseif strcmp(MODEL{i},'grouptime')
        covariates = COVARIATES(COVARIATES < 15);
    end
    for j = 1:length(covariates)
        for k = 1:length(INFERENCE_TYPE)
            RunSwe(INFO_DIR, OUT_DIR, MODEL{i}, covariates(j), ROI_PREP, WILD_BOOT, INFERENCE_TYPE{k}, ONLY_DISPLAY, OVERWRITE);
        end
    end
end

%% ========================================================================
function display_message(COVARIATES)
% displays message according to the COVARIATES, that define the model
if COVARIATES == 11
    fprintf('%s\n',...
    'Model1_group_timef_cage_sex_meanFD: In case of parametric estimation',...
    'please type in the following contrasts for all runs',...
    'group [1 1 1 -1 -1 -1 0 0 0]',...
    'time [-1 0 1 -1 0 1 0 0 0; 0 -1 1 0 -1 1 0 0 0]',...
    'group*time [-1 0 1 1 0 -1 0 0 0; 0 -1 1 0 1 -1 0 0 0]')
elseif COVARIATES == 12
    fprintf('%s\n',...
    'Model2_group_timef_cage_sex: In case of parametric estimation',...
    'please type in the following contrasts for all runs',...
    'group [1 1 1 -1 -1 -1 0 0]',...
    'time [-1 0 1 -1 0 1 0 0; 0 -1 1 0 -1 1 0 0]',...
    'group*time [-1 0 1 1 0 -1 0 0; 0 -1 1 0 1 -1 0 0]')
elseif COVARIATES == 13
    fprintf('%s\n',...
    'Model3_group_timec_cage_sex_meanFD: In case of parametric estimation',...
    'please type in the following contrasts for all runs',...
    'group [0 1 -1 0 0 0 0 0]',...
    'time [0 0 0 0.5 0.5 0 0 0]',...
    'group*time [0 0 0 1 -1 0 0 0]')
elseif COVARIATES == 14
    fprintf('%s\n',...
    'Model4_group_timec_cage_sex: In case of parametric estimation',...
    'please type in the following contrasts for all runs',...
    'group [0 1 -1 0 0 0 0]',...
    'time [0 0 0 0.5 0.5 0 0]',...
    'group*time [0 0 0 1 -1 0 0]')
elseif COVARIATES == 21
    fprintf('%s\n',...
    'Model1_BMI_meanFD_cage_sex: In case of parametric estimation',...
    'please type in the following contrasts for all runs',...
    'avgBMI [0 1 0 0 0 0]',...
    'BMIcgn [0 0 1 0 0 0]')
elseif COVARIATES == 22
    fprintf('%s\n',...
    'Model2_BMI_cage_sex: In case of parametric estimation',...
    'please type in the following contrasts for all runs',...
    'avgBMI [0 1 0 0 0]',...
    'BMIcgn [0 0 1 0 0]')
end

end

%% ------------------------------------------------------------------------
function RunSwe(INFO_DIR, OUT_DIR, MODEL, COVARIATES, ROI_PREP, WILD_BOOT, INFERENCE_TYPE, ONLY_DISPLAY, OVERWRITE)
% Carries out analysis according to specified parameters. 
% RunSwe makes use of three functions: SpecifyModel, RunModel and
% Display Results. If SpecifyModel and RunModel has run once and the
% SwE.mat file is 

if ONLY_DISPLAY && not(WILD_BOOT)
    display_message(COVARIATES)
    % wait so displayed information can be read
    pause(7)
end

nrun = length(ROI_PREP);
for crun = 1:nrun
    % estimate any of the six models for every one_roi_prep in ROI_PREP
    if not(WILD_BOOT)
        % parametric estimation
        wild_con = false;
        % check if analysis has been conducted already
        [out_folder, exist_already] = create_out_folder(OUT_DIR, MODEL, ROI_PREP{crun}, COVARIATES, wild_con, INFERENCE_TYPE);
        clear matlabbatch
        [matlabbatch,location_SwE_mat] = SpecifyModel(ROI_PREP{crun}, MODEL, COVARIATES, out_folder, INFO_DIR, wild_con, INFERENCE_TYPE); 
        if ONLY_DISPLAY
            spm('defaults', 'FMRI'); 
            matlabbatch = DisplayResults(location_SwE_mat);
            spm_jobman('run', matlabbatch);
        elseif not(exist_already) || OVERWRITE
            spm('defaults', 'FMRI'); 
            matlabbatch = RunModel(matlabbatch, location_SwE_mat);
            spm_jobman('run', matlabbatch);
        end
        % manually type in contrasts (issue: see below)
    else
        % non-parametric estimation with Wild-Bootstrap
        for wild_con = 1:3
            % for 'bmi' model, there is only two contrasts
            if not(strcmp(MODEL,'bmi') && wild_con == 3)
                % check if analysis has run already
                [out_folder, exist_already] = create_out_folder(OUT_DIR, MODEL, ROI_PREP{crun}, COVARIATES, wild_con, INFERENCE_TYPE);
                clear matlabbatch
                [matlabbatch,location_SwE_mat] = SpecifyModel(ROI_PREP{crun}, MODEL, COVARIATES, out_folder, INFO_DIR, wild_con, INFERENCE_TYPE); 
                if ONLY_DISPLAY
                    spm('defaults', 'FMRI');
                    matlabbatch = DisplayResults(location_SwE_mat);
                    spm_jobman('run', matlabbatch);
                elseif not(exist_already) || OVERWRITE
                    spm('defaults', 'FMRI');
                    matlabbatch = RunModel(matlabbatch, location_SwE_mat);
                    spm_jobman('run', matlabbatch);
                end
            end
        end
    end
 
% save results of contrasts in folder one_roi_prep and subfolder modelx
% currently not possible to automate contrasts and save results
% https://github.com/NISOx-BDI/SwE-toolbox/issues/135
end
spm('Quit')

end

%% ------------------------------------------------------------------------
function [matlabbatch,location_SwE_mat] = SpecifyModel(one_roi_prep, MODEL, COVARIATES, out_folder, INFO_DIR, wild_con, INFERENCE_TYPE)
% Specifies design matrix and estimation procedure as preparation for the
% estimation to be carried out (by function RunModel); the specified model 
% is stored in SwE.mat file in location_SwE_mat.

%% Directory to get information for model definition
cd(INFO_DIR)

%% Directory to save SwE file
% save directory of (new) folder
matlabbatch{1}.spm.tools.swe.smodel.dir = {out_folder};

%% Load Scans -------------------------------------------------------------
% constructing cell array with the scans (FC maps per subject and time point)
scans_dir = readcell(fullfile(INFO_DIR,'scans.txt'), 'Delimiter',' ','Whitespace',"'");
if contains(one_roi_prep, 'gsr')
    % special case if preprocessing is 'gsr'
    if endsWith(one_roi_prep,'_z')
        FC_files = strcat(lower(one_roi_prep(1:end-5)),'cc_gsr_seed_correlation_z_trans.nii');
    else
        FC_files = strcat(lower(one_roi_prep(1:end-3)),'cc_gsr_seed_correlation_trans.nii');
    end
else
    if endsWith(one_roi_prep,'_z')
        FC_files = strcat(lower(one_roi_prep(1:end-2)),'_seed_correlation_z_trans.nii');
    else
        FC_files = strcat(lower(one_roi_prep),'_seed_correlation_trans.nii');
    end
end
scans_of_roi = fullfile(scans_dir, one_roi_prep, FC_files);

% load scans to matlabbatch
matlabbatch{1}.spm.tools.swe.smodel.scans = scans_of_roi;
                                         
%% SwE type ---------------------------------------------------------------
% .Modified
% .. Define Groups
matlabbatch{1}.spm.tools.swe.smodel.type.modified.groups = readmatrix('group.txt');

% ..Visits
matlabbatch{1}.spm.tools.swe.smodel.type.modified.visits = readmatrix('tp.txt');

% small sample adjustment (4 = type C2)
matlabbatch{1}.spm.tools.swe.smodel.type.modified.ss = 4;
% degrees of freedom type (3 = approx III)
matlabbatch{1}.spm.tools.swe.smodel.type.modified.dof_mo = 3;

%% Subjects ---------------------------------------------------------------
matlabbatch{1}.spm.tools.swe.smodel.subjects = readmatrix('subjNr.txt');

%% Covariates (Design matrix) ---------------------------------------------
% desgin matrix for model estimation
cov = create_design_matrix(MODEL, COVARIATES);
matlabbatch{1}.spm.tools.swe.smodel.cov = cov;
  
%% Masking ----------------------------------------------------------------
%% Multiple Covariates (none)
matlabbatch{1}.spm.tools.swe.smodel.multi_cov = struct('files', {});

%% Masking (none)
% .Threshold masking
matlabbatch{1}.spm.tools.swe.smodel.masking.tm.tm_none = 1;
% ..Implicit Mask (yes)
matlabbatch{1}.spm.tools.swe.smodel.masking.im = 1;
% .. Explicit Mask
matlabbatch{1}.spm.tools.swe.smodel.masking.em = {'/data/pt_02161/Analysis/Project2_resting_state/seed-based/Second_level /MNI_resampled_brain_mask.nii,1'};

%% Non-parametric Wild Bootstrap ------------------------------------------
% . No
if ~wild_con
    matlabbatch{1}.spm.tools.swe.smodel.WB.WB_no = 0;
else
    % . Yes
    % .. Small sample adjustments for WB resampling (4 = type C2)
    matlabbatch{1}.spm.tools.swe.smodel.WB.WB_yes.WB_ss = 4;
    % .. Number of bootstraps
    matlabbatch{1}.spm.tools.swe.smodel.WB.WB_yes.WB_nB = 999;
    % .. Type of SwE (0 = U-SwE (recommended))
    matlabbatch{1}.spm.tools.swe.smodel.WB.WB_yes.WB_SwE = 0;
    % ... T or F contrast (CAVE: only one contrast at a time)
    if strcmp(MODEL,'bmi')
        if wild_con == 1
            matlabbatch{1}.spm.tools.swe.smodel.WB.WB_yes.WB_stat.WB_T.WB_T_con = [0 1 0 0 0 0];
        elseif wild_con == 2
            matlabbatch{1}.spm.tools.swe.smodel.WB.WB_yes.WB_stat.WB_T.WB_T_con = [0 0 1 0 0 0];
        end
    elseif strcmp(MODEL,'grouptime')
        if wild_con == 1
            matlabbatch{1}.spm.tools.swe.smodel.WB.WB_yes.WB_stat.WB_T.WB_T_con = [1 1 1 -1 -1 -1 0 0 0];
        elseif wild_con == 2
            matlabbatch{1}.spm.tools.swe.smodel.WB.WB_yes.WB_stat.WB_F.WB_F_con = [-1 0 1 -1 0 1 0 0 0;
                                                                                    0 -1 1 0 -1 1 0 0 0];
        elseif wild_con == 3
            matlabbatch{1}.spm.tools.swe.smodel.WB.WB_yes.WB_stat.WB_F.WB_F_con = [-1 0 1 1 0 -1 0 0 0; 
                                                                                    0 -1 1 0 1 -1 0 0 0];
        end
    end
    %  .. Inference Type (voxelwise, clusterwise, TFCE)
    if strcmp(INFERENCE_TYPE,'voxel')
        matlabbatch{1}.spm.tools.swe.smodel.WB.WB_yes.WB_infType.WB_voxelwise = 0;
    elseif strcmp(INFERENCE_TYPE,'cluster')
        % cluster-forming threshold (default)
        matlabbatch{1}.spm.tools.swe.smodel.WB.WB_yes.WB_infType.WB_clusterwise.WB_inputType.WB_img = 0;
        matlabbatch{1}.spm.tools.swe.smodel.WB.WB_yes.WB_infType.WB_clusterwise.WB_clusThresh = 0.001;
    elseif strcmp(INFERENCE_TYPE,'tfce')
        % E and H values as default (strongly recommended)
        matlabbatch{1}.spm.tools.swe.smodel.WB.WB_yes.WB_infType.WB_TFCE.WB_TFCE_E = 0.5;
        matlabbatch{1}.spm.tools.swe.smodel.WB.WB_yes.WB_infType.WB_TFCE.WB_TFCE_H = 2;
    end
end

%% Other ------------------------------------------------------------------
% Global calculation - Omit
matlabbatch{1}.spm.tools.swe.smodel.globalc.g_omit = 1;
% Global normalisation
%.Overall grand mean scaling - No
matlabbatch{1}.spm.tools.swe.smodel.globalm.gmsca.gmsca_no = 1;
% . Normalisation - None
matlabbatch{1}.spm.tools.swe.smodel.globalm.glonorm = 1;

%% Output folder ----------------------------------------------------------
location_SwE_mat = fullfile(out_folder, 'SwE.mat');

end

%% ------------------------------------------------------------------------
function [out_folder, exist_already] = create_out_folder(OUT_DIR, MODEL, one_roi_prep, COVARIATES, wild_con, INFERENCE_TYPE)
% Creates output folder and returns path to output folder as string.

if COVARIATES == 11
    model_name = 'Model1_group_timef_cage_sex_meanFD';
elseif COVARIATES == 12
    model_name = 'Model2_group_timef_cage_sex';
elseif COVARIATES == 13
    model_name = 'Model4_group_timec_cage_sex_meanFD';
elseif COVARIATES == 14
    model_name = 'Model5_group_timec_cage';
elseif COVARIATES == 21
    model_name = 'Model1_bmi_cage_sex_meanFD';
elseif COVARIATES == 22
    model_name = 'Model2_bmi_cage_sex';
end

if wild_con
    if strcmp(INFERENCE_TYPE,'voxel')
        model_name = strcat(model_name,'_WB_c0',num2str(wild_con),'_vox');
    elseif strcmp(INFERENCE_TYPE,'cluster')
        model_name = strcat(model_name,'_WB_c0',num2str(wild_con),'_cl');
    elseif strcmp(INFERENCE_TYPE,'tfce')
        model_name = strcat(model_name,'_WB_c0',num2str(wild_con),'_tfce');
    end 
end

if strcmp(MODEL,'bmi')
    parent_folder = 'Models_BMI';
elseif strcmp(MODEL,'grouptime')
    parent_folder = 'Models_GroupTime';
end

out_folder = fullfile(OUT_DIR, parent_folder, one_roi_prep, model_name);
if ~exist(out_folder, 'dir')
    exist_already = false;
    mkdir(out_folder)
else
    exist_already = true;
    return
end

end 

%% ------------------------------------------------------------------------
function cov = create_design_matrix(MODEL, COVARIATES)
% creates a design matrix cov as a struct compatible with the matlabbatch

if strcmp(MODEL,'bmi')
    % Prepare regressors
    [avgBMIc, cgnBMI] = swe_splitCovariate(readmatrix('BMI.txt'), readmatrix('subjNr'));
    if COVARIATES == 21 || COVARIATES == 22 
    % Intercept
    cov(1).c = ones(106,1);
    cov(1).cname = 'Intercept';

    % Covariates of interest
    cov(2).c = avgBMIc; cov(2).cname = 'avgBMI_centered';
    cov(3).c = cgnBMI; cov(3).cname = 'cgnBMI';

    % Nuisance Covariates
    age = readmatrix('Age.txt'); age = age - mean(age); % centered age
    cov(4).c = age; cov(4).cname = 'age';
    cov(5).c = readmatrix('Sex.txt'); cov(5).cname = 'sex';
        
        % Covariates (Design matrix, Model1)
        if COVARIATES == 21
            cov(6).c = readmatrix('logmeanFD.txt'); cov(6).cname = 'meanFD';
        end
    end
 
elseif strcmp(MODEL,'grouptime')
    % Intercept is implicitely modeled by the group means
    if COVARIATES == 11 || COVARIATES == 12 
        %% Covariates (Design matrix for factorial time)
        % IG bl, fu and fu2
        cov(1).c = readmatrix('IG_bl.txt'); cov(1).cname = 'IG bl';
        cov(2).c = readmatrix('IG_fu.txt'); cov(2).cname = 'IG fu';
        cov(3).c = readmatrix('IG_fu2.txt'); cov(3).cname = 'IG fu2';
        % KG bl, fu and fu2
        cov(4).c = readmatrix('KG_bl.txt'); cov(4).cname = 'KG bl';
        cov(5).c = readmatrix('KG_fu.txt'); cov(5).cname = 'KG fu';
        cov(6).c = readmatrix('KG_fu2.txt'); cov(6).cname = 'KG fu2';
        
        % Covariates (Design matrix)
        age = readmatrix('Age.txt'); age = age - mean(age); % centered age
        cov(7).c = age; cov(7).cname = 'age';
        cov(8).c = readmatrix('Sex.txt'); cov(8).cname = 'sex';
        if COVARIATES == 11
            % Covariates (Design matrix, Model1)
            cov(9).c = readmatrix('logmeanFD.txt'); cov(9).cname = 'meanFD';
        end
    elseif COVARIATES == 13 || COVARIATES == 14 
        %% Covariates (Design matrix for continuous time)
        % IG
        cov(1).c = readmatrix('group_IG.txt'); cov(1).cname = 'IG';
        % KG
        cov(2).c = readmatrix('group_KG.txt'); cov(2).cname = 'KG';
        % time points IG
        cov(3).c = readmatrix('tp_IG.txt'); cov(3).cname = 'time_c_IG';
        % time points KG
        cov(4).c = readmatrix('tp_KG.txt'); cov(4).cname = 'time_c_KG'; 
        
        % Covariates (Design matrix)
        cov(5).c = readmatrix('Age.txt'); cov(5).cname = 'age';
        cov(6).c = readmatrix('Sex.txt'); cov(6).cname = 'sex';
        if COVARIATES == 13
            % Covariates (Design matrix, Model3)
            cov(7).c = readmatrix('logmeanFD.txt'); cov(7).cname = 'meanFD';
        end
    end
end

end

%% ------------------------------------------------------------------------
function matlabbatch = RunModel(matlabbatch, location_SwE_mat)
% Estimates Model as specified in matlabbatch.
% Path to SwE.mat needs to be specified by location_SwE_mat.
% If matlabbatch defines a model for non-parametric estimation, the
% function will estimate the contrast specified in matlabbatch.
% If matlabbatch defines a model for parametric estimation, the function
% will only compute beta estimates but not contrasts.
matlabbatch{2}.spm.tools.swe.rmodel.des = {location_SwE_mat};
end

%% ------------------------------------------------------------------------
function [matlabbatch] = DisplayResults(location_SwE_mat)
% Displays results of estimated model as saved in SwE.mat file. 
% Path to SwE.mat needs to be specified by location_SwE_mat.
matlabbatch{1}.spm.tools.swe.results.dis = {location_SwE_mat};
end

