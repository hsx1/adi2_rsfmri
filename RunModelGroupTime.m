function RunModelGroupTime(param)

% loop over all models, only all covariates, all inference types for
% wild_boot = true
for i = 1:length(param.MODEL)
    %covariates = COVARIATES;
    for j = 1:length(param.COVARIATES)
        for k = 1:length(param.INFERENCE_TYPE)
            if param.WILD_BOOT inference = param.INFERENCE_TYPE{k}; else inference = 'parametric'; end
            fprintf('Evaluate model %s with covariates %i, %s inference...\n', param.MODEL{i},param.COVARIATES(j),inference)
            RunSwe(param,i,j,k);
        end
    end
end
spm('Quit')

%% ========================================================================
function display_message(param)
% displays message according to the COVARIATES, that define the model
if strcmp(param.MODEL, 'grouptime')
    if param.COVARIATES == 11
        fprintf('%s\n',...
        'Model1_group_timef_cage_sex_meanFD: In case of parametric estimation',...
        'please type in the following contrasts for all runs',...
        'group [1 1 1 -1 -1 -1 0 0 0]',...
        'time [-1 0 1 -1 0 1 0 0 0; 0 -1 1 0 -1 1 0 0 0]',...
        'group*time [-1 0 1 1 0 -1 0 0 0; 0 -1 1 0 1 -1 0 0 0]')
    elseif param.COVARIATES == 12
        fprintf('%s\n',...
        'Model2_group_timef_cage_sex: In case of parametric estimation',...
        'please type in the following contrasts for all runs',...
        'group [1 1 1 -1 -1 -1 0 0]',...
        'time [-1 0 1 -1 0 1 0 0; 0 -1 1 0 -1 1 0 0]',...
        'group*time [-1 0 1 1 0 -1 0 0; 0 -1 1 0 1 -1 0 0]')
    elseif param.COVARIATES == 13
        fprintf('%s\n',...
        'Model3_group_timec_cage_sex_meanFD: In case of parametric estimation',...
        'please type in the following contrasts for all runs',...
        'group [0 1 -1 0 0 0 0 0]',...
        'time [0 0 0 0.5 0.5 0 0 0]',...
        'group*time [0 0 0 1 -1 0 0 0]')
    elseif param.COVARIATES == 14
        fprintf('%s\n',...
        'Model4_group_timec_cage_sex: In case of parametric estimation',...
        'please type in the following contrasts for all runs',...
        'group [0 1 -1 0 0 0 0]',...
        'time [0 0 0 0.5 0.5 0 0]',...
        'group*time [0 0 0 1 -1 0 0]')
    end
elseif strcmp(param.MODEL, 'grouptime2tp')
    if param.COVARIATES == 11
        fprintf('%s\n',...
        'Model1_group_timef_cage_sex_meanFD: In case of parametric estimation',...
        'please type in the following contrasts for all runs',...
        'group [1 1 -1 -1 0 0 0]',...
        'time [-1 1 -1 1 0 0 0]',...
        'group*time [-1 1 1 -1 0 0 0]')
    elseif param.COVARIATES == 12
        fprintf('%s\n',...
        'Model2_group_timef_cage_sex: In case of parametric estimation',...
        'please type in the following contrasts for all runs',...
        'group [1 1 -1 -1 0 0 0]',...
        'time [-1 1 -1 1 0 0]',...
        'group*time [-1 1 1 -1 0 0]') 
    end
end

end

%% ------------------------------------------------------------------------
function RunSwe(param,i,j,k)
% Carries out analysis according to specified parameters. 
% RunSwe makes use of three functions: SpecifyModel, RunModel and
% Display Results. If SpecifyModel and RunModel has run once and the
% SwE.mat file is 
param.MODEL = param.MODEL{i};
param.COVARIATES = param.COVARIATES(j);
param.INFERENCE_TYPE = param.INFERENCE_TYPE(k);


if param.ONLY_DISPLAY && not(param.WILD_BOOT)
    display_message(param)
    % wait so displayed information can be read
    pause(7)
end

if strcmp(param.MODEL,'grouptime2tp')
    param.INFO_DIR = fullfile(param.INFO_DIR, 'only2tp');
else
    param.INFO_DIR = fullfile(param.INFO_DIR, 'total');
end

nrun = length(param.ROI_PREP);
for crun = 1:nrun
    % estimate any of the six models for every roi_prep in ROI_PREP
    if not(param.WILD_BOOT)
        % parametric estimation
        param.wild_con = false;
        % check if analysis has been conducted already
        [out_folder, exist_already] = create_out_folder(param, crun);
        fprintf('If necessary create output folder \n %s ...\n',out_folder)
        clear matlabbatch
        fprintf('Specify Model...\n')
        [matlabbatch,location_SwE_mat] = SpecifyModel(param, crun, out_folder);
        if param.ONLY_DISPLAY
            fprintf('Display Results...\n')
            spm('defaults', 'FMRI'); 
            matlabbatch = DisplayResults(location_SwE_mat);
            spm_jobman('run', matlabbatch);
            pause(param.VIEWSEC)
            save("xSwE.mat")
        elseif not(exist_already) || param.OVERWRITE
            fprintf('Estimate model...\n')
            spm('defaults', 'FMRI'); 
            matlabbatch = RunModel(matlabbatch, location_SwE_mat);
            spm_jobman('run', matlabbatch);
        else
            fprintf('... model already estimated or error in parameter definitions.\n')
        end
        % manually type in contrasts (issue: see below)
    else
        % non-parametric estimation with Wild-Bootstrap
        for wild_con = 1:3
            param.wild_con = wild_con;
            % check if analysis has run already
            [out_folder, exist_already] = create_out_folder(param, crun);
            fprintf('If necessary create output folder \n %s ...\n',out_folder)
            clear matlabbatch
            fprintf('Specify Model...\n')
            [matlabbatch, location_SwE_mat] = SpecifyModel(param, crun, out_folder); 
            if param.ONLY_DISPLAY
                fprintf('Display Results...\n')
                spm('defaults', 'FMRI');
                matlabbatch = DisplayResults(location_SwE_mat);
                spm_jobman('run', matlabbatch);
                pause(param.VIEWSEC)
                save("xSwE.mat")
            elseif not(exist_already)
                fprintf('Estimate model...\n')
                spm('defaults', 'FMRI');
                matlabbatch = RunModel(matlabbatch, location_SwE_mat);
                spm_jobman('run', matlabbatch);
            else
                fprintf('... model already estimated or error in parameter definitions.\n')
            end
        end
    end
 
% save results of contrasts in folder roi_prep and subfolder modelx
% currently not possible to automate contrasts and save results
% https://github.com/NISOx-BDI/SwE-toolbox/issues/135
end
fprintf('...done.')

end

%% ------------------------------------------------------------------------
function [matlabbatch,location_SwE_mat] = SpecifyModel(param, crun, out_folder)
% Specifies design matrix and estimation procedure as preparation for the
% estimation to be carried out (by function RunModel); the specified model 
% is stored in SwE.mat file in location_SwE_mat.

%% Directory to get information for model definition
cd(param.INFO_DIR)

%% Directory to save SwE file
% save directory of (new) folder
matlabbatch{1}.spm.tools.swe.smodel.dir = {out_folder};

%% Load Scans -------------------------------------------------------------
% constructing cell array with the scans (FC maps per subject and time point)
scans_dir = readcell(fullfile(param.INFO_DIR,'scans.txt'), 'Delimiter',' ','Whitespace',"'");
roi_prep = param.ROI_PREP{crun};
if contains(roi_prep, 'gsr')
    % special case if preprocessing is 'gsr'
    if endsWith(roi_prep,'_z')
        FC_files = strcat(lower(roi_prep(1:end-5)),'cc_gsr_seed_correlation_z_trans.nii');
    else
        FC_files = strcat(lower(roi_prep(1:end-3)),'cc_gsr_seed_correlation_trans.nii');
    end
else
    if endsWith(param.ROI_PREP{crun},'_z')
        FC_files = strcat(lower(roi_prep(1:end-2)),'_seed_correlation_z_trans.nii');
    else
        FC_files = strcat(lower(roi_prep),'_seed_correlation_trans.nii');
    end
end
scans_of_roi = fullfile(scans_dir, roi_prep, FC_files);

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
cov = create_design_matrix(param);
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
if strcmp(param.MASK,'brain')
    mask_path = {fullfile(param.MASK_DIR,'MNI_resampled_brain_mask.nii,1')};
elseif strcmp(param.MASK,'gm')
    mask_path = {fullfile(param.MASK_DIR,'mni_icbm152_gm_tal_nlin_sym_09a.nii,1')};
end

matlabbatch{1}.spm.tools.swe.smodel.masking.em = mask_path;
%% Non-parametric Wild Bootstrap ------------------------------------------
% . No
if ~param.wild_con
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
        
    if strcmp(param.MODEL, 'grouptime2tp')
        c01 = [1 1 -1 -1 0 0 0];
        c02 = [-1 1 -1 1 0 0 0];
        c03 = [-1 1 1 -1 0 0 0];
    else 
        c01 = [1 1 1 -1 -1 -1 0 0 0];
        c02 = [-1 0 1 -1 0 1 0 0 0; 0 -1 1 0 -1 1 0 0 0];
        c03 = [-1 0 1 1 0 -1 0 0 0; 0 -1 1 0 1 -1 0 0 0];
    end
    
    % if model without covariates, shorten contrasts
    if (param.COVARIATES == 12) s = 1; else s = 0; end
    
    if param.wild_con == 1
        matlabbatch{1}.spm.tools.swe.smodel.WB.WB_yes.WB_stat.WB_T.WB_T_con = c01(1:end-s);
    elseif param.wild_con == 2
        matlabbatch{1}.spm.tools.swe.smodel.WB.WB_yes.WB_stat.WB_F.WB_F_con = c02(:,1:end-s);
    elseif param.wild_con == 3
        matlabbatch{1}.spm.tools.swe.smodel.WB.WB_yes.WB_stat.WB_F.WB_F_con = c03(:,1:end-s);
    end
    %  .. Inference Type (voxelwise, clusterwise, TFCE)
    if strcmp(param.INFERENCE_TYPE,'voxel')
        matlabbatch{1}.spm.tools.swe.smodel.WB.WB_yes.WB_infType.WB_voxelwise = 0;
    elseif strcmp(param.INFERENCE_TYPE,'cluster')
        % cluster-forming threshold (default)
        matlabbatch{1}.spm.tools.swe.smodel.WB.WB_yes.WB_infType.WB_clusterwise.WB_inputType.WB_img = 0;
        matlabbatch{1}.spm.tools.swe.smodel.WB.WB_yes.WB_infType.WB_clusterwise.WB_clusThresh = 0.001;
    elseif strcmp(param.INFERENCE_TYPE,'tfce')
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
function [out_folder, exist_already] = create_out_folder(param, crun)
% Creates output folder and returns path to output folder as string.

if strcmp(param.MODEL,'grouptime2tp')
    parent_folder = 'GroupTime_only2tp';
else
    parent_folder = 'GroupTime_total';
end


if strcmp(param.MASK, 'gm')
    mask_def = 'gm';
else
    mask_def = 'brain';
end

if param.COVARIATES == 11
    model_name = 'group-timef-age-sex-meanFD';
elseif param.COVARIATES == 12
    model_name = 'group-timef-age-sex';
elseif param.COVARIATES == 13
    model_name = 'group-timec-age-sex-meanFD';
elseif param.COVARIATES == 14
    model_name = 'group-timec-age';
end

if param.wild_con % Wild Bootstrap
    if strcmp(param.INFERENCE_TYPE,'voxel')
        model_name = strcat(model_name,'_WB-c0',num2str(param.wild_con),'vox');
    elseif strcmp(param.INFERENCE_TYPE,'cluster')
        model_name = strcat(model_name,'_WB-c0',num2str(param.wild_con),'cl');
    elseif strcmp(param.INFERENCE_TYPE,'tfce')
        model_name = strcat(model_name,'_WB-c0',num2str(param.wild_con),'tfce');
    end 
else % Parametric Estimation
    model_name = strcat(model_name,'_PE-all');
end

% Create folder
out_folder = fullfile(param.OUT_DIR, parent_folder, mask_def, param.ROI_PREP{crun}, model_name);
if ~exist(out_folder, 'dir')
    exist_already = false;
    mkdir(out_folder)
else
    exist_already = true;
    return
end

end 

%% ------------------------------------------------------------------------
function cov = create_design_matrix(param)
% creates a design matrix cov as a struct compatible with the matlabbatch

% Prepare regressors (regressor count = r)
r = 1;
% Intercept is implicitely modeled by the group means
if param.COVARIATES == 11|| param.COVARIATES == 12
    %% Covariates (Design matrix for factorial time)
    % IG bl, fu and fu2
    cov(r).c = readmatrix('IG_bl.txt'); cov(r).cname = 'IG bl'; r = r+1;
    cov(r).c = readmatrix('IG_fu.txt'); cov(r).cname = 'IG fu'; r = r+1;
    if strcmp(param.MODEL,'grouptime')
        cov(r).c = readmatrix('IG_fu2.txt'); cov(r).cname = 'IG fu2'; r = r+1;
    end
    % KG bl, fu and fu2
    cov(r).c = readmatrix('KG_bl.txt'); cov(r).cname = 'KG bl'; r = r+1;
    cov(r).c = readmatrix('KG_fu.txt'); cov(r).cname = 'KG fu'; r = r+1;
    if strcmp(param.MODEL,'grouptime')
        cov(r).c = readmatrix('KG_fu2.txt'); cov(r).cname = 'KG fu2'; r = r+1;
    end

    % Covariates (Design matrix)
    age = readmatrix('Age.txt'); age = age - mean(age); % centered age
    cov(r).c = age; cov(r).cname = 'age'; r = r+1;
    cov(r).c = readmatrix('Sex.txt'); cov(r).cname = 'sex'; r = r+1;
    if param.COVARIATES == 11
        % Covariates (Design matrix, Model1)
        cov(r).c = readmatrix('logmeanFD.txt'); cov(r).cname = 'logmeanFD'; r = r+1;
    end
elseif param.COVARIATES == 13 || param.COVARIATES == 14 
    %% Covariates (Design matrix for continuous time)
    % IG and KG
    cov(r).c = readmatrix('group_IG.txt'); cov(r).cname = 'IG'; r = r+1;
    cov(r).c = readmatrix('group_KG.txt'); cov(r).cname = 'KG'; r = r+1;
    % time points IG and time points KG
    cov(r).c = readmatrix('tp_IG.txt'); cov(r).cname = 'time_c_IG'; r = r+1;
    cov(r).c = readmatrix('tp_KG.txt'); cov(r).cname = 'time_c_KG'; r = r+1;
        
    % Covariates (Design matrix)
    cov(r).c = readmatrix('Age.txt'); cov(r).cname = 'age'; r = r+1;
    cov(r).c = readmatrix('Sex.txt'); cov(r).cname = 'sex'; r = r+1;
    if param.COVARIATES == 13
        % Covariates (Design matrix, Model3)
        cov(r).c = readmatrix('logmeanFD.txt'); cov(r).cname = 'meanFD'; r = r+1;
    end
end
r = 0;
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

end

