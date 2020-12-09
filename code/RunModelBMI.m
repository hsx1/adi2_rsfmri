function RunModelBMI(param)

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
fprintf("Enter > spm('Quit') < to exit.")
%spm('Quit')

%% ========================================================================
function display_message(param)
% displays message according to the COVARIATES, that define the model
if param.COVARIATES == 21
    fprintf('%s\n',...
    'Model1_BMI_meanFD_cage_sex: In case of parametric estimation',...
    'please type in the following contrasts for all runs',...
    'avgBMI [0 1 0 0 0 0]',...
    'BMIcgn [0 0 1 0 0 0]')
elseif param.COVARIATES == 22
    fprintf('%s\n',...
    'Model2_BMI_cage_sex: In case of parametric estimation',...
    'please type in the following contrasts for all runs',...
    'avgBMI [0 1 0 0 0]',...
    'BMIcgn [0 0 1 0 0]')
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


if param.ONLY_DISPLAY && not(param.WILD_BOOT) %strcmp(param.ACTION,'display')
    display_message(param)
    % wait so displayed information can be read
    pause(7)
end

if param.EXCLFD==true
    param.INFO_DIR = fullfile(param.INFO_DIR, 'ExclFD');
else
    param.INFO_DIR = fullfile(param.INFO_DIR, 'noExclFD');
end

if strcmp(param.MODEL,'bmiIG')
    param.INFO_DIR = fullfile(param.INFO_DIR, 'IG/total');
elseif strcmp(param.MODEL,'bmi2tp')
    param.INFO_DIR = fullfile(param.INFO_DIR, 'both/BLFU');
else
    param.INFO_DIR = fullfile(param.INFO_DIR, 'both/total');
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
        if param.ONLY_DISPLAY %strcmp(param.ACTION,'display')
            fprintf('Display Results...\n')
            spm('defaults', 'FMRI'); 
            matlabbatch = DisplayResults(location_SwE_mat);
            spm_jobman('run', matlabbatch);
            pause(param.VIEWSEC)
        elseif not(exist_already) || param.OVERWRITE %strcmp(param.ACTION,'overwrite')
            fprintf('Estimate model...\n')
            % delete existing files in folder if existent
            if exist_already 
                rmdir(out_folder,'s'); % delete former dir
                mkdir(out_folder); % create new empty one
            end
            spm('defaults', 'FMRI'); 
            matlabbatch = RunModel(matlabbatch, location_SwE_mat);
            spm_jobman('run', matlabbatch);
        else
            fprintf('... model already estimated or error in parameter definitions.\n')
        end
        % manually type in contrasts (issue: see below)
    else
        % non-parametric estimation with Wild-Bootstrap
        for wild_con = 1:2
            param.wild_con = wild_con;
            % check if analysis has run already
            [out_folder, exist_already] = create_out_folder(param, crun);
            fprintf('If necessary create output folder \n %s ...\n',out_folder)
            clear matlabbatch
            fprintf('Specify Model...\n')
            [matlabbatch, location_SwE_mat] = SpecifyModel(param, crun, out_folder); 
            if param.ONLY_DISPLAY %strcmp(param.ACTION,'display')
                fprintf('Display Results...\n')
                spm('defaults', 'FMRI');
                matlabbatch = DisplayResults(location_SwE_mat);
                spm_jobman('run', matlabbatch);
                pause(param.VIEWSEC)
                save("xSwE.mat")
            elseif not(exist_already)|| param.OVERWRITE %strcmp(param.ACTION,'overwrite')
                fprintf('Estimate model...\n')
                % delete existing files in folder if existent
                if exist_already 
                    rmdir(out_folder,'s'); % delete former dir
                    mkdir(out_folder); % create new empty one
                end
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
fprintf('... done.')

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
smodel.dir = {out_folder};

% cifti + gifti additional information (SwE 2.2.1)
smodel.ciftiAdditionalInfo.ciftiGeomFile = struct('brainStructureLabel', {}, 'geomFile', {}, 'areaFile', {});
smodel.ciftiAdditionalInfo.volRoiConstraint = 0;
smodel.giftiAdditionalInfo.areaFileForGiftiInputs = {};

%% Load Scans -------------------------------------------------------------
% constructing cell array with the scans (FC maps per subject and time point)
scans_dir = readcell(fullfile(param.INFO_DIR,'scans.txt'), 'Delimiter',' ','Whitespace',"'");
roi_prep = param.ROI_PREP{crun};
scans_of_roi = create_scans_list(scans_dir, roi_prep);

% load scans to matlabbatch
smodel.scans = scans_of_roi;
                                         
%% SwE type ---------------------------------------------------------------
% .Modified
% .. Define Groups
smodel.type.modified.groups = readmatrix('group.txt');

% ..Visits
smodel.type.modified.visits = readmatrix('tp.txt');

% small sample adjustment (4 = type C2)
smodel.type.modified.ss = 4;
% degrees of freedom type (3 = approx III)
smodel.type.modified.dof_mo = 3;

%% Subjects ---------------------------------------------------------------
smodel.subjects = readmatrix('subjNr.txt');

%% Covariates (Design matrix) ---------------------------------------------
% desgin matrix for model estimation
cov = create_design_matrix(param);
smodel.cov = cov;
  
%% Masking ----------------------------------------------------------------
%% Multiple Covariates (none)
smodel.multi_cov = struct('files', {});

%% Masking (none)
% .Threshold masking
smodel.masking.tm.tm_none = 1;
% ..Implicit Mask (yes)
smodel.masking.im = 1;
% .. Explicit Mask
if strcmp(param.MASK,'brain')
    mask_path = {fullfile(param.MASK_DIR, param.MASK_B)};
elseif strcmp(param.MASK,'gm')
    mask_path = {fullfile(param.MASK_DIR, param.MASK_GM)};
end

smodel.masking.em = mask_path;
%% Non-parametric Wild Bootstrap ------------------------------------------
% . No
if ~param.wild_con
    smodel.WB.WB_no = 0;
else
    % . Yes
    % .. Small sample adjustments for WB resampling (4 = type C2)
    smodel.WB.WB_yes.WB_ss = 4;
    % .. Number of bootstraps
    smodel.WB.WB_yes.WB_nB = 1000;
    % .. Type of SwE (0 = U-SwE (recommended))
    smodel.WB.WB_yes.WB_SwE = 0;
    % ... T or F contrast (CAVE: only one contrast at a time)
    c01 = [0 1 0 0 0 0];
    c02 = [0 0 1 0 0 0];
    % if model without covariates, shorten contrasts
    s = 0;
    if param.COVARIATES == 22
        s = 1;
    elseif param.COVARIATES == 23
        s = 3;
    end
    
    if param.wild_con == 1
        smodel.WB.WB_yes.WB_stat.WB_T.WB_T_con = c01(1:end-s);
    elseif param.wild_con == 2
        smodel.WB.WB_yes.WB_stat.WB_T.WB_T_con = c02(1:end-s);
    end
    %  .. Inference Type (voxelwise, clusterwise, TFCE)
    if strcmp(param.INFERENCE_TYPE,'voxel')
        smodel.WB.WB_yes.WB_infType.WB_voxelwise = 0;
    elseif strcmp(param.INFERENCE_TYPE,'cluster')
        % cluster-forming threshold (default)
        smodel.WB.WB_yes.WB_infType.WB_clusterwise.WB_clusThresh = 0.001;
        smodel.WB.WB_yes.WB_infType.WB_clusterwise.WB_inputType.WB_img = 0;
    elseif strcmp(param.INFERENCE_TYPE,'tfce')
        % E and H values as default (strongly recommended)
        smodel.WB.WB_yes.WB_infType.WB_TFCE.WB_TFCE_E = 0.5;
        smodel.WB.WB_yes.WB_infType.WB_TFCE.WB_TFCE_H = 2;
    end
end

%% Other ------------------------------------------------------------------
% Global calculation - Omit
smodel.globalc.g_omit = 1;
% Global normalisation
%.Overall grand mean scaling - No
smodel.globalm.gmsca.gmsca_no = 1;
% . Normalisation - None
smodel.globalm.glonorm = 1;

%% save model specification in matlabbatch
matlabbatch{1}.spm.tools.swe.smodel = smodel;

%% Output folder ----------------------------------------------------------
location_SwE_mat = fullfile(out_folder, 'SwE.mat');

end

%% ------------------------------------------------------------------------
function [out_folder, exist_already] = create_out_folder(param, crun)
% Creates output folder and returns path to output folder as string.

if param.EXCLFD==true
    excl = 'ExclFD';
else
    excl = 'noExclFD';
end

if strcmp(param.MODEL,'bmiIG')
    parent_folder = 'BMI_onlyIG';
elseif strcmp(param.MODEL,'bmi2tp')
    parent_folder = 'BMI_only2tp';
else
    parent_folder = 'BMI_total';
end

if strcmp(param.MASK, 'gm')
    mask_def = 'gm';
else
    mask_def = 'brain';
end

if param.COVARIATES == 21
    model_name = 'bmi-age-sex-meanFD';
elseif param.COVARIATES == 22
    model_name = 'bmi-age-sex';
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
out_folder = fullfile(param.OUT_DIR, excl, parent_folder, mask_def, param.ROI_PREP{crun}, model_name);
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
[avgBMIc, cgnBMI] = swe_splitCovariate(readmatrix('BMI.txt'), readmatrix('subjNr'));

% Intercept
length_intercept = length(readmatrix('subjNr'));
cov(r).c = ones(length_intercept,1); cov(r).cname = 'Intercept'; r = r+1;
    
% Covariates of interest    
if param.COVARIATES == 21 || param.COVARIATES == 22 
    cov(r).c = avgBMIc; cov(r).cname = 'avgBMI_centered'; r = r+1;
    cov(r).c = cgnBMI; cov(r).cname = 'cgnBMI'; r = r+1;

    % Nuisance Covariates
    age = readmatrix('Age.txt'); age = age - mean(age); % centered age
    cov(r).c = age; cov(r).cname = 'age'; r = r+1;
    cov(r).c = readmatrix('Sex.txt'); cov(r).cname = 'sex'; r = r+1;
        
    % Covariates (Design matrix, Model1)
    if param.COVARIATES == 21
        cov(r).c = readmatrix('logmeanFD.txt'); cov(r).cname = 'meanFD'; r = r+1;
    end
end
r = 0;
end

%% ------------------------------------------------------------------------
function [scans_of_roi] = create_scans_list(scans_dir, roi_prep)
% creates a list of scans as input for swe model specification "Scans" from
% directory of scans and regions of interest
if contains(roi_prep, 'gsr')
    % special case if preprocessing is 'gsr'
    if endsWith(roi_prep,'_z')
        FC_files = strcat(lower(roi_prep(1:end-5)),'cc_gsr_seed_correlation_z_trans.nii');
    else
        FC_files = strcat(lower(roi_prep(1:end-3)),'cc_gsr_seed_correlation_trans.nii');
    end
else
    if endsWith(roi_prep,'_z')
        FC_files = strcat(lower(roi_prep(1:end-2)),'_seed_correlation_z_trans.nii');
    else
        FC_files = strcat(lower(roi_prep),'_seed_correlation_trans.nii');
    end
end
scans_of_roi = fullfile(scans_dir, roi_prep, FC_files);
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
