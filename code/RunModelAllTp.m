%% ========================================================================
function RunModelAllTp(crun)
% Mixed effect test for network activation.
% Analysis is carried out with SwE toolbox. Results can be displayed as
% from other models (see README).
    
RunSwe(crun);

%% ========================================================================
function display_message(COVARIATES)
% displays message according to the COVARIATES, that define the model
if COVARIATES == 41
    fprintf('%s\n',...
    'fsthdbt',...
    '...')
elseif COVARIATES == 42
    fprintf('%s\n',...
    'jzdizd',...
    '...')
end

end

%% ------------------------------------------------------------------------
function RunSwe(crun)
% Carries out analysis according to specified parameters. 
% RunSwe makes use of three functions: SpecifyModel, RunModel and
% Display Results. If SpecifyModel and RunModel has run once and the
% SwE.mat file is 

if crun.ONLY_DISPLAY && not(crun.WILD_BOOT) %strcmp(crun.ACTION,'display')
    display_message(crun.COVARIATES)
    % wait so displayed information can be read
    pause(3)
end

if crun.EXCLFD==true
    crun.INFO_DIR = fullfile(crun.INFO_DIR, 'ExclFD','both/total/');
else
    crun.INFO_DIR = fullfile(crun.INFO_DIR, 'noExclFD','both/total/');
end

if ~crun.WILD_BOOT 
    crun.INFERENCE_TYPE = 'parametric'; 
    % only parametric estimation
    crun.wild_con = false;
end
fprintf('Evaluate model %s with covariates %i, %s inference...\n', crun.MODEL,crun.COVARIATES,crun.INFERENCE_TYPE)

crun.wild_con = 1;
for i = crun.wild_con
    crun.wild_con = i;
    % check if analysis has run already
    [out_folder, exist_already] = create_out_folder(crun);
    fprintf('If necessary create output folder \n %s ...\n', out_folder)
      clear matlabbatch
    fprintf('Specify Model...\n')
      [matlabbatch, location_SwE_mat] = SpecifyModel(crun, out_folder);
    if crun.ONLY_DISPLAY %strcmp(crun.ACTION,'display')
        fprintf('Display Results...\n')
        cd(out_folder)
        [hReg,xSwE,SwE] = swe_results_ui("Setup");
        spm2csv(hReg,xSwE);
        if crun.VIEW
            h=questdlg('Please press OK to continue or Pause to inspect.','Continue','OK','Pause','OK');
            switch h
                case 'OK'
                    spm('Quit')
                    continue
                case 'Pause'
                    b = inputdlg('How long to you want to pause (sec)?','Pause',1,{'60'}); 
                    b = str2double(b);
                    pause(b)
                otherwise
                    continue
            end 
        end
    elseif not(exist_already) || crun.OVERWRITE %strcmp(crun.ACTION,'overwrite')
        fprintf('Estimate model...\n')
        % delete existing files in folder if existent
        if exist_already
            rmdir(out_folder, 's'); % delete former dir
            mkdir(out_folder); % create new empty one
        end
        spm('defaults', 'FMRI');
        matlabbatch = RunModel(matlabbatch, location_SwE_mat);
        spm_jobman('run', matlabbatch);
    else
        fprintf('... model already estimated.\n')
    end
end 
% save results of contrasts in folder roi_prep and subfolder modelx
% currently not possible to automate contrasts and save results
% https://github.com/NISOx-BDI/SwE-toolbox/issues/135

fprintf('...done.')

end

%% ------------------------------------------------------------------------
function [matlabbatch,location_SwE_mat] = SpecifyModel(crun, out_folder)
% Specifies design matrix and estimation procedure as preparation for the
% estimation to be carried out (by function RunModel); the specified model 
% is stored in SwE.mat file in location_SwE_mat.

%% Directory to get information for model definition
cd(crun.INFO_DIR)

%% Directory to save SwE file
% save directory of (new) folder
smodel.dir = cellstr(out_folder);

% cifti + gifti additional information (SwE 2.2.1)
smodel.ciftiAdditionalInfo.ciftiGeomFile = struct('brainStructureLabel', {}, 'geomFile', {}, 'areaFile', {});
smodel.ciftiAdditionalInfo.volRoiConstraint = 0;
smodel.giftiAdditionalInfo.areaFileForGiftiInputs = {};

%% Load Scans -------------------------------------------------------------
% constructing cell array with the scans (FC maps per subject and time point)
scans_dir = readcell(fullfile(crun.INFO_DIR,'scans.txt'), 'Delimiter',' ','Whitespace',"'");

scans_of_roi = create_scans_list(scans_dir, crun.ROI_PREP);

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
cov = create_design_matrix(crun);
smodel.cov = cov;
  
%% Masking ----------------------------------------------------------------
%% Multiple Covariates (none)
smodel.multi_cov = struct('files', {});

%% Masking (none)
% .Threshold masking
smodel.masking.tm.tm_none = 1;
% ..Implicit Mask (no)
smodel.masking.im = 1;
% .. Explicit Mask
if strcmp(crun.MASK,'brain')
    mask_path = cellstr(fullfile(crun.MASK_DIR, crun.MASK_B));
elseif strcmp(crun.MASK,'gm')
    mask_path = cellstr(fullfile(crun.MASK_DIR, crun.MASK_GM));
end

smodel.masking.em = mask_path;
%% Non-parametric Wild Bootstrap ------------------------------------------
% . No
if ~crun.wild_con
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
    c01 = [1 0 0 0];
    % if model without covariates, shorten contrasts
    s = 0;
    if crun.COVARIATES == 42
        s = 1;
    elseif crun.COVARIATES == 43
        s = 3;
    end
    
    if crun.wild_con == 1
        smodel.WB.WB_yes.WB_stat.WB_T.WB_T_con = c01(1:end-s);
    end
    %  .. Inference Type (voxelwise, clusterwise, TFCE)
    if strcmp(crun.INFERENCE_TYPE,'voxel')
        smodel.WB.WB_yes.WB_infType.WB_voxelwise = 0;
    elseif strcmp(crun.INFERENCE_TYPE,'cluster')
        % cluster-forming threshold (default)
        smodel.WB.WB_yes.WB_infType.WB_clusterwise.WB_clusThresh = 0.001;
        smodel.WB.WB_yes.WB_infType.WB_clusterwise.WB_inputType.WB_img = 0;
    elseif strcmp(INFERENCE_TYPE,'tfce')
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
smodel.globalm.gmsca.gmsca_no = 0;
% . Normalisation - None
smodel.globalm.glonorm = 1;

%% save model specification in matlabbatch
matlabbatch{1}.spm.tools.swe.smodel = smodel;

%% Output folder ----------------------------------------------------------
location_SwE_mat = fullfile(out_folder, 'SwE.mat');

end

%% ------------------------------------------------------------------------
function [out_folder, exist_already] = create_out_folder(crun)
% Creates output folder and returns path to output folder as string.

if crun.EXCLFD==true
    excl = 'ExclFD';
else
    excl = 'noExclFD';
end

parent_folder = 'Network_0-6-12';

% if strcmp(crun.MASK, 'gm')
%     mask_def = 'gm';
% else
%     mask_def = 'brain';
% end

if crun.COVARIATES == 41
    model_name = 'network-fd-age-sex'; %previously 'network-age-sex'
elseif crun.COVARIATES == 42
    model_name = 'network-age-sex';%previously 'network-fd-age-sex'
elseif crun.COVARIATES == 43
    model_name = 'network';
end

if crun.wild_con % Wild Bootstrap
    if strcmp(crun.INFERENCE_TYPE,'voxel')
        model_name = strcat(model_name,'_WB-c0',num2str(crun.wild_con),'vox');
    elseif strcmp(crun.INFERENCE_TYPE,'cluster')
        model_name = strcat(model_name,'_WB-c0',num2str(crun.wild_con),'cl');
    elseif strcmp(crun.INFERENCE_TYPE,'tfce')
        model_name = strcat(model_name,'_WB-c0',num2str(crun.wild_con),'tfce');
    end 
else % Parametric Estimation
    model_name = strcat(model_name,'_PE-all');
end

% Create folder
out_folder = fullfile(crun.OUT_DIR, excl, parent_folder, crun.MASK, crun.ROI_PREP, model_name);
if ~exist(out_folder, 'dir')
    exist_already = false;
    mkdir(out_folder)
else
    exist_already = true;
    return
end

end 

%% ------------------------------------------------------------------------
function cov = create_design_matrix(crun)
% Creates a design matrix cov as a struct compatible with the matlabbatch

% Prepare regressors

% Prepare regressors (regressor count = r)
r = 1;
% Intercept
length_intercept = length(readmatrix('subjNr'));
cov(r).c = ones(length_intercept,1); cov(r).cname = 'Intercept'; r = r+1;

% Covariates of interest  
if crun.COVARIATES == 41
    % Nuisance Covariates
    cov(r).c = readmatrix('logmeanFD.txt'); cov(r).cname = 'logmeanFD'; r = r+1;
end

if crun.COVARIATES == 41 || crun.COVARIATES == 42
    age = readmatrix('Age.txt'); age = age - mean(age); % centered age
    cov(r).c = age; cov(r).cname = 'age'; r = r+1;
    cov(r).c = readmatrix('Sex.txt'); cov(r).cname = 'sex'; r = r+1;
end


r = 0;

end

%% ------------------------------------------------------------------------
function [scans_of_roi] = create_scans_list(scans_dir, roi_prep)
% creates a list of scans as input for swe model specification "Scans" from
% directory of scans and regions of interest
if contains(roi_prep, "gsr")
    % special case if preprocessing is 'gsr'
    if endsWith(roi_prep,"_z")
        FC_files = strcat(lower(extractBefore(roi_prep,strlength(roi_prep)-4)),"cc_gsr_seed_correlation_z_trans.nii");
    else
        FC_files = strcat(lower(extractBefore(roi_prep,strlength(roi_prep)-2)),"cc_gsr_seed_correlation_trans.nii");
    end
else
    if endsWith(roi_prep,"_z")
        FC_files = strcat(lower(extractBefore(roi_prep,strlength(roi_prep)-1)),"_seed_correlation_z_trans.nii");
    else
        FC_files = strcat(lower(roi_prep),"_seed_correlation_trans.nii");
    end
end
scans_of_roi = cellstr(fullfile(scans_dir, roi_prep, FC_files));
end

%% ------------------------------------------------------------------------
function matlabbatch = RunModel(matlabbatch, location_SwE_mat)
% Estimates Model as specified in matlabbatch.
% Path to SwE.mat needs to be specified by location_SwE_mat.
% If matlabbatch defines a model for non-parametric estimation, the
% function will estimate the contrast specified in matlabbatch.
% If matlabbatch defines a model for parametric estimation, the function
% will only compute beta estimates but not contrasts.
matlabbatch{2}.spm.tools.swe.rmodel.des = cellstr(location_SwE_mat);

end


end