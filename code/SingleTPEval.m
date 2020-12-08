function SingleTPEval(param)
% t-test for network activity seperate for each time point.
% Normal SPM analysis - results can be viewed as is common for SPM.

% if not otherwise specified
param.tp = readmatrix(fullfile(param.INFO_DIR,'total','tp.txt'));

nrun = length(param.ROI_PREP);
ntp = length(unique(param.tp));
for j = 1:length(param.COVARIATES)
    for crun = 1:nrun % for every roi_prep
        for ctp = 1 %:ntp % for every time point
            fprintf('If necessary create output folder...\n')
            RunModel(param,j,crun,ctp)
        end
    end
end

end

%% Function ---------------------------------------------------------------
function RunModel(param,j,crun,ctp)
% Run function to evaluate models for specific time points.

param.COVARIATES = param.COVARIATES(j);
if ctp == 1
    param.wave = 'bl'; 
elseif ctp == 2
    param.wave = 'fu'; 
elseif ctp == 3 
    param.wave = 'fu2'; 
end

[out_folder, exist_already] = create_out_folder(param, crun);
if ~exist(out_folder, 'dir')
    mkdir(out_folder)
end
clear matlabbatch
fprintf('Specify model...\n')
[matlabbatch] = SpecifyModel(param, crun, ctp, out_folder);
fprintf('Run model and display results...\n')
% run SPM 
spm('defaults', 'FMRI'); 
spm_jobman('run', matlabbatch);

% [gathered_files] = MoveTmaps(param, out_folder,crun);
end 

%% ------------------------------------------------------------------------
function [matlabbatch] = SpecifyModel(param, crun, ctp, out_folder)

% Directory to get information for model definition
cd(param.INFO_DIR)

if param.ONLY_DISPLAY == false
    %% Factorial Design
    factorial_design.dir = {out_folder};

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
    % filter by current time point
    scans_of_roi(param.tp == ctp)

    %% Design matrix
    [cov, des] = ConstructDesignMatrix(param, scans_of_roi);
    factorial_design.des = des;
    factorial_design.cov = cov;

    %% Further information for estimation
    factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    factorial_design.masking.tm.tm_none = 1;
    factorial_design.masking.im = 1;

    if strcmp(param.MASK,'brain')
        mask_path = {fullfile(param.MASK_DIR,'MNI_resampled_brain_mask.nii,1')};
    elseif strcmp(param.MASK,'gm')
        mask_path = {fullfile(param.MASK_DIR,'mni_icbm152_gm_tal_nlin_sym_09a.nii,1')};
    end

    factorial_design.masking.em = mask_path;
    factorial_design.globalc.g_omit = 1;
    factorial_design.globalm.gmsca.gmsca_no = 1;
    factorial_design.globalm.glonorm = 1;
    matlabbatch{1}.spm.stats.factorial_design = factorial_design;

    %% Estimate
    fmri_est.spmmat = {fullfile(out_folder,'SPM.mat')};
    fmri_est.write_residuals = 0;
    fmri_est.method.Classical = 1;
    matlabbatch{2}.spm.stats.fmri_est = fmri_est;
    
    %% Contrast manager
    con.spmmat = {fullfile(out_folder,'SPM.mat')};
    con.consess{1}.tcon.name = 'mean';
    con.consess{1}.tcon.weights = 1;
    con.consess{1}.tcon.sessrep = 'none';
    con.delete = 0;
    matlabbatch{3}.spm.stats.con = con;

elseif param.ONLY_DISPLAY == true

    %% Results
    results.spmmat = {fullfile(out_folder,'SPM.mat')};
    results.conspec.titlestr = '';
    results.conspec.contrasts = 1;
    results.conspec.threshdesc = 'FWE';
    results.conspec.thresh = 0.05;
    results.conspec.extent = 0;
    results.conspec.conjunction = 1;
    results.conspec.mask.none = 1;
    results.units = 1;
    results.export{1}.ps = true;
    matlabbatch{1}.spm.stats.results = results; % matlabbatch{4} if done in one step with analysis

end

end

%% ------------------------------------------------------------------------
function [cov, des] = ConstructDesignMatrix(param, scans_of_roi)
    if param.COVARIATES == 41
        % t-test
        des.t1.scans = scans_of_roi;
        cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    elseif param.COVARIATES == 42
        % multiple regression
        des.mreg.scans = scans_of_roi;
            
        des.mreg.mcov(1).c = readmatrix('logmeanFD.txt');
        des.mreg.mcov(1).cname = 'logmeanFD';
        des.mreg.mcov(1).iCC = 1; % mean centering
        des.mreg.mcov(2).c = readmatrix('Age.txt');
        des.mreg.mcov(2).cname = 'age';
        des.mreg.mcov(2).iCC = 1; % mean centering
        des.mreg.mcov(3).c = readmatrix('Sex.txt');
        des.mreg.mcov(3).cname = 'sex';
        des.mreg.mcov(3).iCC = 5; % no mean centering
        des.mreg.incint = 1;
        
        %cov(1).c = readmatrix('logmeanFD.txt');
        %cov(1).cname = 'logmeanFD';
        %cov(1).iCFI = 1; % interactions none
        %cov(1).iCC = 1; % mean centering
        %cov(2).c = readmatrix('Age.txt');
        %cov(2).cname = 'age';
        %cov(2).iCFI = 1;
        %cov(2).iCC = 1;
        %cov(3).c = readmatrix('Sex.txt');
        %cov(3).cname = 'sex';
        %cov(3).iCFI = 1;
        %cov(3).iCC = 5; % no mean centering
        % overwrite cov
        cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    end
end

%% ------------------------------------------------------------------------
function [out_folder, exist_already] = create_out_folder(param, crun)
% Creates output folder and returns path to output folder as string.

if param.EXCLFD==true
    param.INFO_DIR = fullfile(param.INFO_DIR, 'ExclFD');
else
    param.INFO_DIR = fullfile(param.INFO_DIR, 'noExclFD');
end

parent_folder = 'Network_singleTP';

if param.COVARIATES == 41
    model_name = 'network';
elseif param.COVARIATES == 42
    model_name = 'network-fd-age-sex';
end

model_name = strcat(model_name,'_PE-all');

% Create folder
out_folder = fullfile(param.OUT_DIR, excl, parent_folder, param.MASK, param.ROI_PREP{crun}, param.wave, model_name);
if ~exist(out_folder, 'dir')
    exist_already = false;
    mkdir(out_folder)
else
    exist_already = true;
    return
end

end

%% ------------------------------------------------------------------------
function [gathered_files] = MoveTmaps(param, out_folder,crun)
    % extract contrast images of simple ttests and put into 
    listing = dir(out_folder);
    oldps = fullfile(out_folder,listing(endsWith({listing.name},'.ps')).name);
    newps = fullfile(param.OUT_DIR,'Network_tests',strcat(param.ROI_PREP{crun},'_',wave,'.ps'));
    movefile(oldps,newps);
    oldbeta = fullfile(out_folder,listing(startsWith({listing.name},'beta')).name);
    newbeta = fullfile(param.OUT_DIR,'Network_tests',strcat('beta_',param.ROI_PREP{crun},'_',wave,'.nii'));
    movefile(oldbeta,newbeta);
end
