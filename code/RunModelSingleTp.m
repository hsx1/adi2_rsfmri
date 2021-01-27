function RunModelSingleTp(crun)
% t-test for network activity seperate for each time point.
% Normal SPM analysis - results can be viewed as is common for SPM.

if crun.EXCLFD==true
    crun.INFO_DIR = fullfile(crun.INFO_DIR, 'ExclFD');
else
    crun.INFO_DIR = fullfile(crun.INFO_DIR, 'noExclFD');
end
crun.INFO_DIR = fullfile(crun.INFO_DIR,'both/total/');
crun.tp = readmatrix(fullfile(crun.INFO_DIR,'tp.txt'));

for ctp = 1 %:ntp % for every time point
    fprintf('If necessary create output folder...\n')
    RunModel(crun,ctp)
end

end

%% Function ---------------------------------------------------------------
function RunModel(crun,ctp)
% Run function to evaluate models for specific time points.

if ctp == 1
    crun.wave = 'bl'; 
elseif ctp == 2
    crun.wave = 'fu'; 
elseif ctp == 3 
    crun.wave = 'fu2'; 
end


[out_folder, exist_already] = create_out_folder(crun);
% delete existing files in folder if existent

if exist_already
    if not(crun.OVERWRITE)&& not(crun.ONLY_DISPLAY) 
        return
    elseif crun.OVERWRITE && not(crun.ONLY_DISPLAY) 
        rmdir(out_folder, 's'); % delete former dir
        mkdir(out_folder); % create new empty one
    end
end

fprintf('Output folder: \n %s ...\n',out_folder)
clear matlabbatch
fprintf('Specify model...\n')
[matlabbatch] = SpecifyModel(crun, ctp, out_folder);
fprintf('Run model and display results...\n')
%fprintf('Evaluate network test with covariates %i...\n',crun.COVARIATES(j))
% run SPM 
spm('defaults', 'FMRI'); 
spm_jobman('run', matlabbatch);

% [gathered_files] = MoveTmaps(crun, out_folder);
end 

%% ------------------------------------------------------------------------
function [matlabbatch] = SpecifyModel(crun, ctp, out_folder)

% Directory to get information for model definition
cd(crun.INFO_DIR)

if crun.ONLY_DISPLAY == false
    %% Factorial Design
    factorial_design.dir = cellstr(out_folder);

    scans_dir = readcell(fullfile(crun.INFO_DIR,'scans.txt'), 'Delimiter',' ','Whitespace',"'");
    scans_of_roi = create_scans_list(scans_dir, crun.ROI_PREP);
    % filter by current time point (check)
    % scans_of_roi(crun.tp == ctp)

    %% Design matrix
    [cov, des] = ConstructDesignMatrix(crun, scans_of_roi);
    factorial_design.des = des;
    factorial_design.cov = cov;

    %% Further information for estimation
    factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    factorial_design.masking.tm.tm_none = 1;
    factorial_design.masking.im = 1;
    
    if strcmp(crun.MASK,'brain')
        mask_path = cellstr(fullfile(crun.MASK_DIR, crun.MASK_B));
    elseif strcmp(crun.MASK,'gm')
        mask_path = cellstr(fullfile(crun.MASK_DIR, crun.MASK_GM));
    end

    factorial_design.masking.em = mask_path;
    factorial_design.globalc.g_omit = 1;
    factorial_design.globalm.gmsca.gmsca_no = 1;
    factorial_design.globalm.glonorm = 1;
    matlabbatch{1}.spm.stats.factorial_design = factorial_design;

    %% Estimate
    fmri_est.spmmat = cellstr(fullfile(out_folder,'SPM.mat'));
    fmri_est.write_residuals = 0;
    fmri_est.method.Classical = 1;
    matlabbatch{2}.spm.stats.fmri_est = fmri_est;
    
    %% Contrast manager
	con.spmmat = cellstr(fullfile(out_folder,'SPM.mat'));
    % {'/data/pt_02161/Results/Project2_resting_state/connectivity/Analysis/noExclFD/Network_baseline/brain/Nacc_cc_z/network-fd-age-sex_PE-all/SPM.mat'};
	con.consess{1}.tcon.name = 'mean';
	con.consess{1}.tcon.weights = 1;
	con.consess{1}.tcon.sessrep = 'none';
	con.delete = 0;
	matlabbatch{3}.spm.stats.con = con;

elseif crun.ONLY_DISPLAY == true

    %% Results
    results.spmmat = cellstr(fullfile(out_folder,"SPM.mat"));
    % {'/data/pt_02161/Results/Project2_resting_state/connectivity/Analysis/noExclFD/Network_baseline/brain/PCC_cc_z/network-fd-age-sex_PE-all/SPM.mat'}; 
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
function [cov, des] = ConstructDesignMatrix(crun, scans_of_roi)
    if crun.COVARIATES == 41
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
        
        cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
        
    elseif crun.COVARIATES == 42
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
     
        cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
        
    elseif crun.COVARIATES == 43
        % simple t-test (without covariates
        des.t1.scans = scans_of_roi;
        cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    end
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
function [out_folder, exist_already] = create_out_folder(crun)
% Creates output folder and returns path to output folder as string.

if crun.EXCLFD==true
    excl = 'ExclFD';
else
    excl = 'noExclFD';
end

if strcmp(crun.wave, 'bl')
    parent_folder = 'Network_baseline';
elseif strcmp(crun.wave, 'fu')
    parent_folder = 'Network_6months';
elseif strcmp(crun.wave, 'fu2')
    parent_folder = 'Network_12months';
end


if crun.COVARIATES == 41
    model_name = 'network-fd-age-sex';
elseif crun.COVARIATES == 42
    model_name = 'network-age-sex';
elseif crun.COVARIATES == 43
    model_name = 'network';
end

model_name = strcat(model_name,'_PE-all');

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

