function SingleTPEval(param)
% t-test for network activity seperate for each time point.
% Normal SPM analysis - results can be viewed as is common for SPM.

param.tp = readmatrix(fullfile(param.INFO_DIR,'total','tp.txt'));
nrun = length(param.ROI_PREP);
ntp = length(unique(param.tp));
for crun = 1:nrun % for every roi_prep
    for ctp = 1:ntp % for every time point
        fprintf('If necessary create output folder...\n')
        if (ctp == 1) wave = 'bl'; elseif (ctp == 2) wave = 'fu'; elseif (ctp == 3) wave = 'fu2'; end
        parent_folder = 'Network_tests';
        out_folder = fullfile(param.OUT_DIR, parent_folder, param.MASK, param.ROI_PREP{crun}, wave);
        if ~exist(out_folder, 'dir')
            mkdir(out_folder)
        else
            fprintf('The analysis has been carried out before.\n')
        end
        clear matlabbatch
        fprintf('Specify model...\n')
        [matlabbatch] = SpecifyModel(param, crun, ctp, out_folder);
        fprintf('Run model and display results...\n')
        % run SPM 
        spm('defaults', 'FMRI'); 
        spm_jobman('run', matlabbatch);
        % extract picture
        listing = dir(out_folder);
        oldps = fullfile(out_folder,listing(endsWith({listing.name},'.ps')).name);
        newps = fullfile(param.OUT_DIR,'t_tests',strcat(param.ROI_PREP{crun},'_',wave,'.ps'));
        movefile(oldps,newps);
        oldbeta = fullfile(out_folder,listing(startsWith({listing.name},'beta')).name);
        newbeta = fullfile(param.OUT_DIR,'t_tests',strcat('beta_',param.ROI_PREP{crun},'_',wave,'.nii'));
        movefile(oldbeta,newbeta);
    end
end

%% Function ---------------------------------------------------------------
function [matlabbatch] = SpecifyModel(param, crun, ctp, out_folder)
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

factorial_design.des.t1.scans = scans_of_roi;

% Design matrix
factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
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
matlabbatch{4}.spm.stats.results = results;

end

end
