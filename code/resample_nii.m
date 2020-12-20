result_dir = "/data/pt_02161/Results/Project2_resting_state/connectivity/Analysis";
BMImodel_dirs =  [
    "noExclFD/BMI_total/brain/PCC_cc_z/bmi-age-sex-meanFD_WB-c01cl",% deact
    "noExclFD/BMI_total/brain/PCC_cc_z/bmi-age-sex_WB-c01cl" % deact
    ]; 
FDBMImodel_dirs = [
    "noExclFD/FD_total/brain/Nacc_cc_z/bmi-fd-age-sex_WB-c02cl", % deact
    "noExclFD/FD_total/brain/Nacc_cc_z/bmi-fd-age-sex_WB-c03cl", % act
    "noExclFD/FD_total/brain/Nacc_gsr_z/bmi-fd-age-sex_WB-c03cl", % act
    "noExclFD/FD_total/brain/PCC_cc_z/bmi-fd-age-sex_WB-c01cl" % deact
    ]; 
FDmodel_dirs = [
    "noExclFD/FD_total/brain/Nacc_cc_z/fd-age-sex_WB-c01cl", % act
    "noExclFD/FD_total/brain/Nacc_gsr_z/fd-age-sex_WB-c01cl" % act
    ];

folder = fullfile(result_dir, BMImodel_dirs);
nii = "swe_vox_zTstat_c01.nii";
V = niftiread(fullfile(folder, nii));
Vres = imresize3(V, 3, "nearest");
niftiwrite(Vres,fullfile(folder, strcat(nii{1}(1:end-4),"_res.nii")))