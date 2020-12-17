#!/bin/bash


#This script is to extract agg FC values for all participants from DMN/reward network masks

output_dir="/data/pt_02161/Results/Project2_resting_state/connectivity/Analysis/aggFC/"


# DMN
# threshold: log10(p_FWE=0.05)=-1.30103

in="/data/pt_02161/Results/Project2_resting_state/connectivity/Analysis/noExclFD/Network_0-6-12/brain/PCC_gsr_z/network-age-sex_WB-c01cl/swe_clusternorm_Tstat_lpFWE-WB_c01.nii"

fslmaths $in -thr 1.30103 -bin "$output_dir/DMN/DMN_swe_clusternorm_Tstat_lpFWE0.05-WB_c01_mask.nii"

fslmerge -t "$output_dir/DMN/DMN_merged.nii.gz" `ls /data/pt_02161/Results/Project2_resting_state/connectivity/Analysis/Individual_scans/PCC_gsr_z/*`

ls /data/pt_02161/Results/Project2_resting_state/connectivity/Analysis/Individual_scans/PCC_gsr_z/* >> $output_dir/DMN/agg_FC_DMN_ID.txt

fslstats -t $output_dir/DMN/DMN_merged.nii.gz -k $output_dir/DMN/DMN_swe_clusternorm_Tstat_lpFWE0.05-WB_c01_mask.nii  -m -s >> $output_dir/DMN/agg_FC_DMN.txt

# Rew
# threshold: log10(p_FWE=0.05)=-1.30103

in="/data/pt_02161/Results/Project2_resting_state/connectivity/Analysis/noExclFD/Network_0-6-12/brain/Nacc_gsr_z/network-age-sex_WB-c01cl/swe_clusternorm_Tstat_lpFWE-WB_c01.nii"

fslmaths $in -thr 1.30103 -bin "$output_dir/Rew/Rew_swe_clusternorm_Tstat_lpFWE0.05-WB_c01_mask.nii"

fslmerge -t "$output_dir/Rew/Rew_merged.nii.gz" `ls /data/pt_02161/Results/Project2_resting_state/connectivity/Analysis/Individual_scans/Nacc_gsr_z/*`

ls /data/pt_02161/Results/Project2_resting_state/connectivity/Analysis/Individual_scans/Nacc_gsr_z/* >> $output_dir/Rew/agg_FC_Rew_ID.txt

fslstats -t $output_dir/Rew/Rew_merged.nii.gz -k $output_dir/Rew/Rew_swe_clusternorm_Tstat_lpFWE0.05-WB_c01_mask.nii  -m -s >> $output_dir/Rew/agg_FC_Rew.txt

##Also make figures of networks for manuscript

#fslmaths /data/pt_02161/Results/Project2_resting_state/connectivity/Analysis/noExclFD/Network_0-6-12/brain/Nacc_gsr_z/network-age-sex_WB-c01cl/swe_vox_Tstat_lpFWE-WB_c01.nii -mas /data/pt_02161/Results/Project2_resting_state/connectivity/Analysis/aggFC/Rew/Rew_swe_clusternorm_Tstat_lpFWE0.05-WB_c01_mask.nii.gz /data/pt_02161/Results/Project2_resting_state/connectivity/Analysis/aggFC/Rew/Rew_swe_vox_Tstat_lpFWE-WB_c01_thresholded_for_clusterFWE0.05.nii

#fslmaths /data/pt_02161/Results/Project2_resting_state/connectivity/Analysis/noExclFD/Network_0-6-12/brain/PCC_gsr_z/network-age-sex_WB-c01cl/swe_vox_Tstat_lpFWE-WB_c01.nii -mas /data/pt_02161/Results/Project2_resting_state/connectivity/Analysis/aggFC/DMN/DMN_swe_clusternorm_Tstat_lpFWE0.05-WB_c01_mask.nii.gz /data/pt_02161/Results/Project2_resting_state/connectivity/Analysis/aggFC/DMN/DMN_swe_vox_Tstat_lpFWE-WB_c01_thresholded_for_clusterFWE0.05.nii
