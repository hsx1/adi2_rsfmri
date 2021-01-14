# Plot DMN and Rew networks on glass brain

# use python 3 

from nilearn import plotting
stat_img="/data/pt_02161/Results/Project2_resting_state/connectivity/Analysis/aggFC/DMN/DMN_swe_vox_Tstat_lpFWE-WB_c01_thresholded_for_clusterFWE0.05.nii.gz"
plotting.plot_glass_brain(stat_img, output_file="/data/pt_02161/Results/Project2_resting_state/connectivity/Analysis/aggFC/DMN/DMN_thresholded_for_clusterFWE0.05.png", display_mode='ortho', colorbar=True, figure=None, axes=None, title=None, threshold='auto', annotate=True, black_bg=False, cmap='plasma', alpha=0.7, vmin=None, vmax=None, plot_abs=True, symmetric_cbar='auto', resampling_interpolation='continuous')

stat_img="/data/pt_02161/Results/Project2_resting_state/connectivity/Analysis/aggFC/Rew/Rew_swe_vox_Tstat_lpFWE-WB_c01_thresholded_for_clusterFWE0.05.nii.gz"
plotting.plot_glass_brain(stat_img, output_file="/data/pt_02161/Results/Project2_resting_state/connectivity/Analysis/aggFC/Rew/Rew_thresholded_for_clusterFWE0.05.png", display_mode='ortho', colorbar=True, figure=None, axes=None, title=None, threshold='auto', annotate=True, black_bg=False, cmap='plasma', alpha=0.7, vmin=None, vmax=None, plot_abs=True, symmetric_cbar='auto', resampling_interpolation='continuous')



