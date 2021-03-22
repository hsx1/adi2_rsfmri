# Plot DMN and Rew networks on glass brain

# use python 3

from nilearn import plotting
stat_img="/data/pt_02161/Results/Project2_resting_state/connectivity/Analysis/noExclFD/Network_0-6-12/brain/PCC_gsr_z/network-age-sex_WB-c01cl/thresholded.nii"
plotting.plot_glass_brain(stat_img, output_file="/data/pt_02161/Analysis/Project2_resting_state/seed-based/Second_level /code_and_manuscript/report/fig/DMN_thresholded_for_clusterFWE0.05.pdf", display_mode='ortho', colorbar=True, figure=None, axes=None, title=None, threshold='auto', annotate=True, black_bg=False, cmap='plasma', alpha=0.7, vmin=None, vmax=None, plot_abs=True, symmetric_cbar='auto', resampling_interpolation='continuous')

stat_img="/data/pt_02161/Results/Project2_resting_state/connectivity/Analysis/noExclFD/Network_0-6-12/brain/Nacc_gsr_z/network-age-sex_WB-c01cl/thresholded.nii"
plotting.plot_glass_brain(stat_img, output_file="/data/pt_02161/Analysis/Project2_resting_state/seed-based/Second_level /code_and_manuscript/report/fig/Rew_thresholded_for_clusterFWE0.05.pdf", display_mode='ortho', colorbar=True, figure=None, axes=None, title=None, threshold='auto', annotate=True, black_bg=False, cmap='plasma', alpha=0.7, vmin=None, vmax=None, plot_abs=True, symmetric_cbar='auto', resampling_interpolation='continuous')
