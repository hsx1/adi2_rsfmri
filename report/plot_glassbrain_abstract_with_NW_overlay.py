# Plot DMN and Rew networks on glass brain

# use python >= 3.8

from nilearn import plotting
from os import path
import pandas as pd


def main():

    OUT_DIR = "/data/pt_02161/Analysis/Project2_resting_state/seed-based/Second_level /code_and_manuscript/report"
    INFO_FILE = "/data/pt_02161/Results/Project2_resting_state/connectivity/Analysis/noExclFD/result_report.txt"

    # dataframe with filenames, path list, cluster threshold
    df = pd.read_table(INFO_FILE, header=None)
    df[1] = df[1].str.replace('txt', 'nii') # root_paths  = df[1].tolist()

    
	
    for i in range(len(df[1])):
        IN_FILE = df[1][i]
        #NW = "/data/pt_02161/Results/Project2_resting_state/connectivity/Analysis/aggFC/DMN/DMN_swe_vox_Tstat_lpFWE-WB_c01_thresholded_for_clusterFWE0.05.nii.gz"
        if ('PCC' in df[0][i]):
         NW="/data/pt_02161/Results/Project2_resting_state/connectivity/Analysis/aggFC/DMN/DMN_swe_vox_Tstat_lpFWE-WB_c01_thresholded_for_clusterFWE0.05_bin.nii.gz"
        else:
         NW="/data/pt_02161/Results/Project2_resting_state/connectivity/Analysis/aggFC/Rew/Rew_swe_vox_Tstat_lpFWE-WB_c01_thresholded_for_clusterFWE0.05_bin.nii.gz"
        OUT_FILE = path.join(OUT_DIR, (df[0][i] + "_overlay.pdf"))

        # do for each result
        #print(f'Processing:\n- in:  {IN_FILE}\n- out: {OUT_FILE}')
        plot(IN_FILE, NW, OUT_FILE)

    print('done!')


def plot(infile: str, NW:str, outfile: str):
    display=plotting.plot_glass_brain(
        stat_map_img=infile,
        display_mode='ortho',
        colorbar=True,
        figure=None,
        axes=None,
        title=None,
        threshold='auto',
        annotate=True,
        black_bg=False,
        cmap='plasma',
        alpha=0.7,
        vmin=None,
        vmax=None,
        plot_abs=True,
        symmetric_cbar='auto',
        resampling_interpolation='continuous')
    display.add_contours(NW, linewidths=.3, alpha=0.9, colors='r')
    display.savefig(outfile)



if __name__ == '__main__':
    main()
