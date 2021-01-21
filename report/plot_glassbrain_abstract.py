# Plot DMN and Rew networks on glass brain

# use python >= 3.8

from nilearn import plotting
from os import path


def main():
    ROOT_DIR = '/data/pt_02161/Results/Project2_resting_state/connectivity/Analysis/noExclFD/FD_total/brain/Nacc_cc_z/bmi-fd-age-sex_WB-c02cl'

    IN_FILE = path.join(ROOT_DIR, 'thresholded_cgnBMI_decrease.nii')
    OUT_FILE = path.join(ROOT_DIR, 'thresholded_cgnBMI_decrease.pdf')

    print(f'Processing:\n- in:  {IN_FILE}\n- out: {OUT_FILE}')
    plot(IN_FILE, OUT_FILE)
    print('done!')


def plot(infile: str, outfile: str):
    plotting.plot_glass_brain(
        stat_map_img=infile,
        output_file=outfile,
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


if __name__ == '__main__':
    main()
