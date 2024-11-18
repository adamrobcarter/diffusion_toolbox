import matplotlib.pyplot as plt
import numpy as np
import common
import scipy.integrate, scipy.special, scipy.optimize, scipy.signal
import countoscope_theory.nmsd, countoscope_theory.structure_factor
import box_counting.D_of_L

SHOW_THEORY = False
SHOW_TIMESCALEINTEGRAL_FIT = True

LATE_CN_ALPHA = 0.2

MAX_NUM_BOXES = 10

UNC_INCLUDES_NMSD_UNC = False

NOFIT_CROP_THRESHOLD = 1e-4

RESCALE_X_L2 = 1
# RESCALE_X = RESCALE_X_L2
RESCALE_X = 0

PLATEAU_SOURCE = 'var'

# import warnings
# warnings.filterwarnings('ignore')


if __name__ == '__main__':
    fig, ax = plt.subplots(1, 1)

    files = common.files_from_argv('box_counting/data', 'counted_')
    for i, file in enumerate(files):
        color = ['tab:blue', 'tab:orange'][i]
        box_counting.D_of_L.go(file, ax=ax, plateau_source=PLATEAU_SOURCE, plot_color=color, save_data=False)

    filenames = '_'.join(files)
    common.save_fig(fig, f'box_counting/figures_png/integrand_mult_{filenames}.png')