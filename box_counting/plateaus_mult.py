import matplotlib.pyplot as plt
import numpy as np
import common
import countoscope_theory.nmsd, countoscope_theory.structure_factor
from .msd_single import get_plateau, PLATEAU_SOURCES, PLATEAU_SOURCE_NAMES
import scipy.optimize
import visualisation.Ds_overlapped
import tqdm

from box_counting.plateaus import go


if __name__ == '__main__':
    
    fig, ax = plt.subplots(1, 1)

    cols = ['tab:orange', 'tab:green', 'tab:blue']

    for file_i, file in enumerate(files := common.files_from_argv('box_counting/data/', 'counted_')):

        go(file, ax, ['var', 'sDFT'], rescale_window=False, label_prefix=f'{file} ', rescale_sigma=False,
           colors=[cols[file_i], cols[file_i]])
        # go(file, ax, PLATEAU_SOURCES)

    files = '_'.join(files)
    common.save_fig(fig, f'box_counting/figures_png/plateaus_mult_{files}.png', dpi=200)