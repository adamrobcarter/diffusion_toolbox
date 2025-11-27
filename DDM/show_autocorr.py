import common
import matplotlib.pyplot as plt
import numpy as np

if __name__ == '__main__':
    for file in common.files_from_argv('DDM/data', 'autocorr_'):
        data = common.load(f'DDM/data/autocorr_{file}.npz')
        
        autocorr_avg = data['autocorr_avg']
        autocorr_std = data['autocorr_std']
        t             = data['t']

        fig, ax = plt.subplots(1, 1)

        plot = t >= 0

        ax.errorbar(t[plot], autocorr_avg[plot], yerr=autocorr_std[plot]/np.sqrt(2048**2))

        common.save_fig(fig, f'DDM/figures_png/autocorr_{file}.png')