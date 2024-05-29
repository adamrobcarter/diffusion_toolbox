import common
import matplotlib.pyplot as plt

for file in common.files_from_argv('intensity_correlation/data', 'correlated_'):
    data = common.load(f'intensity_correlation/data/correlated_{file}.npz')
    corr = data['I_corr']
    bins = data['bins']
    x = (bins[1:] + bins[:-1])/2

    plt.scatter(x, corr)
    plt.loglog()
    common.save_fig(plt.gcf(), f'intensity_correlation/figures_png/corr_{file}.png')