import common
import matplotlib.pyplot as plt
import matplotlib.cm
import numpy as np

for file in common.files_from_argv('preprocessing/data', 'stack_'):
    data = common.load(f'preprocessing/data/stack_{file}.npz')
    stack = data['stack']
    frame = stack[0, :, :] - stack[-1, :, :]
    vlim = np.abs(frame).max()

    fig, ax = plt.subplots(1, 1)
    im = ax.imshow(frame, cmap=matplotlib.cm.seismic, vmin=-vlim, vmax=vlim)
    fig.colorbar(im)
    common.save_fig(fig, f'preprocessing/figures_png/firstlast_{file}.png', dpi=300)