import common
import matplotlib.pyplot as plt

METHOD_NONE = 0
METHOD_DIFF = 1

# METHOD = METHOD_DIFF
METHOD = METHOD_NONE

CUMULATIVE_MEAN = False

for file in common.files_from_argv('preprocessing/data/', 'stack_'):
    data = common.load(f'preprocessing/data/stack_{file}.npz')
    stack = data['stack']
    intensities = stack.mean(axis=(1, 2))

    if METHOD == METHOD_DIFF:
        intensities = intensities[1:] - intensities[:-1]

    print('total average intensity', intensities.mean())

    if CUMULATIVE_MEAN:
        for i in range(intensities.shape[0]-1, 0, -1):
            intensities[i] = intensities[:i].mean()

    fig, ax = plt.subplots(1, 1)
    ax.plot(intensities)
    ax.set_xlabel('frame')
    ax.set_ylabel('average intensity')
    # ax.set_ylim(0, intensities.max()*1.05)
    if not CUMULATIVE_MEAN and METHOD == METHOD_DIFF:
        ax.set_title('intensity of frame differences')

    filename = f'intensities_vs_time'
    if CUMULATIVE_MEAN:
        filename += 'cumulmean'
    filename += f'_{file}'
    common.save_fig(fig, f'preprocessing/figures_png/{filename}.png')