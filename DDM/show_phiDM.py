import common
import matplotlib.pyplot as plt

if __name__ == '__main__':
    for file in common.files_from_argv('DDM/data/', 'phiDM_'):
        data = common.load(f'DDM/data/phiDM_{file}.npz')
        R = data['R']

        R = R.cumsum(axis=0)

        fig, (ax_x, ax_y) = plt.subplots(2, 1)

        print(R)

        ax_x.plot(R[:, 0])
        ax_y.plot(R[:, 1])

        common.save_fig(fig, f'DDM/figures_png/phiDM_{file}.png')