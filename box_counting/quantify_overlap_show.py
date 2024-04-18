import matplotlib.pyplot as plt
import common
import sys
import numpy as np

NUM_SPLITS = 10
SPACINGS = [100, 30, 10, 3]

for file in sys.argv[1:]:

    datas = []

    for ax_i, spacing in enumerate(SPACINGS):

        datas.append([])

        for split in range(NUM_SPLITS):
            datas[ax_i].append(common.load(f'box_counting/data/counted_{file}_qo_split{split}_spacing{spacing}.npz'))
            

    fig, axs = plt.subplots(1, 2, figsize=(6, 3))

    for ax_i, spacing in enumerate(SPACINGS):
        print(f'spacing = {spacing}')
        ax = axs[ax_i]

        N2_means = []

        for split_i, split in enumerate(range(NUM_SPLITS)):
            # data = common.load(f'box_counting/data/counted_{file}_qo_split{split}_spacing{spacing}.npz')
            data = datas[ax_i][split_i]
            N2_mean        = data['N2_mean']
            N2_means.append(N2_mean)
            N2_std         = data['N2_std']
            N_stats        = data['N_stats']
            phi            = data['pack_frac']
            sigma          = data['particle_diameter']
            time_step      = data['time_step']

            box_sizes = N_stats[:, 0]
            N_mean    = N_stats[:, 1]
            N_var     = N_stats[:, 2]
            num_of_boxes = N_stats[:, 5]
            ax.set_title(f'{np.mean(num_of_boxes):.0f} boxes')

            num_timesteps = N2_mean.shape[1]
            num_boxes     = N2_mean.shape[0]
            t_all = np.arange(0, num_timesteps) * time_step

            for box_size_index in range(0, num_boxes, 4):

                L = box_sizes[box_size_index]

                delta_N_sq     = N2_mean[box_size_index, :]
                delta_N_sq_err = N2_std [box_size_index, :]
                t = np.copy(t_all)
                t_theory = np.logspace(np.log10(t_all[1] / 2), np.log10(t_all.max()*2), 100)

                anomalous = delta_N_sq < 1e-14
                anomalous[0] = False # don't want to remove point t=0 as it could legit be zero
                if np.any(anomalous):
                    # print(f'found {anomalous.sum()/delta_N_sq.size*100:.3f}% anomalous')
                    delta_N_sq     = delta_N_sq    [~anomalous]
                    delta_N_sq_err = delta_N_sq_err[~anomalous]
                    t              = t             [~anomalous]
                
                # N2_fu
                label = rf'$L={L:.1f}\mathrm{{\mu m}}$'
                # label += f', $D={D0:.3f}±{np.sqrt(pcov[0][0]):.3f}$'


                ax.plot(t[1:], delta_N_sq[1:], linestyle='none', marker='o', markersize=5)

        N2_means_from_diff_splits_combined = np.stack(N2_means, axis=2)
        # ^ axes are box x timestep x split

        for box_size_index in range(0, num_boxes):
            N2_at_this_box = N2_means_from_diff_splits_combined[box_size_index, :, :]
            N2_stds = N2_at_this_box.std(axis=1) # axis 1 is now split
            # avg_std = N2_std.mean() # only remaining axis is timestamp
            # avg_std = N2_stds[100]
            avg_std = np.median(N2_stds)
            # print(avg_std)
            print(f'L={box_sizes[box_size_index]}, {avg_std:.5f}')

        ax.loglog()
    common.save_fig(fig, f'box_counting/figures_png/quantify_overlap_{file}.png')