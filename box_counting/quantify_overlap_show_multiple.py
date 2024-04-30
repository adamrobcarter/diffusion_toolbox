import matplotlib.pyplot as plt
import common
import sys
import numpy as np

NUM_SPLITS = 8
# SPACINGS =  [100.5, 30.5, 10.5, 3.5]
# SPACINGS = [256.5, 128.5, 64.5, 32.5, 16.5, 8.5, 4.5, 2.5, 1.5, 0.5]
SPACING = 16.5

# on the alice0.02/34/66, need to check that the splits have the same time length!
# if not, crop them in time (before splitting/counting) so they do!

fig_summary, ax_summary = plt.subplots(1, 1)

std_devs = []
pack_fracs = []

for file in sys.argv[1:]:

    datas = []

    datas.append([])

    for split in range(NUM_SPLITS):
        datas[0].append(common.load(f'box_counting/data/counted_{file}_qo_split{split}_spacing{SPACING}.npz'))


    ax_i = 0
    # print(f'spacing = {spacing}')

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
        pack_frac      = data['pack_frac']

        box_sizes = N_stats[:, 0]
        N_mean    = N_stats[:, 1]
        N_var     = N_stats[:, 2]
        num_of_boxes = N_stats[:, 5]

        num_timesteps = N2_mean.shape[1]
        num_box_sizes = N2_mean.shape[0]
        t_all = np.arange(0, num_timesteps) * time_step

        for box_size_index in range(0, num_box_sizes, 4):

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
            label = rf'$L={L:.1f}\mathrm{{\mu m}}$, $n={num_of_boxes[box_size_index]}'
            # label += f', $D={D0:.3f}±{np.sqrt(pcov[0][0]):.3f}$'


    N2_means_from_diff_splits_combined = np.stack(N2_means, axis=2)
    # ^ axes are box x timestep x split

    std_devs_at_spacing = []

    for box_size_index in range(0, num_box_sizes):
        N2_at_this_box = N2_means_from_diff_splits_combined[box_size_index, :, :]
        N2_stds = N2_at_this_box.std(axis=1) # axis 1 is now split
        N2_stds /= N2_at_this_box.mean(axis=1)
        # avg_std = N2_std.mean() # only remaining axis is timestamp
        # avg_std = N2_stds[100]
        avg_std = np.median(N2_stds)
        # print(avg_std)
        print(f'L={box_sizes[box_size_index]}, {avg_std:.5f}')

        std_devs_at_spacing.append(avg_std)

    std_devs.append(std_devs_at_spacing)
    pack_fracs.append(pack_frac)
    print('adding pack frac', pack_frac)

    
std_devs  = np.array(std_devs)
# std_devs is now (num spacings) x (num box sizes)
pack_fracs = np.array(pack_fracs)
print(pack_fracs)


for box_size_index in range(num_box_sizes):
    print(pack_fracs.shape)
    print(std_devs[:, box_size_index].shape)
    label = fr'${box_sizes[box_size_index]}\mathrm{{\mu m}}$, $n={num_of_boxes[box_size_index]:.0f}$'
    ax_summary.scatter(pack_fracs, std_devs[:, box_size_index],
                        label=label)

# ax_summary.semilogx()
# ax_summary.semilogy()
ax_summary.set_xlabel('packing fraction')
ax_summary.set_ylabel('average std. dev. / mean')
ax_summary.legend()

files = '_'.join(sys.argv[1:])
common.save_fig(fig_summary, f'box_counting/figures_png/quantify_overlap_multiple_{files}.png')
    
