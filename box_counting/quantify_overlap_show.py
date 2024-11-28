import matplotlib.pyplot as plt
import common
import sys
import numpy as np
import matplotlib.cm

# SPACINGS = np.round(np.logspace(np.log10(2.5), np.log10(156.5), 50), decimals=2) # extra4

# on the alice0.02/34/66, need to check that the splits have the same time length!
# if not, crop them in time (before splitting/counting) so they do!

PERIMETER = 0
NUM_BOXES = 1
SEP = 2
TOTAL_AREA = 3
X_AXIS = NUM_BOXES

SHOW_SLOPE = False

def go(file, ax):
    data = common.load(f'box_counting/data/quantify_overlap_{file}.npz')
    spacings   = data['SPACINGS']
    num_splits = data['NUM_SPLITS']
    sigma      = data['spacing0_split0_particle_diameter']

    # fig, axs = plt.subplots(1, 2, figsize=(6, 3))

    std_devs = []
    num_boxes = []
    seps = []
    # box_sizes

    ax_i = 0

    for spacing_i, spacing in enumerate(spacings):


        # print(f'spacing = {spacing}')

        N2_means = []

        # plot = spacing_i in [2, 7]
        # # print('#################', plot, ax_i, spacing_i)
        # ax = axs[ax_i] if plot else None
        # ax_i += 1 if plot else 0

        for split_i, split in enumerate(range(num_splits)):
            # data = common.load(f'box_counting/data/counted_{file}_qo_split{split}_spacing{spacing}.npz')
            # data = data[spacing_i][split_i]
            key = f'spacing{spacing_i}_split{split_i}'
            N2_mean        = data[f'{key}_N2_mean']
            N2_means.append(N2_mean)
            N2_std         = data[f'{key}_N2_std']
            phi            = data[f'{key}_pack_frac']
            sigma          = data[f'{key}_particle_diameter']
            time_step      = data[f'{key}_time_step']
            box_sizes      = data[f'{key}_box_sizes']
            sep_sizes      = data[f'{key}_sep_sizes']
            num_of_boxes   = data[f'{key}_num_boxes']
            
            # N_mean    = N_stats[:, 1]
            # N_var     = N_stats[:, 2]
            # num_of_boxes = data[f'{key}_num_boxes']
            # ax.set_title(f'{np.mean(num_of_boxes):.0f} boxes') if plot else None

            num_timesteps = N2_mean.shape[1]
            num_box_sizes = N2_mean.shape[0]
            t_all = np.arange(0, num_timesteps) * time_step

            for box_size_index in range(0, num_box_sizes, 1):

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

                color =  matplotlib.cm.afmhot(np.interp(box_size_index, (0, len(box_sizes)), (0.2, 0.75)))
                # color = matplotlib.cm.tab10(box_size_index)
                # ax.plot(t[1:], delta_N_sq[1:], linestyle='none', marker='o', markersize=1, color=color) if plot else None

        # [print(a.shape) for a in N2_means]
        N2_means = [a[:, :7198] for a in N2_means]
        N2_means_from_diff_splits_combined = np.stack(N2_means, axis=2)
        # ^ axes are box x timestep x split

        std_devs_at_spacing = []

        for box_size_index in range(0, num_box_sizes):
            N2_at_this_box = N2_means_from_diff_splits_combined[box_size_index, :, :]
            # print('N2_at_this_box', common.nanfrac(N2_at_this_box))
            # N2_stds = N2_at_this_box.std(axis=1) # axis 1 is now split
            # N2_stds /= N2_at_this_box.mean(axis=1)
            N2_stds =  np.nanstd (N2_at_this_box, axis=1) # axis 1 is now split
            N2_stds /= np.nanmean(N2_at_this_box, axis=1)
            N2_stds /= np.sqrt(num_splits)
            # avg_std = N2_std.mean() # only remaining axis is timestamp
            # avg_std = N2_stds[100]
            # print('N2_stds', common.nanfrac(N2_stds))
            avg_std = np.median(N2_stds)
            # print(avg_std)
            # print(f'L={box_sizes[box_size_index]}, {avg_std:.5f}')

            std_devs_at_spacing.append(avg_std)

        std_devs.append(std_devs_at_spacing)
        num_boxes.append(num_of_boxes)
        seps.append(sep_sizes)

        # ax.loglog() if plot else None
    # common.save_fig(fig, f'box_counting/figures_png/quantify_overlap_{file}.png', dpi=200)

    
    std_devs  = np.array(std_devs)
    # std_devs is now (num spacings) x (num box sizes)
    num_boxes = np.array(num_boxes)
    seps = np.array(seps)

    total_area = num_boxes * box_sizes**2
    # area_fraction = total_area / ()
    total_perim = num_boxes * box_sizes * 4

    
    for box_size_index in range(num_box_sizes):
        # print(num_boxes.shape)
        # print(std_devs[:, box_size_index].shape)
        # color = matplotlib.cm.afmhot((box_size_index+1)/(len(box_sizes)+5))
        color =  matplotlib.cm.afmhot(np.interp(box_size_index, (0, len(box_sizes)), (0.15, 0.75)))
        # x = num_boxes[:, box_size_index]
        # x = total_area[:, box_size_index]
        if X_AXIS == PERIMETER:
            x = total_perim[:, box_size_index]
        elif X_AXIS == NUM_BOXES:
            x = num_boxes[:, box_size_index]
        elif X_AXIS == SEP:
            x = seps[:, box_size_index]
        elif X_AXIS == TOTAL_AREA:
            x = total_area[:, box_size_index]
        else:
            raise
        # ax_summary.scatter(x, std_devs[:, box_size_index],
        #                    label=fr'${box_sizes[box_size_index]}\mathrm{{\mu m}}$',
        #                    color=color)
        ax.plot(x, std_devs[:, box_size_index],
                        marker='o',
                        label=fr'$L={box_sizes[box_size_index]/sigma:.2g}σ$',
                        color=color)

    
    if X_AXIS == PERIMETER:
        ax.set_xlabel('total perimeter ($\mathrm{\mu m}$)')
    elif X_AXIS == NUM_BOXES:
        ax.set_xlabel('number of boxes')
    elif X_AXIS == SEP:
        ax.set_xlabel('sep ($\mathrm{\mu m}$)')
    elif X_AXIS == TOTAL_AREA:
        total_area[:, box_size_index]
        ax.set_xlabel('total area ($\mathrm{\mu m^2}$)')
    ax.legend()
    ax.set_ylabel(r'error on $\langle \Delta N^2(t) \rangle$')
    
    # common.save_fig(fig_summary, f'box_counting/figures_png/quantify_overlap_summary_{file}.png')

    ax.semilogy()
    # func = lambda x, m, c: m*x + c
    # log_func = lambda x, m, c: np.log10(func(x, m, c))
    # popt, pcov = scipy.optimize.curve_fit(log_func, np.log10())
    if X_AXIS is not SEP:
        ax.semilogx()

        if SHOW_SLOPE:
            x_ind = np.logspace(0.7, 2.5)
            m = -1/2
            c = -0.5
            y_ind = 10**(m * np.log10(x_ind) + c)
            ax.plot(x_ind, y_ind, color='black', linewidth=2)
            print(y_ind)
            ax.text(0.3e2, 2e-2, r'$n^{-\frac{1}{2}}$', fontsize=14)



if __name__ == '__main__':
        
    for file in sys.argv[1:]:
        fig, ax = plt.subplots(1, 1)

        go(file, ax)
        
        common.save_fig(fig, f'box_counting/figures_png/quantify_overlap_summary_{file}_logy.png')
    