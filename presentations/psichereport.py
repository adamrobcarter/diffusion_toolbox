import common
import DDM.static_fourier_show
import DDM.static_fourier
import DDM.show_mult
import matplotlib.pyplot as plt
import numpy as np

path = '/home/acarter/presentations/psiche/figures'

def go_static(set):
    files = [f'psiche{i}' for i in set]
    files_small = [f'psiche{i}_small' for i in set]
    files_small_diff = [f'psiche{i}_small_diff1' for i in set]

    ## recalculate statics
    for file_small in files_small:
        DDM.static_fourier.do_static_fourier(file_small, remove_bkg=False, frame_diff=False)
    
    ##### statics
    static_fig, static_ax = plt.subplots(1, 1, figsize=(5, 5))
    
    ymin = 1e100
    ymax = 0

    used_files = files_small
    # used_files = files_small_diff

    for i, file in enumerate(used_files):
        color = f'C{i}'
        particle_diameter, k_bin_mids, v = DDM.static_fourier_show.show_static_fourier(file, static_ax, color)

        ymin = np.nanmin([ymin, v[k_bin_mids>0.1].min()])
        ymax = np.nanmax([ymax, v[k_bin_mids>0.1].max()])
    
    static_ax.set_ylim(ymin, ymax)
    # common.add_exponential_index_indicator(ax, exponent=-4, anchor=(1, 1), xlabel='k')
    static_ax.legend(fontsize=9)

    static_fig.suptitle(', '.join(files))

    filenames = '_'.join(used_files)
    filename = f'static_fourier_av_{filenames}'
    common.save_fig(static_fig, f'{path}/{filename}.pdf', hide_metadata=True)

    plt.close(static_fig) # to stop plt warnings about more than 20 figs open


def go_dynamic(set, do_fit=False, fit_use_flow=False):
    files = [f'psiche{i}' for i in set]
    dynamics_fig = DDM.show_mult.go(files, discrete_colors=True, do_fit=do_fit,
                                    fit_use_flow=fit_use_flow, log_y=True,
                                    save_data=False, legend_fontsize=7,
                                    num_displayed_ks=5, k_index_end_early=1)

    # fig.suptitle(f'{common.name(file)}, $\sigma={sigma}$, pixel$={pixel}$')
    
    dynamics_fig.suptitle(', '.join(files))

    filestring = '_'.join(files)
    filename = f'{path}/ddm_mult_{filestring}.pdf'
    common.save_fig(dynamics_fig, filename, hide_metadata=True)

    plt.close(dynamics_fig) # to stop plt warnings about more than 20 figs open


go_static(['009'])
go_static(['010'])
go_static(['020'])
go_static(['023'])
go_static(['026'])
go_static(['030'])
go_static(['034'])
go_static(['037'])
go_static(['189', '190'])
go_static(['052', '053'])
go_static(['057', '058'])
go_static(['074', '075'])
go_static(['079', '080'])
go_static(['086', '089'])
go_static(['107', '108'])
go_static(['114', '115'])
go_static(['137', '138'])
go_static(['142', '144'])
go_static(['149', '150',])
go_static(['162', '163'])
# go_static(['176', '179'])
go_static(['185', '187'])
go_static(['191', '192'])
go_static(['194', '195'])

# go_dynamic(['009'])
# go_dynamic(['010'])
# go_dynamic(['020'], do_fit=True)
# go_dynamic(['023'], do_fit=True)
# go_dynamic(['026'], do_fit=True)
# go_dynamic(['030'], do_fit=True)
# go_dynamic(['034'], do_fit=True)
# go_dynamic(['037'], do_fit=True, fit_use_flow=True)
# go_dynamic(['024'], do_fit=True, fit_use_flow=True)
# go_dynamic(['189', '190'])
# go_dynamic(['052', '053'])
# go_dynamic(['057', '058'])
# go_dynamic(['074', '075'])
# go_dynamic(['079', '080'])
# go_dynamic(['086', '089'])
# go_dynamic(['107', '108'])
# go_dynamic(['114', '115'])
# go_dynamic(['137', '138'])
# go_dynamic(['142', '144'])
# go_dynamic(['149', '150'])
# go_dynamic(['162', '163'])
# # go_dynamic(['176', '179'])
# go_dynamic(['185', '187'])
# go_dynamic(['191', '192'])
# go_dynamic(['194', '195'])