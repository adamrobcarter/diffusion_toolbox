import common
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize
import scipy.special
import DDM.show

DISCRETE_COLORS = True

def go(files, discrete_colors, num_displayed_ks=14, do_fit=False, k_index_offset=0, **kwargs):

    fig, axs = plt.subplots(1, num_displayed_ks, figsize=(num_displayed_ks*3.5, 3.7))
    [ax.semilogy() for ax in axs]

    all_ymins = np.full(num_displayed_ks, np.nan)
    all_ymaxs = np.full(num_displayed_ks, np.nan)

    for file_index, file in enumerate(files):
        # we load this stuff here, not in `show()` as we normally do, so that we can use `show()` in the live situation
        data = common.load(f'DDM/data/ddm_{file}.npz')
        k          = data['k']
        F_D_sq     = data['F_D_sq']
        F_D_sq_unc = data['F_D_sq_unc']
        t          = data['t']

        sigma = data.get('particle_diameter')
        pixel = data['pixel_size']
        NAME  = data.get('NAME')
        channel = data.get('channel')

        k_index_offset = 0
        if 'psiche' in file:
            k_index_offset = k.size // 2

        if discrete_colors:
            color = f'C{file_index}'
        else:
            color = common.colormap(file_index, 0, len(files))

        label = data.get('NAME', file)
    
        ymins, ymaxs, tmin, tmax = DDM.show.show(file, axs, k, F_D_sq, F_D_sq_unc, t, sigma, pixel, NAME=NAME, channel=channel,
                num_displayed_ks=num_displayed_ks, k_index_offset=k_index_offset, particle_diameter=data.get('particle_diameter'),
                color=color, do_fit=do_fit, plot_label=label, particle_material=data.get('particle_material'), **kwargs)
        
    #     assert np.isfinite(ymins).all()
    #     assert np.isfinite(ymaxs).all()
    #     all_ymins = np.nanmin([all_ymins, ymins], axis=0)
    #     all_ymaxs = np.nanmax([all_ymaxs, ymaxs], axis=0)


    # for i, ax in enumerate(axs):
    #     ax.set_ylim(all_ymins[i]*0.95, all_ymaxs[i]*1.05)
        
                
    return fig

if __name__ == '__main__':

    files = common.files_from_argv('DDM/data', 'ddm_')

    num_displayed_ks = 14
    if 'psiche' in files[0]:
        num_displayed_ks = 5
        
    fig = go(files, discrete_colors=DISCRETE_COLORS, num_displayed_ks=num_displayed_ks)

    # fig.suptitle(f'{common.name(file)}, $\sigma={sigma}$, pixel$={pixel}$')

    filestring = '_'.join(files)
    filename = f'DDM/figures_png/ddm_mult_{filestring}.png'
    common.save_fig(fig, filename)