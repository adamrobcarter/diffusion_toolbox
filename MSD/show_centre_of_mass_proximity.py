import matplotlib.pyplot as plt
import common
import numpy as np
import scipy.optimize
import MSD.show

if __name__ == '__main__':
    for file in common.files_from_argv('MSD/data', 'msd_centre_of_mass_proximity_'):
        fig, ax = plt.subplots(1, 1)
        data = common.load(f'MSD/data/msd_centre_of_mass_proximity_{file}.npz')
        t_orig = np.arange(0, data['msds'].shape[1]) * data['time_step']

        num_used_particles = data['num_used_particles']
        density            = data['density']
        particle_diameter  = data['particle_diameter']
        msds               = data['msds']
        msd_uncs           = data['msd_uncs']

        Ds = np.full_like(num_used_particles, np.nan)
        D_uncs = np.full_like(num_used_particles, np.nan)

        for box_index, box_size in enumerate(groupsizes := data['box_sizes']):
            # if box_index % 5 != 0:
            #     continue

            msd               = msds[box_index, :]
            msd_unc           = msd_uncs[box_index, :]

            isnan = np.isnan(msd)
            msd     = msd    [~isnan]
            msd_unc = msd_unc[~isnan]
            t       = t_orig [~isnan]

            if msd.size == 0:
                print(f'skipping L={box_size:.3g}')
                continue

            # print(f'predicted N {density*box_size**2:.3g} used {num_used_particles[box_index]:.3g}')

            color = common.colormap(box_index, 0, len(groupsizes))

            normalisation = density*box_size**2
            normalisation = num_used_particles[box_index]

            real_L = np.sqrt(num_used_particles[box_index] / density)

            print(f'N={num_used_particles[box_index]:.3g} is L={real_L/particle_diameter:.3g}σ')

            label = fr'$L \approx {real_L/particle_diameter:.1f} \sigma$'
            ax.errorbar(t, msd*normalisation, yerr=msd_unc*normalisation, marker='.', markersize=5, linestyle='none', color=color, label=label)

            fits = MSD.show.fit_msd(t, msd*normalisation, msd_unc*normalisation)
            ax.plot(fits['short']['t'], fits['short']['MSD'], color='grey')
            Ds    [box_index] = fits['short']['D']
            D_uncs[box_index] = fits['short']['D_unc']

        data_msd = common.load(f'MSD/data/msd_{file}.npz')
        msd     = data_msd['msd']
        t_indexes = np.arange(0, msd.size)
        t = np.arange(0, msd.size) * data_msd['time_step']

        # t   = t  [msd > 1e-4]
        # msd = msd[msd > 1e-4]

        color='tab:blue'
        ax.plot(t[1:], msd[1:], marker='.', markersize=8, linestyle='none', color=color, label='regular MSD')
        

        ax.loglog()

        ax.set_ylabel(r'$L^2\rho\langle r(t)^2 \rangle$ ($\mathrm{\mu m}$)')
        ax.set_xlabel('$t$ (s)')
        ax.legend()

        ax.set_ylim(0.8e-1, 1e3)

        nans = np.isnan(Ds)
        Ds = Ds[~nans]
        D_uncs = D_uncs[~nans]
        num_used_particles = num_used_particles[~nans]

        # common.save_fig(fig, f'/home/acarter/presentations/intcha24/figures/msd_{file}.pdf', hide_metadata=True)
        common.save_fig(fig, f'MSD/figures_png/msd_centre_of_mass_proximity_{file}.png')
        common.save_data(f'visualisation/data/Ds_from_MSD_centre_of_mass_proximity_{file}',
                Ds=Ds, D_uncs=D_uncs,
                Ns=num_used_particles, density=density,
                particle_diameter=data.get('particle_diameter'),
                window_size_x=data['window_size_x'], window_size_y=data['window_size_y']
        )