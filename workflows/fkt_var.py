import common
import matplotlib.pyplot as plt
import numpy as np
import isf.calc_f_first
import MSD.calc
import tqdm

sets = [
    # 'sim_nohydro_0.114_L320_seed*_t1m_0.5',
    # 'sim_nohydro_0.114_L320_seed*_t4m_0.5',
    # 'sim_nohydro_0.114_L320_seed*_t15m_0.5',
    'sim_nohydro_0.114_L320_seed*_t1h_0.5',
    # 'sim_nohydro_0.114_L320_seed*_t4h_0.5',
    # 'sim_nohydro_0.114_L320_seed*_t16h_0.5',
    # 'sim_nohydro_0.114_L1280_seed*_t15m_0.5',
    'sim_nohydro_0.114_L1280_seed*_t1h_0.5',
]

# died at particle_linking/data/trajs_sim_nohydro_0.114_L1280_seed10_t1h_0.5_unwrap_div8.npz

## do the calculation
for set in sets:
    files = common.files_from_filenames('particle_detection/data', 'particles_', set)

    # for file in tqdm.tqdm(files, desc='calcing f(k, t)'):
    #     isf.calc_f_first.go(file, quiet=True)

    # for file in tqdm.tqdm(files, desc='showing f(k, t)'):
    #     from isf.show_both import show_single_F_type
    #     # ^^^ rn we do this hack cause the figure is created at import inside show_both
    #     num_displayed_ks = 2
    #     fig, axes = plt.subplots(4, num_displayed_ks)
    #                                     # hack      ^^
    #     show_single_F_type(
    #         file,
    #         type_index=0,
    #         Ftype='f_first',
    #         # Ftype='F_first32',
    #         fig=fig,
    #         axes=axes,
    #         num_displayed_ks=num_displayed_ks,
    #         quiet=True,
    #     )

    # for file in tqdm.tqdm(files, desc='calcing MSDs'):
    #     # if 'seed8' not in file:
    #     #     continue
    #     msd_file = f'{file}_unwrap_div8'
    #     if 't4m' in file or 't1m' in file:
    #         msd_file = f'{file}_unwrap'
    #     MSD.calc.go(msd_file, quiet=True)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(5, 7))

for set_i, set in enumerate(sets):
    files = common.files_from_filenames('visualisation/data', 'Ds_from_f_first_first_', set)
    print(f'found {len(files)} files in {set}')

    all_Ds = []
    all_ks = None
    for file in files:
        data = common.load(f'visualisation/data/Ds_from_f_first_first_{file}.npz')
        Ds = data['Ds']
        ks = data['ks']
        all_Ds.append(Ds)
        if all_ks is None:
            all_ks = ks
        else:
            assert np.all(ks == all_ks)

    variance = np.nanvar (all_Ds, axis=0)
    mean     = np.nanmean(all_Ds, axis=0)

    ax1.scatter(all_ks, variance/mean, label=set, color=common.colormap(set_i, 0, len(sets)))

    # ax2.errorbar(all_ks, mean, yerr=np.sqrt(variance), label=set, color=common.colormap(set_i, 0, len(sets)), linestyle='none', marker='o')
    ax2.scatter(all_ks, mean, label=set, color=common.colormap(set_i, 0, len(sets)))
    ax2.set_ylim(np.nanquantile(mean, 0.03), np.nanquantile(mean, 0.97))

    # all_Ds = np.array(all_Ds)
    # print(all_Ds.shape)

ax1.legend()
ax1.set_xlabel('$k$')
ax1.set_ylabel('$\mathrm{Var}(D)/D$')
ax1.semilogx()
ax1.semilogy()
ax1.grid()

ax2.legend()
ax2.set_xlabel('$k$')
ax2.set_ylabel('$D$')
ax2.semilogx()
common.save_fig(fig, 'workflows/figures/fkt_var.png')