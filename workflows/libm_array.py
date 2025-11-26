import visualisation.Ds_overlapped_mult
import matplotlib.pyplot as plt
import common

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(5, 8))

visualisation.Ds_overlapped_mult.go(
    (
        dict(
            file='ld_nbody_flat_0.114_singlewall_L640_t8h_32s_dt125',
            source='f_t1024',
            msd_file='ld_nbody_flat_0.114_singlewall_L640_t1h_1s_dt125_2d',
            label = '8h'
        ),
        dict(
            file='array8_ld_nbody_flat_0.114_singlewall_L640_t1h_32s_dt125',
            source='f_t1024',
            msd_file='ld_nbody_flat_0.114_singlewall_L640_t1h_1s_dt125_2d',
            label='8 * 1h'
        ),
        dict(
            file='array28_ld_nbody_flat_0.114_singlewall_L640_t17m_32s_dt125',
            source='f_t1024',
            msd_file='ld_nbody_flat_0.114_singlewall_L640_t1h_1s_dt125_2d',
            label='28 * 17min'
        ),
        dict(
            file='ld_nbody_flat_0.114_singlewall_L640_t1h_32s_dt125',
            source='D0Sk_theory',
            label='theory, no hydrodynamics',
            color='black',
            msd_file='ld_nbody_flat_0.114_singlewall_L640_t1h_1s_dt125_2d',
        ),
    ),
    ax = ax1,
    plot_against_k=True,
    # allow_rescale_y=False,
    # labels=('hydro above')
    legend_fontsize=7,
    show_twin_k_axis=False,
    discrete_colors=True,
    allow_rescale_x=True,
)
# ax1.set_ylim(0.5, 4)
ax1.set_ylim(0.8, 5)
ax1.set_xlim(1e-2, 1e1)
ax1.semilogy()
ax1.set_title('$L=640\mathrm{\mu m}$')

visualisation.Ds_overlapped_mult.go(
    (
        dict(
            file='ld_nbody_flat_0.114_singlewall_L1280_t8h_32s_dt50',
            source='f_t1024',
            msd_file='ld_nbody_flat_0.114_singlewall_L1280_t1h_1s_dt50_2d',
            label = '8h'
        ),
        dict(
            file='array8_ld_nbody_flat_0.114_singlewall_L1280_t1h_32s_dt125',
            source='f_t1024',
            msd_file='ld_nbody_flat_0.114_singlewall_L1280_t1h_1s_dt50_2d',
            label='8 * 1h'
        ),
        dict(
            file='array27_ld_nbody_flat_0.114_singlewall_L1280_t17m_32s_dt125',
            source='f_t1024',
            msd_file='ld_nbody_flat_0.114_singlewall_L1280_t1h_1s_dt50_2d',
            label='27 * 17m'
        ),
        dict(
            file='ld_nbody_flat_0.114_singlewall_L640_t1h_32s_dt125',
            source='D0Sk_theory',
            label='theory, no hydrodynamics',
            color='black',
            msd_file='ld_nbody_flat_0.114_singlewall_L1280_t1h_1s_dt50_2d',
        ),
    ),
    ax = ax2,
    plot_against_k=True,
    # allow_rescale_y=False,
    # labels=('hydro above')
    legend_fontsize=7,
    show_twin_k_axis=False,
    discrete_colors=True,
    allow_rescale_x=True,
)
# ax2.set_ylim(0.5, 4)
ax2.set_ylim(0.8, 5)
ax2.set_xlim(1e-2, 1e1)
ax2.semilogy()
ax2.set_title('$L=1280\mathrm{\mu m}$')

common.save_fig(fig, 'workflows/figures/libm_array.png')