import visualisation.Ds_overlapped_mult
import matplotlib.pyplot as plt
import common
import numpy as np

fig, ax = plt.subplots(figsize=(5, 4))

data = common.load('visualisation/data/Ds_from_boxcounting_first_quad_eleanorsmall016short')
D0_firstquad = data['Ds'][15]
print('DDs')
print(data['Ds'])
# raise Exception(f'D0 from first quad: {D0_firstquad}')

phi = 0.01
sigma = 1.395*2
z_mult = 1.06

# rescale_y = True
rescale_y = False

z_mult_str = f'{z_mult:.2f}a'
a = sigma/2
z = z_mult * a

hydrofile_wall          = '/data2/acarter/direct_to_H/' f'data/H_of_k_L300_numkbins60_single_wall_walllub_phi0.181_maxk10.0_montecarlo_full_zpos{z_mult_str}_zwidth0.0a_numis200.npz'
hydrofile_wall_particle = '/data2/acarter/direct_to_H/' f'data/H_of_k_L500_numkbins60_single_wall_walllub_particlelub_phi0.181_maxk10.0_montecarlo_full_zpos{z_mult_str}_zwidth0.0a_numis200.npz'

shortfile = 'ld_nbody_flat_0.01_singlewall_L2560_t1h_1s_dt125'
longfile  = 'ld_nbody_flat_0.01_singlewall_L2560_t8h_32s_dt125'

visualisation.Ds_overlapped_mult.go(
    (
        dict(
            file=shortfile,
            source='timescaleint_nofit_cropped_var',
            plot_index=np.index_exp[:16],
            msd='eleanorsmall016short',
            # msd=D0_firstquad,
            label='experiment',
        ),
        dict(
            file=shortfile,
            source='MSD_first',
            msd=shortfile,
            # msd=D0_firstquad,
            label='MSD',
            linestyle='-'
        ),
        dict(
            file=longfile,
            source='timescaleint_nofit_cropped_var',
            plot_index=np.index_exp[16:-5],
            msd=shortfile,
            label=None,
            color='tab:blue'
        ),
        dict(
            file=shortfile,
            source='D_of_L_theory',
            # msd='eleanorsmall016short',
            msd = common.stokes_einstein_D(sigma),
            label = 'theory, no hydro',
            phi = phi,
            sigma = sigma,
        ),
        # dict(
        #     file='eleanorsmall016short',
        #     source='D_of_L_theory_hydro',
        #     # plot_index=np.index_exp[30:-30],
        #     # msd_file='eleanorsmall016short',
        #     hydro_source_file=hydrofile_wall,
        #     label = f'theory, hydro + wall lubrication z={z_mult_str}',
        #     phi = phi,
        #     sigma = sigma,
        #     Lmax = 500,
        # ),
        # dict(
        #     file='eleanorsmall016short',
        #     source='D_of_L_theory_hydro',
        #     # plot_index=np.index_exp[30:-30],
        #     # msd_file='eleanorsmall016short',
        #     hydro_source_file=hydrofile_wall_particle,
        #     label = f'theory, hydro + wall + particle lubrication $z={z_mult_str}$',
        #     phi = phi,
        #     sigma = sigma,
        #     Lmax = 500,
        # ),
    ),
    ax = ax,
    plot_against_k=False,
    discrete_colors=True,
    allow_rescale_y=rescale_y,
    legend_fontsize=8,
    allow_rescale_x=True,
)
ax.semilogy()

if not rescale_y:
    z = a*z_mult
    import rigiddynamics
    D_SE = common.stokes_einstein_D(sigma)
    ax.axhline(D_SE, color='black', linestyle=':', label='stokes-einstein D')
    D_SE_wall = D_SE * rigiddynamics.theory.Dt_correction_para(a/z)
    ax.axhline(D_SE_wall, color='gray', linestyle=':', label=f'stokes-einstein D, wall at z={z/a}a')
    D_SE_wall_particle = D_SE_wall * (1 - 0.85 * phi)
    ax.axhline(D_SE_wall_particle, color='gray', linestyle='--', label=f'stokes-einstein D , wall at z={z/a}a, particle lub')

ax.legend(fontsize=7)
ax.set_xlim(1e-1, 1e2)

common.save_fig(fig, 'workflows/figures/donev_phi0p01.pdf', hide_metadata=True)
common.save_fig(fig, 'workflows/figures/donev_phi0p01.png')

# python run.py --wall=single_wall_walllub_particlelub --phi=0.181 --method=montecarlo_full --zoffset=2 --a=1.005 --L=200 --num_is=100