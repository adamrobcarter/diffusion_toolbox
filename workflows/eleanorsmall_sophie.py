import visualisation.Ds_overlapped_mult
import matplotlib.pyplot as plt
import common
import numpy as np

fig, ax = plt.subplots(figsize=(4, 3.2))

data = common.load('visualisation/data/Ds_from_boxcounting_first_quad_eleanorsmall016short')
D0_firstquad = data['Ds'][15]
print('DDs')
print(data['Ds'])
# raise Exception(f'D0 from first quad: {D0_firstquad}')

phi = 0.181
sigma = 2.01
z_mult = 2.00

rescale_y = True
# rescale_y = False

z_mult_str = f'{z_mult:.2f}a'
a = sigma/2
z = z_mult * a

hydrofile_wall          = '/data2/acarter/direct_to_H/' f'data/H_of_k_L300_numkbins60_single_wall_walllub_phi0.181_maxk10.0_montecarlo_full_zpos{z_mult_str}_zwidth0.0a_numis200.npz'
hydrofile_wall_particle  = '/data2/acarter/direct_to_H/' f'data/H_of_k_L200_numkbins60_single_wall_walllub_particlelub_phi0.181_maxk10.0_montecarlo_full_zpos1.20a_zwidth0.0a_numis100.npz'
hydrofile_wall_particle2 = '/data2/acarter/direct_to_H/' f'data/H_of_k_L200_numkbins60_single_wall_walllub_particlelub_phi0.181_maxk10.0_montecarlo_full_zpos1.05a_zwidth0.0a_numis100.npz'

visualisation.Ds_overlapped_mult.go(
    (
        dict(
            file='eleanorsmall016short',
            source='timescaleint_nofit_cropped_var',
            plot_index=np.index_exp[:16],
            msd='eleanorsmall016short',
            # msd=D0_firstquad,
            label='experiment',
        ),
        # dict(
        #     file='eleanorsmall016short',
        #     source='MSD_first',
        #     msd='eleanorsmall016short',
        #     # msd=D0_firstquad,
        #     label='experiment MSD',
        #     linestyle='-'
        # ),
        dict(
            file='eleanorsmall016long',
            source='timescaleint_nofit_cropped_var',
            plot_index=np.index_exp[16:-5],
            msd='eleanorsmall016short',
            label=None,
            color='tab:blue'
        ),
        dict(
            file='eleanorsmall016short',
            source='D_of_L_theory',
            # msd='eleanorsmall016short',
            msd = common.stokes_einstein_D(sigma),
            label = 'theory, no hydro',
            phi = phi,
            sigma = sigma,
        ),
        dict(
            file='eleanorsmall016short',
            source='D_of_L_theory_hydro',
            # plot_index=np.index_exp[30:-30],
            # msd_file='eleanorsmall016short',
            hydro_source_file=hydrofile_wall_particle,
            label = f'theory, hydro + lubrication $z=1.2a$',
            phi = phi,
            sigma = sigma,
            Lmax = 500,
        ),
        dict(
            file='eleanorsmall016short',
            source='D_of_L_theory_hydro',
            # plot_index=np.index_exp[30:-30],
            # msd_file='eleanorsmall016short',
            hydro_source_file=hydrofile_wall_particle2,
            label = f'theory, hydro + lubrication $z=1.05a$',
            phi = phi,
            sigma = sigma,
            Lmax = 500,
        ),
    ),
    ax = ax,
    plot_against_k=False,
    discrete_colors=True,
    allow_rescale_y=rescale_y,
    legend_fontsize=8
)
# ax.semilogy()

ax.legend(fontsize=7, loc='upper left')
ax.set_xlim(1e-1, 1e2)

common.save_fig(fig, 'workflows/figures/eleanorsmall_sophie.pdf', hide_metadata=True)
common.save_fig(fig, 'workflows/figures/eleanorsmall_sophie.png')

# python run.py --wall=single_wall_walllub_particlelub --phi=0.181 --method=montecarlo_full --zoffset=2 --a=1.005 --L=200 --num_is=100