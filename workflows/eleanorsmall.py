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

phi = 0.181
sigma = 2.01
z_mult = 2.00

rescale_y = True
# rescale_y = False

z_mult_str = f'{z_mult:.2f}a'
a = sigma/2
z = z_mult * a

hydrofile_wall          = '/data2/acarter/direct_to_H/' f'data/H_of_k_L300_numkbins60_single_wall_walllub_phi0.181_maxk10.0_montecarlo_full_zpos{z_mult_str}_zwidth0.0a_numis200.npz'
hydrofile_wall_particle = '/data2/acarter/direct_to_H/' f'data/H_of_k_L500_numkbins60_single_wall_walllub_particlelub_phi0.181_maxk10.0_montecarlo_full_zpos{z_mult_str}_zwidth0.0a_numis200.npz'

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
        # dict(
        #     file='eleanorsmall016short',
        #     source='boxcounting_first_quad',
        #     msd=D0_firstquad,
        #     label = 'boxcounnting first quad',
        # ),
        # dict(
        #     file='eleanorsmall016short',
        #     source='D_of_L_theory_hydro',
        #     # plot_index=np.index_exp[30:-30],
        #     # msd_file='eleanorsmall016short',
        #     hydro_source_file='/data2/acarter/direct_to_H/' 'data/H_of_k_L640_numkbins60_single_wall_walllub_phi0.181_maxk10.0_montecarlo_zpos1.05a_zwidth0.0a_numis10.npz',
        #     label = 'theory, hydro + wall lubrication old',
        #     phi = 0.181,
        #     sigma = 2.01,
        # ),
        dict(
            file='eleanorsmall016short',
            source='D_of_L_theory_hydro',
            # plot_index=np.index_exp[30:-30],
            # msd_file='eleanorsmall016short',
            hydro_source_file=hydrofile_wall,
            label = f'theory, hydro + wall lubrication z={z_mult_str}',
            phi = phi,
            sigma = sigma,
            Lmax = 500,
        ),
        # dict(
        #     file='eleanorsmall016short',
        #     source='D_of_L_theory_hydro',
        #     # plot_index=np.index_exp[30:-30],
        #     # msd_file='eleanorsmall016short',
        #     hydro_source_file='/data2/acarter/direct_to_H/' 'data/H_of_k_L100_numkbins60_single_wall_walllub_phi0.181_maxk10.0_montecarlo_full_zpos1.10a_zwidth0.0a_numis10.npz',
        #     label = 'theory, hydro + wall lubrication full L100',
        #     phi = 0.181,
        #     sigma = 2.01,
        # ),
        # dict(
        #     file='eleanorsmall016short',
        #     source='D_of_L_theory_hydro',
        #     # plot_index=np.index_exp[30:-30],
        #     # msd_file='eleanorsmall016short',
        #     hydro_source_file='/data2/acarter/direct_to_H/' 'data/H_of_k_L200_numkbins60_single_wall_walllub_particlelub_phi0.181_maxk10.0_montecarlo_full_zpos1.05a_zwidth0.0a_numis10.npz',
        #     label = 'theory, hydro + wall + particle lubrication $z=1.05a$',
        #     phi = 0.181,
        #     sigma = 2.01,
        # ),
        dict(
            file='eleanorsmall016short',
            source='D_of_L_theory_hydro',
            # plot_index=np.index_exp[30:-30],
            # msd_file='eleanorsmall016short',
            hydro_source_file=hydrofile_wall_particle,
            label = f'theory, hydro + wall + particle lubrication $z={z_mult_str}$',
            phi = phi,
            sigma = sigma,
            Lmax = 500,
        ),
        # dict(
        #     file='eleanorsmall016short',
        #     source='D_of_L_theory_hydro',
        #     # plot_index=np.index_exp[30:-30],
        #     # msd_file='eleanorsmall016short',
        #     hydro_source_file='/data2/acarter/direct_to_H/' 'data/H_of_k_L200_numkbins60_nohydro_phi0.181_maxk10.0_montecarlo_full_zpos1.10a_zwidth0.0a_numis10.npz',
        #     label = 'theory, nohydro H',
        #     linestyle='dotted',
        #     phi = 0.181,
        #     sigma = 2.01,
        # ),
        # dict(
        #     file='eleanorsmall016short',
        #     source='D_of_L_theory_hydro',
        #     # plot_index=np.index_exp[30:-30],
        #     # msd_file='eleanorsmall016short',
        #     hydro_source_file='/data2/acarter/direct_to_H/' 'data/H_of_k_L200_numkbins60_nohydro_walllub_phi0.181_maxk10.0_montecarlo_full_zpos1.10a_zwidth0.0a_numis10.npz',
        #     label = 'theory, nohydro + wall lubrication',
        #     phi = 0.181,
        #     sigma = 2.01,
        # ),
        # dict(
        #     file='eleanorsmall016short_bhwindow',
        #     source='f_t1',
        #     plot_index=np.index_exp[5:],
        #     msd='eleanorsmall016short',
        #     # msd=D0_firstquad,
        #     label='experiment',
        # ),
        # dict(
        #     file='eleanorsmall016long_bhwindow',
        #     source='f_t40',
        #     plot_index=np.index_exp[5:],
        #     msd='eleanorsmall016short',
        #     # msd=D0_firstquad,
        #     label='$f(k, t=40)$',
        # ),
        # dict(
        #     file='eleanorsmall016long_bhwindow',
        #     source='f_t160',
        #     plot_index=np.index_exp[5:],
        #     msd='eleanorsmall016short',
        #     # msd=D0_firstquad,
        #     label='$f(k, t=160)$',
        # ),
        # dict(
        #     file='eleanorsmall016long_bhwindow',
        #     source='f_t640',
        #     plot_index=np.index_exp[5:],
        #     msd='eleanorsmall016short',
        #     # msd=D0_firstquad,
        #     label='$f(k, t=640)$',
        # ),
        # dict(
        #     file='eleanorsmall016long_bhwindow',
        #     source='f_t2560',
        #     plot_index=np.index_exp[5:],
        #     msd='eleanorsmall016short',
        #     # msd=D0_firstquad,
        #     label='$f(k, t=2560)$',
        # ),
        # dict(
        #     file='eleanorsmall016short',
        #     source='D0Sk_theory',
        #     # plot_index=np.index_exp[30:-30],
        #     msd='eleanorsmall016short',
        #     label = f'$D(k)$ no hydro, D0 from MSD',
        #     phi = phi,
        #     sigma = sigma,
        # ),
        # dict(
        #     file='eleanorsmall016short',
        #     source='Dk_theory_hydro',
        #     # plot_index=np.index_exp[30:-30],
        #     msd='eleanorsmall016short',
        #     # msd = 5,
        #     hydro_source_file=hydrofile_wall_particle,
        #     label = f'$D(k)$ hydro, wall+particle lub, D0 from MSD',
        #     phi = phi,
        #     sigma = sigma,
        #     mult = 1.5
        # ),
    ),
    ax = ax,
    plot_against_k=False,
    discrete_colors=True,
    allow_rescale_y=rescale_y,
    legend_fontsize=8
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

common.save_fig(fig, 'workflows/figures/eleanorsmall.pdf', hide_metadata=True)
common.save_fig(fig, 'workflows/figures/eleanorsmall.png')

# python run.py --wall=single_wall_walllub_particlelub --phi=0.181 --method=montecarlo_full --zoffset=2 --a=1.005 --L=200 --num_is=100