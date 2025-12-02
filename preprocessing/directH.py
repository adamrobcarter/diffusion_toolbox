import common
import numpy as np
import matplotlib.pyplot as plt

def go_mesu(file_remote):
    file_base = file_remote.split('/')[-1][:-4]

    common.rsync(file_remote, f'preprocessing/directH/{file_base}.npz')

    go(file_base, f'preprocessing/directH/{file_base}.npz', remote_source_file=file_remote)

def go_isis(file):
    file_base = file.split('/')[-1][:-4]

    go(file_base, file, source_file=file)

def go(file_base, local_file, **to_save):
    data = common.load(local_file)
    """
        k=k_binned,
        H=H_binned,
        L=L,
        a=a,
        phi=phi,
        wall=wall,
    """
    k = data['k']
    H = data['H']
    particle_diameter = data['a'] * 2
    phi = data['phi']

    S = common.structure_factor_2d_hard_spheres(k, phi, particle_diameter)
    D0 = common.stokes_einstein_D(particle_diameter)

    # D_over_D0 = H / S
    D = D0 * H / S

    common.save_data(f'visualisation/data/Ds_from_H_theory_{file_base}.npz',
                    #  D_over_D0=D_over_D0,
                     Ds = D,
                     D_uncs = np.zeros_like(D),
                     ks=k,
                     D0 = D0,
                     particle_diameter=particle_diameter, phi=phi,
                     **to_save)
    
    fig, ax = plt.subplots()
    ax.scatter(k, H, label='H(k)')
    ax.scatter(k, S, label='S(k)')
    ax.scatter(k, D, label='D(k)')
    ax.set_xlabel('$k$')
    ax.semilogx()
    ax.semilogy()
    ax.legend()
    ax.grid()
    # ax.semilogy()
    common.save_fig(fig, f'preprocessing/figures_png/directH_{file_base}.png')

if __name__ == '__main__':
    # go_mesu('/store/cartera/direct_to_H/H_of_k_L2560_numkbins60_single_wall_phi0.114_maxk21.81661564992912_individualD_numis100.npz')
    # go_mesu('/store/cartera/direct_to_H/H_of_k_L2560_numkbins60_single_wall_phi0.114_maxk21.8_myRPY_numis10.npz')
    # go_mesu('/store/cartera/direct_to_H/H_of_k_L2560_numkbins60_single_wall_phi0.114_maxk21.8_myRPY_numis1000.npz')
    # go('/store/cartera/direct_to_H/H_of_k_L2560_numkbins60_open_phi0.114_maxk21.8_myRPY_numis1000.npz')

    # for z in range(10):
        # go_mesu(f'/store/cartera/direct_to_H/H_of_k_L2560_numkbins60_single_wall_phi0.114_maxk21.8_myRPY_z{z:.1f}a_numis1000.npz')

    go_isis('/data2/acarter/direct_to_H/' 'data/H_of_k_L2560_numkbins100_single_wall_'          'phi0.114_maxk10.0_theory_zpos1.15a_zwidth0.0a_rcutoff1.0_stepmult1.00_order3.npz')
    # go_isis('/data2/acarter/direct_to_H/' 'data/H_of_k_L2560_numkbins100_single_wall_lubrication_phi0.114_maxk10.0_theory_zpos1.06a_zwidth0.0a_rcutoff1.0_stepmult1.00_order3.npz')
    # go_isis('/data2/acarter/direct_to_H/' 'data/H_of_k_L2560_numkbins60_open_phi0.114_maxk10.0_theory_zpos0.00a_zwidth1.0a_rcutoff1.0_stepmult1.00_order3.npz')