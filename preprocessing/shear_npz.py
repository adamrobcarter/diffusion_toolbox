import numpy as np
import common
import time, os
import tqdm
import subprocess
import itertools
from scipy.spatial.transform import Rotation
import matplotlib.pyplot as plt
import numpy.linalg
import warnings

def go(infile, outfile, nth_timestep=1, max_time=None,
       extra_source_file=None, skiprows=None, multiply_time_by=None, **kwargs):

    print(f'loading npz, last modified {common.get_last_modified_time(infile)} ago')
    t0 = time.time()
    data = common.load(infile)
    particles = data['particles']
    t1 = time.time()
    print(particles.shape, particles.dtype,  f'loaded {common.format_bytes(particles.nbytes)} in {t1-t0:.1f}s at {common.format_bytes(particles.nbytes/(t1-t0))}/sec')

    TIME_COLUMN = 3
    ID_COLUMN = 4

    if np.isnan(particles).any():
        print('NAN FOUND')
        nans = np.isnan(particles).any(axis=1)
        particles = particles[~nans, :]
    
    print('xy min max', particles[:, 0].min(), particles[:, 0].max(), particles[:, 1].min(), particles[:, 1].max())
    
    times = np.unique(particles[:, TIME_COLUMN])
    num_timesteps = times.size
    last_timestep = times[-1]

    particles[:, TIME_COLUMN] -= particles[:, TIME_COLUMN].min()
    assert particles[:, TIME_COLUMN].min() == 0

    frame_time_deltas = particles[1:, TIME_COLUMN] - particles[:-1, TIME_COLUMN]
    assert np.all(frame_time_deltas >= 0)

    assert particles.shape[1] == 9 # we have quaternion for rotation

    last_rot = np.eye(3) # identity matrix
    last_angles = np.array([0.0, 0.0, 0.0])

    angles = np.full((particles.shape[0], 3), np.nan)

    assert np.all(particles[:, ID_COLUMN] == 0)

    for i in tqdm.trange(particles.shape[0]):
        rot = Rotation.from_quat(particles[i, [5, 6, 7, 8]], scalar_first=True).as_matrix() # needs scipy >= 1.15
    
        rot_change = rot @ np.linalg.inv(last_rot) # rotation change from last frame
        # rot_change = np.linalg.inv(last_rot) @ rot # rotation change from last frame

        last_rot = rot
        
        rot_x = np.array([1, 0, 0]) @ rot_change
        rot_y = np.array([0, 1, 0]) @ rot_change
        rot_z = np.array([0, 0, 1]) @ rot_change
        # we now have 3D vectors that we can project

        about_x_from_y = - np.arctan2(rot_y[2], rot_y[1])
        about_x_from_z =   np.arctan2(rot_z[1], rot_z[2])
        if not np.isclose(about_x_from_y, about_x_from_z, atol=1e-1):
            warnings.warn(f'about_x_from_y-about_x_from_z={about_x_from_y-about_x_from_z}')
        # assert np.isclose(about_x_from_y, about_x_from_z, atol=1e-1), f'about_x_from_y-about_x_from_z={about_x_from_y-about_x_from_z}'

        about_y_from_z = - np.arctan2(rot_z[0], rot_z[2])
        about_y_from_x =   np.arctan2(rot_x[2], rot_x[0])
        if not np.isclose(about_y_from_z, about_y_from_x, atol=1e-1):
            warnings.warn(f'about_y_from_z-about_y_from_x={about_y_from_z-about_y_from_x}')
        # assert np.isclose(about_y_from_z, about_y_from_x, atol=1e-1), f'about_y_from_z-about_y_from_x={about_y_from_z-about_y_from_x}'

        about_z_from_x = - np.arctan2(rot_x[1], rot_x[0])
        about_z_from_y =   np.arctan2(rot_y[0], rot_y[1])
        if not np.isclose(about_z_from_x, about_z_from_y, atol=1e-1):
            warnings.warn(f'about_z_from_x-about_z_from_y={about_z_from_x-about_z_from_y}')
        # assert np.isclose(about_z_from_x, about_z_from_y, atol=1e-1), f'about_z_from_x-about_z_from_y={about_z_from_x-about_z_from_y}'

        angles[i, 0] = last_angles[0] + about_x_from_y
        angles[i, 1] = last_angles[1] + about_y_from_z
        angles[i, 2] = last_angles[2] + about_z_from_x

        last_angles[0] = angles[i, 0]
        last_angles[1] = angles[i, 1]
        last_angles[2] = angles[i, 2]



    # rot = Rotation.from_quat(particles[:, [5, 6, 7, 8]], scalar_first=True).as_matrix() # needs scipy >= 1.15
    
    # rot_x = np.array([1, 0, 0]) @ rot
    # rot_y = np.array([0, 1, 0]) @ rot
    # rot_z = np.array([0, 0, 1]) @ rot
    # # we now have 3D vectors that we can project
    # about_z_from_x =   np.arctan2(rot_x[:, 1], rot_x[:, 0])
    # about_y_from_x = - np.arctan2(rot_x[:, 2], rot_x[:, 0])

    # about_x_from_y =   np.arctan2(rot_y[:, 2], rot_y[:, 1])
    # about_z_from_y = - np.arctan2(rot_y[:, 0], rot_y[:, 1])

    # about_y_from_z =   np.arctan2(rot_z[:, 0], rot_z[:, 2])
    # about_x_from_z = - np.arctan2(rot_z[:, 1], rot_z[:, 2])

    # all_angles = np.column_stack((about_z_from_x, about_z_from_y, about_y_from_x, about_y_from_z, about_x_from_y, about_x_from_z, particles[:, 3], particles[:, 4])) # add time and ID columns
    # all_angles = common.periodic_unwrap(all_angles, 6, [0, 1, 2, 3, 4, 5], [2*np.pi, 2*np.pi, 2*np.pi, 2*np.pi, 2*np.pi, 2*np.pi])
    
    # fig, ax = plt.subplots(1, 1)
    # ax.set_xlim(0, 1000)
    # ax.set_ylim(-10, 10)
    # ax.plot(times, all_angles[:, 0], label='about z from x')
    # ax.plot(times, all_angles[:, 1], label='about z from y')
    # ax.plot(times, angles[:, 2], label='about y from x')    
    # ax.plot(times, angles[:, 3], label='about y from z')
    # ax.plot(times, angles[:, 4], label='about x from y')
    # ax.plot(times, angles[:, 5], label='about x from z')
    # ax.plot(times, about_x_from_y, label='about x from y') # bad tz
    # ax.plot(times, about_x_from_z, label='about x from z') # bad ty
    # ax.plot(times, about_y_from_z, label='about y from z') # bad tx
    # ax.plot(times, about_y_from_x, label='about y from x') # bad tz
    # ax.plot(times, about_z_from_x, label='about z from x') # bad ty
    # ax.plot(times, about_z_from_y, label='about z from y') # bad tx
    # ax.plot(times, rot_x[:, 0], label=r'$r_x{}_x$')
    # ax.plot(times, rot_x[:, 1], label=r'$r_x{}_y$')
    # ax.plot(times, rot_x[:, 2], label=r'$r_x{}_z$')
    # ax.plot(times, rot_y[:, 0], label=r'$r_y{}_x$')
    # ax.plot(times, rot_y[:, 1], label=r'$r_y{}_y$')
    # ax.plot(times, rot_y[:, 2], label=r'$r_y{}_z$')
    # ax.plot(times, rot_z[:, 0], label=r'$r_z{}_x$')
    # ax.plot(times, rot_z[:, 1], label=r'$r_z{}_y$')
    # ax.plot(times, rot_z[:, 2], label=r'$r_z{}_z$')
    # ax.legend()
    # common.save_fig(fig, f'particle_linking/figures_png/rot_{outfile}.png', dpi=300)

    particles[:, 5] = angles[:, 0]
    particles[:, 6] = angles[:, 1]
    particles[:, 7] = angles[:, 2]
    particles = particles[:, [0, 1, 2, 3, 4, 5, 6, 7]] # remove extra column

    # # unwrap the angles
    # # particles = common.periodic_unwrap(particles, 3, [5, 6, 7], [2*np.pi, np.pi, 2*np.pi])

    newdata = common.copy_not_particles(data) # copy so we can modify
    newdata['particles'] = particles

    common.save_data(f'particle_detection/data/particles_{outfile}.npz',
        **newdata,
        max_time_hours=round(last_timestep/60/60, 2),
        source_file=infile, extra_source_file=extra_source_file,
        dimension=3,
        **kwargs
    )

    common.save_data(f'particle_linking/data/trajs_{outfile}.npz',
        **newdata,
        max_time_hours=round(last_timestep/60/60, 2),
        source_file=infile, extra_source_file=extra_source_file,
        dimension=3,
        **kwargs
    )

    # end_timestep = data[:, 2].max() // 8
    # data_small = data[data[:, 2] < end_timestep, :]
    # common.save_data(f'particle_linking/data/trajs_{outfile}_div8.npz',
    #     particles=data_small,
    #     time_step=dt, particle_diameter=particle_diameter,
    #     max_time_hours=round(end_timestep/60/60, 2),
    #     source_file=infile, extra_source_file=extra_source_file,
    #     **kwargs
    # )
        
    # end_timestep = data[:, 2].max() // 64
    # data_small = data_small[data_small[:, 2] < end_timestep, :]
    # common.save_data(f'particle_linking/data/trajs_{outfile}_div64.npz',
    #     particles=data_small,
    #     time_step=dt, particle_diameter=particle_diameter,
    #     max_time_hours=round(end_timestep/60/60, 2),
    #     source_file=infile,extra_source_file=extra_source_file,
    #     **kwargs
    # )
    print()
    

def go_mesu_shear(filepath, suffix='', skip_rsync=False, **kwargs):
    filename = filepath.split('/')[-1]
    if not skip_rsync:
        rsync_command = ['rsync', f'cartera@login.mesu.sorbonne-universite.fr:{filepath}', f'raw_data/mesu/{filename}', '--progress']
        print('launching', ' '.join(rsync_command))
        subprocess.run(rsync_command, check=True)

    internalname = filename.split('.npz')[0]

    go(
        f'raw_data/mesu/{filename}',
        f'sim_{internalname}{suffix}',
        extra_source_file=filepath,
    )

if __name__ == '__main__':
    # go_mesu_shear('/store/cartera/shear/shear0.0_RHS.npz')
    # go_mesu_shear('/store/cartera/shear/shear0.0_theta45_RHS.npz')

    # go_mesu_shear('/store/cartera/shear/shear0.0_T0_theta45_RHS.npz')
    # go_mesu_shear('/store/cartera/shear/shear0.0_T0_theta15_RHS.npz')
    # go_mesu_shear('/store/cartera/shear/shear0.0_T0_theta30_RHS.npz')
    
    
    
    # for nblobs in [2562]:
    # # for nblobs in [42, 162, 642, 2562]:
    #     for theta in [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85]:
    #         theta = int(theta)
    #         go_mesu_shear(f'/store/cartera/shear/shear0.0_T0_theta{theta}_RHS_nblobs{nblobs}.npz')
    
    # for nblobs in [42, 162, 642, 2562]:
    #     for theta in [0, 2, 4, 6, 8, 10, 12, 14]:
    #         theta = int(theta)
    #         go_mesu_shear(f'/store/cartera/shear/shear0.0_T0_theta{theta}_RHS_nblobs{nblobs}.npz')

    # for shear in [0.0, 0.25, 0.5, 0.75, 1.0]:
    #     go_mesu_shear(f'/store/cartera/shear/shear{shear}_T0_RHS_nblobs642.npz')

    # go_mesu_shear('/store/cartera/shear/shear0.0_theta0_RHS_nblobs642.npz')
    # go_mesu_shear('/store/cartera/shear/shear0.0_theta0_RHS_nblobs42.npz')
    
    # for nblobs in [42, 162, 642, 2562]:
    #     go_mesu_shear(f'/store/cartera/shear/shear0.0_theta0_RFD_nblobs{nblobs}.npz')

        
    # for n_blobs in [42, 162]:
    #     for T in [296, 148]:
    #         go_mesu_shear(f'/store/cartera/shear/shear0.0_T{T}_theta0_RFD_nblobs{n_blobs}.npz')


        
    # t_max   = [10000]
    # t_save  = [0.01]
    # dt      = [0.01]
    # # dt      = [0.01, 0.005]
    # shear   = [0.080357]
    # # shear   = [0.0, 0.2, 0.4, 0.6, 0.8]
    # gravity = [True]
    # wall    = [True]
    # T       = [296]
    # # method = 'RFDp'
    # # method  = ['RHS', 'RFDp', 'RFDm']
    # method  = ['RHS']
    # suffix  = ['']
    # # n_blobs = [42, 162]
    # n_blobs = [162]
    # theta   = [10]
    # # theta   = [0, 2, 4, 6, 8, 10]


    # for combo in itertools.product(t_max, t_save, dt, shear, gravity, wall, theta, T, method, suffix, n_blobs):
    #     t_max, t_save, dt, shear, gravity, wall, theta, T, method, suffix, n_blobs = combo
    #     go_mesu_shear(f'/store/cartera/shear/shear{shear}_T{T}_theta{theta}_{method}_nblobs{n_blobs}_dt{dt}_tmax{t_max}.npz')


    # testing rotational diffuision
    # go_mesu_shear(f'/store/cartera/shear/shear0_T296_nograv_theta0_nowall_EMRFD_nblobs162_dt0.01_tmax1000.npz', time_step=0.01) # time_step not needed in future
    # go_mesu_shear(f'/store/cartera/shear/shear0.080357_T296_theta10_EMRFD_nblobs162_dt0.01_tmax1000.npz', time_step=0.01) # time_step not needed in future
    # go_mesu_shear(f'/store/cartera/shear/shear0_T296_theta0_EMRFD_nblobs162_dt0.01_tmax1000.npz', time_step=0.01) # time_step not needed in future

    # monolayer
    # go_mesu_shear('/store/cartera/shear/shear0_T296_theta10_L50_phi0.5_Trap_nblobs42_dt0.01_tmax10.npz')
    # go_mesu_shear('/store/cartera/shear/shear0_T296_theta10_L50_phi0.01_Trap_nblobs42_dt0.01_tmax10.npz')

    
    # go_mesu_shear('/store/cartera/shear/shear0.080357_T298_theta10_EMRFD_nblobs42_dt0.005_tmax1000.npz')
    # go_mesu_shear('/store/cartera/shear/shear0.080357_T298_theta0_EMRFD_nblobs42_dt0.005_tmax1000.npz')
    # go_mesu_shear('/store/cartera/shear/shear0_T298_theta10_EMRFD_nblobs42_dt0.005_tmax1000.npz')
    # go_mesu_shear('/store/cartera/shear/shear0_T298_theta0_EMRFD_nblobs42_dt0.005_tmax1000.npz')
    # go_mesu_shear('/store/cartera/shear/shear0_T298_nograv_nowall_EMRFD_nblobs42_dt0.005_tmax1000.npz')
    # be good to compare the wall one with brenner
    # these with tmax10000 are calculating

    # for theta in [0, 2, 4, 6, 8, 10]:
    #     go_mesu_shear(f'/store/cartera/shear/shear0_T298_theta{theta}_EMRFD_nblobs42_dt0.005_tmax1000_torque.npz')

    # arbitrary torques nowall nograv
    # go_mesu_shear(f'/store/cartera/shear/shear0_T0_nograv_nowall_EMRFD_nblobs42_dt0.005_tmax100_torquexy.npz')

    # single particle rotational stuff
    # go_mesu_shear('/store/cartera/shear/shear0_T296_nograv_nowall_EMRFD_nblobs42_dt0.005_tmax10000.npz')
    # go_mesu_shear('/store/cartera/shear/shear0_T296_theta0_EMRFD_nblobs42_dt0.005_tmax10000.npz')
    # go_mesu_shear('/store/cartera/shear/shear0_T296_theta10_EMRFD_nblobs42_dt0.005_tmax10000.npz')
    # go_mesu_shear('/store/cartera/shear/shear0.080357_T296_theta0_EMRFD_nblobs42_dt0.005_tmax10000.npz')
    # go_mesu_shear('/store/cartera/shear/shear0.080357_T296_theta10_EMRFD_nblobs42_dt0.005_tmax10000.npz')
    # go_mesu_shear('/store/cartera/shear/shear0_T296_nograv_nowall_trap_nblobs42_dt0.005_tmax10000.npz')
    # go_mesu_shear('/store/cartera/shear/shear0_T296_nograv_nowall_EMRFD_nblobs42_dt0.005_tmax1000.npz')
    # go_mesu_shear('/store/cartera/shear/shear0_T296_nograv_nowall_EMRFD_nblobs42_dt0.005_tmax10000.npz')

    # dec 2025
    # single particle
    # go_mesu_shear('/store/cartera/shear/shear0.080357_T300_theta10_EMmid_nblobs42_dt0.005_tmax36000.npz')
    # go_mesu_shear('/store/cartera/shear/shear0_T300_theta0_EMmid_nblobs42_dt0.005_tmax36000.npz')
    # go_mesu_shear('/store/cartera/shear/shear0_T300_theta10_EMmid_nblobs42_dt0.005_tmax36000.npz')
    go_mesu_shear('/store/cartera/shear/shear0.080357_T300_theta0_EMmid_nblobs42_dt0.005_tmax36000.npz')