import numpy as np
import common
import time, os
import tqdm

# eleanorlong is
# x min max 0.6982900415999999 293.91458112
# y min max 0.7108168319999999 367.6492512
# num_timesteps 71997.0

# def read_file_skip_timesteps(file, nth_timesteps, max_frames, expected_particles_per_frame):
#     arraysize = 0
#     if max_frames:
#         progress = tqdm.tqdm(total=max_frames)
#     else:
#         progress = tqdm.tqdm()

#     # num_lines = 0

#     with open(file) as f:
#         for line in f:
#             # num_lines += 1

#             if nth_timesteps == 1:
#                 yield line

#                 if max_frames == None:
#                     progress.update()

#                 else:
#                     # we can use the density to guess the current frame number
#                     # which will make it quicker as we won't have to parse the data
#                     # actually it doesn't seem to be quicker rip
#                     # t = num_lines / expected_particles_per_frame
                    
#                     # 4.04it/s by the end

#                     # old option - is this slower?
#                     parts = line.split()
#                     t = int(parts[2])
#                     progress.n = t # don't do progress.refresh here! it's very slow
                    
#                     if t > max_frames:
#                         break

#                 arraysize += 4 * 4
            
#             else:
#                 parts = line.split()
#                 t = int(parts[2])

#                 if t % nth_timesteps == 0:
#                     yield line
#                     arraysize += 4 * 4
                    
#                 if max_frames == None:
#                     progress.update()
#                 else:
#                     progress.n = t # don't do progress.refresh here! it's very slow
                    
#                     if t > max_frames:
#                         break

#             progress.set_description(f'loaded {common.format_bytes(arraysize)}')

# def read_file_tqdm(file):
#     arraysize = 0
#     with open(file) as f:
#         for line in (progress := tqdm.tqdm(f, total=1e10)):
#             yield line
            
#             arraysize += 4 * 4
#             progress.set_description(f'loaded {common.format_bytes(arraysize)}')


def go(infile, outfile, orig_width, out_width, dt, pack_frac_given, particle_diameter, nth_timestep=1, max_time=None):
    print(f'loading raw file, last modified {common.get_last_modified_time(infile)} ago')

    density = 4/np.pi * pack_frac_given / particle_diameter**2
    expected_particles_per_frame = density * orig_width**2
    if max_time:
        max_frame = max_time * dt
        max_rows = int((max_frame + 1) * expected_particles_per_frame)
        print('max_rows', max_rows)

        max_items = max_rows * 3
    else:
        max_rows = None
        max_items = -1

    t0 = time.time()
    # data = np.loadtxt(infile)
    # data = np.loadtxt(infile, dtype=np.float32)
    # data = np.loadtxt(read_file_skip_timesteps(infile, nth_timesteps, max_frames, expected_particles_per_frame), dtype=np.float32)
    if infile.endswith('.txt'):
        data = np.loadtxt(infile, max_rows=max_rows, dtype=np.float32)
    elif infile.endswith('.bin'):
        data = np.fromfile(infile, dtype=np.float32, count=max_items) # there is a parameter to this function for not loading the whole array
        data = data.reshape((-1, 3))
    t1 = time.time()
    print(data.shape, data.dtype,  f'loaded {common.format_bytes(data.nbytes)} in {t1-t0:.1f}s at {common.format_bytes(data.nbytes/(t1-t0))}/sec')
    
    # discard the non-nth step rows
    if nth_timestep > 1:
        print('doing nth-timestep')
        data = data[data[:, 2] % nth_timestep == 0, :]
        dt *= nth_timestep
        # we gotta readjust the timestep number now we discarded some timesteps
        data[:, 2] /= nth_timestep

    # we should delete the last timestep cause it may be incomplete
    last_timestep = data[:, 2].max()
    data = data[data[:, 2] != last_timestep, :]
    
    # print(f'loaded in {t1-t0:.0f}s. shape', data.shape, common.arraysize(data))
    # num_timesteps = data[:, 2].max()+1
    times = np.unique(data[:, 2])
    num_timesteps = times.size
    print(times, num_timesteps)
    print(f'{num_timesteps:.0f} timesteps, {data[:, 2].max()*dt/60/60:.1f} hours')
    assert num_timesteps > 30

    
    avg_particles_per_frame = data.shape[0] / num_timesteps
    density = avg_particles_per_frame / orig_width**2
    print('avg part per frame', avg_particles_per_frame, 'L^2', orig_width**2)
    pack_frac_calced = np.pi/4 * density * particle_diameter**2
    print(f'pack_frac_calced={pack_frac_calced:.4f}, pack_frac_given={pack_frac_given:.4f}')
    assert np.isclose(pack_frac_calced, pack_frac_given, rtol=0.1), f'pack frac calced {pack_frac_calced}, given {pack_frac_given}'

    # print('removing some timesteps')
    # keep_rows = data[:, 2] % nth_timesteps == 0
    # data = data[keep_rows, :]
    # print(f'removed {keep_rows.sum()/keep_rows.size:.3f}')

    data[:, 2] -= data[:, 2].min()
    assert data[:, 2].min() == 0

    # num_timesteps_before = int(data[:, 2].max() + 1)
    # print('num timesteps before crop', num_timesteps_before)


    print('time1:')
    # common.term_hist(data[:, 2])

    # not_too_long = data[:, 2] < 26477
    # print(f'keeping {not_too_long.sum()/not_too_long.size:.2f} (too long)')
    # data = data[not_too_long, :]

    assert orig_width == out_width, 'at the moment'

    # print('raw histogram:')
    # common.term_hist(data[:, 1])
    data[:, 0] = data[:, 0] % orig_width
    data[:, 1] = data[:, 1] % orig_width

    assert orig_width == out_width # for the mo

    # print('modded into window:')
    # common.term_hist(data[:, 1])
    print('x,y min,max', data[:, 0].min(), data[:, 0].max(), data[:, 1].min(), data[:, 1].max())
    keep = ( data[:, 0] >= 0 ) & ( data[:, 0] <= out_width ) &  ( data[:, 1] >= 0 ) & ( data[:, 1] <= out_width ) # I think they should be "<" not "<=" but for some reason this is failing on sim_nohydro_034_L1280
    assert keep.sum() == keep.size

    # if orig_width > out_width:
    #     num_timesteps = int(data[:, 2].max() + 1)
    #     common.save_data(f'particle_detection/data/particles_{outfile}_nocrop.npz', particles=data,
    #                 time_step=0.5, particle_diameter=2.79, pack_frac_given=0.34,
    #                 window_size_x=orig_width, window_size_y=orig_width,
    #                 num_timesteps=num_timesteps)

    # print('cropped into window:')
    # common.term_hist(data[keep, 1])
    # print(f'keeping {keep.sum()/keep.size:.2f} (inside crop)')
    # data = data[keep, :]

    # print('time2:')
    # common.term_hist(data[:, 2])

    common.save_data(f'particle_detection/data/particles_{outfile}.npz',
        particles=data,
        time_step=dt, particle_diameter=particle_diameter, pack_frac_given=pack_frac_given, pack_frac=pack_frac_calced,
        window_size_x=out_width, window_size_y=out_width, max_time_hours=round(last_timestep*dt/60/60, 2),
        source_file=infile, density=density,
    )

    if data.shape[1] == 4:
        common.save_data(f'particle_linking/data/trajs_{outfile}.npz',
            particles=data,
            time_step=dt, particle_diameter=particle_diameter, pack_frac_given=pack_frac_given, pack_frac=pack_frac_calced,
            window_size_x=out_width, window_size_y=out_width, max_time_hours=round(last_timestep*dt/60/60, 2),
            source_file=infile, density=density,
        )
    
    if data.size > 5e7:
        end_timestep = data[:, 2].max() // 8
        data_small = data[data[:, 2] < end_timestep, :]

        common.save_data(f'particle_detection/data/particles_{outfile}_div8.npz',
            particles=data_small,
            time_step=dt, particle_diameter=2.79, pack_frac_given=pack_frac_given, pack_frac=pack_frac_calced,
            window_size_x=out_width, window_size_y=out_width, max_time_hours=round(end_timestep*dt/60/60, 2),
            source_file=infile, density=density,
        )

        if data.shape[1] == 4:
            common.save_data(f'particle_linking/data/trajs_{outfile}_div8.npz',
                particles=data_small,
                time_step=dt, particle_diameter=particle_diameter, pack_frac_given=pack_frac_given, pack_frac=pack_frac_calced,
                window_size_x=out_width, window_size_y=out_width, max_time_hours=round(end_timestep*dt/60/60, 2),
                source_file=infile, density=density,
            )

            
        end_timestep = data[:, 2].max() // 64
        data_small = data_small[data_small[:, 2] < end_timestep, :]

        common.save_data(f'particle_detection/data/particles_{outfile}_div64.npz',
            particles=data_small,
            time_step=dt, particle_diameter=2.79, pack_frac_given=pack_frac_given, pack_frac=pack_frac_calced,
            window_size_x=out_width, window_size_y=out_width, max_time_hours=round(end_timestep*dt/60/60, 2),
            source_file=infile,
        )

        if data.shape[1] == 4:
            common.save_data(f'particle_linking/data/trajs_{outfile}_div8.npz',
                particles=data_small,
                time_step=dt, particle_diameter=particle_diameter, pack_frac_given=pack_frac_given, pack_frac=pack_frac_calced,
                window_size_x=out_width, window_size_y=out_width, max_time_hours=round(end_timestep*dt/60/60, 2),
                source_file=infile,
            )
    print()
    
# 0.34
# go('raw_data/brennan/spec_softetakt_long_run_dtau_0.025_nsave_2.suspension_phi_0.34_L_320_modified.txt', 'brennan_hydro_034', 320, 320)
# go('raw_data/brennan/noHydro2D_Leim_run_dt_0.0125_nsave_40.suspension_phi_0.34_L_640_modified.txt', 'brennan_nohydro_034', 640, 320)
# 0.66
# go('/data2/acarter/Spectral_Sophie_Boxes/data/spec_softetakt_long_run_dtau_0.025_nsave_4.suspension_phi_0.66_L_288_modified.txt', 'brennan_hydro_066', 288, 288)
# go('/data2/acarter/Spectral_Sophie_Boxes/data/noHydro2D_Leim_run_dt_9.765625e5_nsave_256_long.suspension_phi_0.66_L_640_eq_modified.txt', 'brennan_nohydro_066', 640, 288)

# go('/data2/acarter/sim/RigidMultiblobsWall/Lubrication/Lubrication_Examples/Monolayer/data/nohydro2D_L1280_dt4.suspension_phi_0.34_L_1280_modified.txt', 'sim_nohydro_034_L1280', 1280, 1280, 4)
# go('/data2/acarter/sim/RigidMultiblobsWall/Lubrication/Lubrication_Examples/Monolayer/data/nohydro2D_L640_dt1.suspension_phi_0.34_L_640_modified.txt',   'sim_nohydro_034_L640',  640,  640,  1)

# from brennnan just before paper
# go('raw_data/brennan/spec_softetakt_dt_0.25_nsave_2.suspension_phi_0.016_L_320_modified.txt', 'brennan_hydro_002_L320', 320, 320, 0.5, 0.016, 2.972)
# go('raw_data/brennan/spec_softetakt_dt_0.25_nsave_2.suspension_phi_0.016_L_640_modified.txt', 'brennan_hydro_002_L640', 640, 640, 0.5, 0.016, 2.972)
# go('raw_data/brennan/spec_softetakt_dt_0.25_nsave_2.suspension_phi_0.114_L_320_modified.txt', 'brennan_hydro_011_L320', 320, 320, 0.5, 0.114, 2.972)
# go('raw_data/brennan/spec_softetakt_dt_0.25_nsave_64.suspension_phi_0.016_L_320_modified.txt', 'brennan_hydro_002_L320_longer', 320, 320, 16, 0.016, 2.972)

# from anubis
# go('raw_data/anubis/nohydro2D_L640_dt0.5_s2.972.suspension_phi0.016_L640_s2.972.bin', 'sim_nohydro_002_L640',        640, 640, 0.5, 0.016, 2.972)
# go('raw_data/anubis/nohydro2D_L640_dt16_s2.972.suspension_phi0.016_L640_s2.972.bin',  'sim_nohydro_002_L640_longer', 640, 640, 16,  0.016, 2.972)
# go('raw_data/anubis/nohydro2D_L160_dt0.5_s2.972.suspension_phi0.114_L160_s2.972.bin', 'sim_nohydro_011_L160',        160, 160, 0.5, 0.114, 2.972)
# go('raw_data/anubis/nohydro2D_L160_dt16_s2.972.suspension_phi0.114_L160_s2.972.bin',  'sim_nohydro_011_L160_longer', 160, 160, 16,  0.114, 2.972)
# from anubis after first paper submission
# go('raw_data/anubis/nohydro2D_L1280_dt64_s2.972.suspension_phi0.114_L1280_s2.972.bin',  'sim_nohydro_011_L1280_longer', 1280, 1280, 64,  0.114, 2.972)

datas = [
    # L  dt  phi
    # (320, 1, 0.01),
    # (640, 1, 0.01),
    # (1280, 4, 0.01),
    # (160, 1, 0.01),
    # (160, 4, 0.1),
    # (320, 4, 0.1),
    # (640, 8, 0.1, '', 1, None), # this was the one that was done again
    # (640, 0.5, 0.1, '_short'),
    # (320, 8, 0.1, '', 1, 24*60*60), # this is 24h
    # (320, 8, 0.01, '', 1, 24*60*60),
    # (640, 8, 0.01, '', 1, 24*60*60),
    # (320, 8, 0.1, '', 1, 24*60*60),
    # (640, 8, 0.1, '', 1, 24*60*60),
    # (544, 8, 0.1, '', 1, 24*60*60),
    # (544, 0.5, 0.1, '_dt2', 4, 24*60*60),
    # (800, 8, 0.02, '', 1, 24*60*60),
    # (320, 8, 0.114, '_long', 1, None, 2.972),
    # (640, 0.5, 0.10, '', 1, None, 2.79),
    # (640, 8, 0.10, '_longer', 2, None, 2.79),

    # (320, 0.5, 0.114, '',        1, None, 2.972), # sim_nohydro_011_L320
    # (320, 0.5, 0.016, '',        1, None, 2.972), # sim_nohydro_002_L320
    # (320, 16, 0.114, '_longer',  1, None, 2.972), # sim_nohydro_011_L320_longer
    # (320, 32, 0.114, '_longest', 1, None, 2.972), # sim_nohydro_011_L320_longest
    # (320, 16, 0.016, '_longer',  1, None, 2.972), # sim_nohydro_002_L320_longer
    # (640, 0.5, 0.016, '',        1, None, 2.972), # sim_nohydro_002_L640
    # (640, 0.5, 0.114, '',        1, None, 2.972), # sim_nohydro_011_L640
    # (640, 16, 0.114, '_longer',  1, None, 2.972), # sim_nohydro_011_L640_longer

    # after first submission of paper
    # (1280, 0.5, 0.114, '',        1, None, 2.972), # sim_nohydro_011_L1280
    
    # reimport
    # (544, 0.5, 0.1, '',        1, None, 2.79), # sim_nohydro_011_L320
]

for L, dt, phi, suffix, nth_timestep, max_time, particle_diameter in datas:
    phistr = f'{phi*100:.0f}'.zfill(3)
    # go(f'/data2/acarter/sim/RigidMultiblobsWall/Lubrication/Lubrication_Examples/Monolayer/data/nohydro2D_L{L}_dt{dt}.suspension_phi_{phi}_L_{L}_modified.txt', f'sim_nohydro_{phistr}_L{L}{suffix}',  L, L, dt, phi, nth_timestep, max_time)
    
    if particle_diameter == 2.972:
        # new ones
        file = f'/data2/acarter/sim/RigidMultiblobsWall/Lubrication/Lubrication_Examples/Monolayer/data/nohydro2D_L{L}_dt{dt}_s2.972.suspension_phi{phi}_L{L}_s2.972.bin'
    else:
        # old ones
        file = f'/data2/acarter/sim/RigidMultiblobsWall/Lubrication/Lubrication_Examples/Monolayer/data/nohydro2D_L{L}_dt{dt}.suspension_phi_{phi}_L_{L}_modified.txt'
       
    go(
        infile = file,
        outfile = f'sim_nohydro_{phistr}_L{L}{suffix}',
        orig_width=L, 
        out_width=L, 
        dt=dt, 
        pack_frac_given=phi, nth_timestep=nth_timestep, max_time=max_time, particle_diameter=particle_diameter
    )


# new times (not frames) method
datas2 = [
    # (1280, 64, 0.5, 0.114, '_mixt', 1, None, 2.972),
    # (640, 16, 0.5, 0.114, '_mixt', 1, None, 2.972),
]

for L, t1, t2, phi, suffix, nth_timestep, max_time, particle_diameter in datas2:
    phistr = f'{phi*100:.0f}'.zfill(3)
    # go(f'/data2/acarter/sim/RigidMultiblobsWall/Lubrication/Lubrication_Examples/Monolayer/data/nohydro2D_L{L}_dt{dt}.suspension_phi_{phi}_L_{L}_modified.txt', f'sim_nohydro_{phistr}_L{L}{suffix}',  L, L, dt, phi, nth_timestep, max_time)
    go(
        # new ones:
       f'/data2/acarter/sim/RigidMultiblobsWall/Lubrication/Lubrication_Examples/Monolayer/data/nohydro2D_L{L}_t{t1}_{t2}.suspension_phi{phi}_L{L}_s2.972.bin',
       # old ones:
    #    f'/data2/acarter/sim/RigidMultiblobsWall/Lubrication/Lubrication_Examples/Monolayer/data/nohydro2D_L{L}_dt{dt}.suspension_phi_{phi}_L_{L}_modified.txt',
       f'sim_nohydro_{phistr}_L{L}{suffix}',
       orig_width=L, 
       out_width=L, 
       dt=1,
       pack_frac_given=phi, nth_timestep=nth_timestep, max_time=max_time, particle_diameter=particle_diameter)


###### mesu single t
datas = [
    ('/data2/acarter/toolbox/raw_data/mesu/hydro_dt0.2_phi0.114_L320.bin', '_dt0.2'),
]

for filepath, suffix in datas:
    nth_timestep = 1
    max_time = None
    particle_diameter = 2.972

    filename = filepath.split('/')[-1]
    # print('use the following in a new shell')
    # print(f'rsync cartera@login.mesu.sorbonne-universite.fr:{filepath} raw_data/mesu/{filename}')
    # input()

    filename_no_ext = filename.split('.bin')[0]

    L   = int  (filename_no_ext.split('_L' )[1].split('_')[0])
    phi = float(filename_no_ext.split('phi')[1].split('_')[0])
    print(L, phi)

    phistr = f'{phi*100:.0f}'.zfill(3)

    go(
        f'raw_data/mesu/{filename}',
        f'sim_hydro_{phistr}_L{L}{suffix}',
        orig_width=L, 
        out_width=L, 
        dt=1,
        pack_frac_given=phi, nth_timestep=nth_timestep, max_time=max_time, particle_diameter=particle_diameter
    )




# infile, outfile, orig_width, out_width, dt, pack_frac_given, particle_diameter
# go('/data2/acarter/Spectral_Sophie_Boxes/data/spec_softetakt_long_run_dtau_0.025_nsave_2.suspension_phi_0.1_L_1280_modified.txt', 'brennan_hydro_010_L1280', 1280, 1280, 0.5, 0.1, 1)
# go(
#     infile = '/data2/acarter/Spectral_Sophie_Boxes/data/spec_softetakt_long_run_dtau_0.025_nsave_2.suspension_phi_0.1_L_544_modified.txt',
#     outfile = 'brennan_hydro_010_L544',
#     orig_width = 544, 
#     out_width = 544,
#     dt = 0.5,
#     pack_frac_given = 0.1,
#     particle_diameter = 2.79
# )
# go('/data2/acarter/Spectral_Sophie_Boxes/data/spec_softetakt_long_run_dtau_0.025_nsave_2.suspension_phi_0.02_L_1600_modified.txt', 'brennan_hydro_002_L1600', orig_width=1600, out_width=1600, dt=0.5, pack_frac_given=0.02, particle_diameter=2.79)
# go('/data2/acarter/Spectral_Sophie_Boxes/data/spec_softetakt_long_run_dtau_0.025_nsave_2.suspension_phi_0.02_L_800_modified.txt',  'brennan_hydro_002_L800',   800,  800, 0.5, 0.02, 1)


"""
#################################
#################################
##########  Hydro Data ##########
#################################
#################################

#################################
# this is the longest and best run that I did for phi = 0.34 *with* hydro
# The naming is sort of wierd, but I used a timestep of 0.25(s) and saved every other step to a file,
# so the effective timestep (the time betwen data frames in the file) is dt=0.5(s).  
# it's quated in the file name as a dimensionless timestep dtau=0.025. This is a stupid way to name files and I
# stopped doing it later and just quoted dt
#################################
Lx = 320
Ly = 320
infile = "./data/spec_softetakt_long_run_dtau_0.025_nsave_2.suspension_phi_0.34_L_320_modified.txt"
outfile = "./Count_Data_Cpp/New_Py_Test_phi_0.34" #
Nframes = 26478 # number of data frames (this might be less than there actually are)
a = 1.395 #radius of particles


#################################
# this is the longest and best run that I did for phi = 0.66 *with* hydro
# The naming is sort of wierd, but I used a timestep of 0.125(s) and saved every 4th step to a file,
# so the effective timestep (the time betwen data frames in the file) is dt=0.5(s).  
# it's quated in the file name as a dimensionless timestep dtau=0.025 *which is wrong in this case and I copy/pasted wrong* 
# it should've been dtau=0.0125. but dt=0.5(s) is right. 
#################################
Lx = 288
Ly = 288
infile = "./data/spec_softetakt_long_run_dtau_0.025_nsave_4.suspension_phi_0.66_L_288_modified.txt"
outfile = "./Count_Data_Cpp/Py_Test_phi_0.66"
Nframes = 8828 # number of data frames


#################################
#################################
########  No Hydro Data #########
#################################
#################################

# this is no-hydro data with an effective dt = (sim. timestep)*(save frequency) = 0.0125*40
# with packing fraction phi=0.34 and Lx=Ly=640
# to see how many frames are in the file just look at the third column of the last row
# $ tail -1 noHydro2D_Leim_run_dt_0.0125_nsave_40.suspension_phi_0.34_L_640_modified.txt
# 445.3523336427782 390.71188213483447 40000
infile = noHydro2D_Leim_run_dt_0.0125_nsave_40.suspension_phi_0.34_L_640_modified.txt


# this is no-hydro data with an effective dt = (sim. timestep)*(save frequency) = 9.765625e-5*32
# with packing fraction phi=0.34 and Lx=Ly=640
infile = noHydro2D_Leim_run_dt_9.765625e5_nsave_32.suspension_phi_0.34_L_640_modified.txt

# this is no-hydro data with an effective dt = (sim. timestep)*(save frequency) = 9.765625e-5*256
# with packing fraction phi=0.66 and Lx=Ly=640
noHydro2D_Leim_run_dt_9.765625e5_nsave_256_long.suspension_phi_0.66_L_640_eq_modified.txt
"""