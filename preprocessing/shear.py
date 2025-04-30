import numpy as np
import common
import time, os
import tqdm
import subprocess


def go(infile, outfile, particle_diameter, dt=None, nth_timestep=1, max_time=None,
       extra_source_file=None, skiprows=None, multiply_time_by=None, **kwargs):
    print(f'loading raw file, last modified {common.get_last_modified_time(infile)} ago')

    if max_time:
        raise NotImplemented()
        max_frame = max_time * dt
        max_rows = int((max_frame + 1) * expected_particles_per_frame)
        print('max_rows', max_rows)

        max_items = max_rows * 3
    else:
        max_rows = None
        max_items = -1

    t0 = time.time()

    data = np.fromfile(infile, dtype=np.float32, count=max_items) # there is a parameter to this function for not loading the whole array
    data = data.reshape((-1, 4))

    t1 = time.time()
    print(data.shape, data.dtype,  f'loaded {common.format_bytes(data.nbytes)} in {t1-t0:.1f}s at {common.format_bytes(data.nbytes/(t1-t0))}/sec')

    # # dt = None, times = [0, 0.5, 1, 1.5, ...]
    # times = np.unique(data[:, 2])
    # multiply_time_by = 1/times[1]
    # dt = times[1]
    # # dt = 0.5, times = [0, 1, 2, 3, ...]
    
    # discard the non-nth step rows
    if nth_timestep > 1:
        print('doing nth-time')
        data = data[data[:, 2] % nth_timestep == 0, :]
        dt *= nth_timestep
        # we gotta readjust the timestep number now we discarded some timesteps
        data[:, 2] /= nth_timestep

    if skiprows:
        #### there is probably a way of doing this with loadtxt/fromfile that would be quicker
        print(f'doing skiprows ({skiprows})')
        print(data.shape)
        data = data[skiprows:, :]
        print(data.shape)
        # now need to reset the time column
        data[:, 2] -= data[:, 2].min()
        print()

    # we should delete the last timestep cause it may be incomplete
    print('trimming end')
    last_timestep = data[:, 2].max()
    data = data[data[:, 2] != last_timestep, :]
    
    print('xy min max', data[:, 0].min(), data[:, 0].max(), data[:, 1].min(), data[:, 1].max())
    
    # print(f'loaded in {t1-t0:.0f}s. shape', data.shape, common.arraysize(data))
    # num_timesteps = data[:, 2].max()+1
    times = np.unique(data[:, 2])
    num_timesteps = times.size
    print(times, num_timesteps)
    # print(f'{num_timesteps:.0f} timesteps, {data[:, 2].max()*dt/60/60:.1f} hours')


    data[:, 2] -= data[:, 2].min()
    assert data[:, 2].min() == 0

    common.save_data(f'particle_linking/data/trajs_{outfile}.npz',
        particles=data,
        time_step=dt, particle_diameter=particle_diameter,
        max_time_hours=round(last_timestep/60/60, 2),
        source_file=infile, extra_source_file=extra_source_file,
        **kwargs
    )

    end_timestep = data[:, 2].max() // 8
    data_small = data[data[:, 2] < end_timestep, :]
    common.save_data(f'particle_linking/data/trajs_{outfile}_div8.npz',
        particles=data_small,
        time_step=dt, particle_diameter=particle_diameter,
        max_time_hours=round(end_timestep/60/60, 2),
        source_file=infile, extra_source_file=extra_source_file,
        **kwargs
    )
        
    end_timestep = data[:, 2].max() // 64
    data_small = data_small[data_small[:, 2] < end_timestep, :]
    common.save_data(f'particle_linking/data/trajs_{outfile}_div64.npz',
        particles=data_small,
        time_step=dt, particle_diameter=particle_diameter,
        max_time_hours=round(end_timestep/60/60, 2),
        source_file=infile,extra_source_file=extra_source_file,
        **kwargs
    )
    print()
    

def go_mesu_shear(filepath, suffix='', skip_rsync=False, **kwargs):
    particle_diameter = 2

    filename = filepath.split('/')[-1]
    if not skip_rsync:
        rsync_command = ['rsync', f'cartera@login.mesu.sorbonne-universite.fr:{filepath}', f'raw_data/mesu/{filename}', '--progress']
        print('launching', ' '.join(rsync_command))
        subprocess.run(rsync_command, check=True)

    kwargs = dict()
    
    for param in ['T', 'theta', 'shear']:
        if param in filename:
            print(param, filename, filename.split(param))
            val = float(filename.split(param)[1].split('_')[0])
            kwargs[param] = val
    # hydro = filename_no_ext.split('_')[0]

    # phistr = f'{phi*100:.0f}'.zfill(3)

    go(
        f'raw_data/mesu/{filename}',
        f'sim_shear{suffix}',
        particle_diameter=particle_diameter,
        extra_source_file=filepath,
        **kwargs
    )

if __name__ == '__main__':
    # shear/gravity balance
    # go_mesu_shear('/store/cartera/shear/shear_test_nograv.bin', suffix='_test_nograv')
    # go_mesu_shear('/home/cartera/New_Rigid_Code/shear0_test.bin', suffix='_shear0')
    # go_mesu_shear('/store/cartera/shear/shear_shear0_nowall_nograv.bin', suffix='_shear0_nowall_nograv')
    # go_mesu_shear('/store/cartera/shear/shear_shear0_nowall_nograv_T10.bin', suffix='_shear0_nowall_nograv_T10')
    # go_mesu_shear('/store/cartera/shear/shear_shear0_wall_grav_T10.bin', suffix='_shear0_wall_grav_T10')
    # go_mesu_shear('/store/cartera/shear/shear1.0_dt0.001_T0_grav_theta0.6_wall.bin', suffix='_shear1.0_wall_grav_T0_theta0.6')
    # go_mesu_shear('/store/cartera/shear/shear1.0_dt0.001_T0_grav_theta0.61_wall.bin', suffix='_shear1.0_wall_grav_T0_theta0.61')
    # go_mesu_shear('/store/cartera/shear/shear1.0_dt0.001_T0_grav_theta0.62_wall.bin', suffix='_shear1.0_wall_grav_T0_theta0.62')
    # go_mesu_shear('/store/cartera/shear/shear1.0_dt0.001_T0_grav_theta0.63_wall.bin', suffix='_shear1.0_wall_grav_T0_theta0.63')
    # go_mesu_shear('/store/cartera/shear/shear1.0_dt0.001_T0_grav_theta0.64_wall.bin', suffix='_shear1.0_wall_grav_T0_theta0.64')
    # go_mesu_shear('/store/cartera/shear/shear1.0_dt0.001_T0_grav_theta0.65_wall.bin', suffix='_shear1.0_wall_grav_T0_theta0.65')
    # go_mesu_shear('/store/cartera/shear/shear1.0_dt0.001_T0_grav_theta0.66_wall.bin', suffix='_shear1.0_wall_grav_T0_theta0.66')
    
    go_mesu_shear('/store/cartera/shear/shear0.0_RHS.bin', suffix='0.0_RHS')
    
    

