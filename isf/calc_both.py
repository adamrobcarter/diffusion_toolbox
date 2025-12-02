import numpy as np
# import scattering_functions.scattering_functions_nonumba as scattering_functions
import scattering_functions
import common
import time
import warnings
import isf.show_both

def less_than_or_close(a, b):
    return a <= b or np.isclose(a, b)
def greater_than_or_close(a, b):
    return a >= b or np.isclose(a, b)

def load_and_setup(file, F_type='F', d_frames=None, drift_removed=False, max_K=None, quiet=False):
    if F_type == 'F_s':
        filepath = f"particle_linking/data/trajs_{file}.npz"
    else:
        filepath = f"particle_detection/data/particles_{file}.npz"
    data = common.load(filepath, quiet=quiet)

    return setup(file, data['particles'], data, d_frames, F_type, drift_removed, max_K, quiet) # passing particles and metadata separately might seem weird but we need it when calculating on arrays
    
def setup(file, particles, metadata, d_frames, F_type='F', drift_removed=False, max_K=None, quiet=False):
    particles     = particles.astype(np.float32) # rows of x,y,t
    pixel_size    = metadata.get('pixel_size')
    window_size_x = metadata['window_size_x']
    window_size_y = metadata['window_size_y']
    dimension     = metadata.get('dimension', 2)
    # assert 'density' in data
    if '_pot' not in file and 'nbody' not in file:
        assert particles[:, 0].min() >= 0
        assert particles[:, 1].min() >= 0
        assert less_than_or_close(particles[:, 0].max(), window_size_x), f'particles[:, 0].max() = {particles[:, 0].max()}, window_size_x={window_size_x}'
        assert less_than_or_close(particles[:, 1].max(), window_size_y), f'particles[:, 1].max() = {particles[:, 1].max()}, window_size_y={window_size_y}'
    
    times = int(particles[:, 2].max() + 1)

    if drift_removed:
        particles = common.remove_drift(particles)
    else:
        if not quiet: print('not removing drift')

    if False:
        print('adding drift')
        particles = common.add_drift(particles, 0.05, 0.05)
    
    min_K = 2*np.pi / np.array(common.get_used_window(file, window_size_x, window_size_y))
    if not quiet: print(f'min k =', min_K)

    if not max_K:
        if pixel_size:
            max_K = 2 * np.pi / pixel_size
        else:
            if not quiet: print('Pixel size not given, using max_K = 21.81661564992912')
            max_K = 21.81661564992912 # was 10
            # we use this cause it's the same as used by eleanorlong

    if dimension == 2:
        columns = {'x': 0, 'y': 1, 't': 2}
    elif dimension == 3:
        columns = {'x': 0, 'y': 1, 't': 3}

    particles_at_frame, times = scattering_functions.get_particles_at_frame(F_type, particles, columns=columns)

    if d_frames is None:
        d_frames = times[:10]
        # d_frames = common.exponential_integers(1, min(isf.show_both.LONG_END, times-1)) - 1
    else:
       d_frames = np.array(d_frames) # user provided

    return particles_at_frame, times, d_frames, min_K, max_K, metadata

def calc_for_f_type(
        file,
        F_type,
        log=True,
        max_K=None,
        drift_removed=False,
        num_k_bins=50, # computation time is proportional to this squared
        file_suffix='',
        file_prefix='',
        cores=16,
        max_time_origins=100, # this sets (approx) how many different time origins will be averaged over
                              # computation time is directly proportional
        use_zero=True,
        save=True,
        d_frames=None,
        use_doublesided_k=False,
        window=None,
        quiet=False,
    ):

    t0 = time.time()

    warnings.warn('is d_frames now d_times?')

    particles_at_frame, times, d_frames, min_K, max_K, data = load_and_setup(file, F_type, d_frames, drift_removed, max_K, quiet=quiet)

    if not quiet: print('starting calculation')
    if not quiet: print('particles_at_frame:', common.arraysize(particles_at_frame))

    if not quiet: print(f'going with min k = ({min_K[0]:.4f}, {min_K[1]:.4f})')
    if not quiet: print('going with d_frames =', d_frames)

    Fs, F_unc, ks, F_unbinned, F_unc_unbinned, k_unbinned, k_x, k_y, d_frames = scattering_functions.intermediate_scattering(
        F_type, num_k_bins, max_time_origins, d_frames, 
        particles_at_frame, times, max_K, min_K, cores=cores,
        use_doublesided_k=use_doublesided_k, window=window,
        Lx=data['window_size_x'], Ly=data['window_size_y'],
        quiet=quiet,
    )

    t1 = time.time()
        
    filename = f"isf/data/{F_type}_{file_prefix}{file}{file_suffix}"

    assert Fs.shape[0] == d_frames.size

    to_save = dict(F=Fs, F_unc=F_unc, k=ks, k_x=k_x, k_y=k_y,
        F_unbinned=F_unbinned, k_unbinned=k_unbinned, F_unc_unbinned=F_unc_unbinned,
        t=d_frames*data['time_step'], min_K=min_K,
        num_k_bins=num_k_bins, max_time_origins=max_time_origins, computation_time=t1-t0, log=log,
        particle_diameter=data['particle_diameter'], drift_removed=drift_removed,
        pixel_size=data.get('pixel_size'), pack_frac_given=data.get('pack_frac_given'), pack_frac=data.get('pack_frac'),
        window_size_x=data.get('window_size_x'), window_size_y=data.get('window_size_y'),
        NAME=data.get('NAME'), channel=data.get('channel'), window=window,
        num_timesteps=times, max_time_hours=data.get('max_time_hours'), density=data.get('density'),
    )

    if save:
        common.save_data(filename, quiet=quiet, **to_save)

    if not quiet: print()

    return to_save

if __name__ == '__main__':
    for file in common.files_from_argv('particle_detection/data', 'particles_'):
        calc_for_f_type(file, 'F_s')
        calc_for_f_type(file, 'f')
