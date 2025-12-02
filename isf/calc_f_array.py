import common
import isf.calc_both
import warnings
import scattering_functions
import time
import numpy as np
import tqdm

max_time_origins = 150
num_k_bins = 60
d_frames = [0, 32]

if __name__ == '__main__':
    for file in common.files_from_argv('particle_detection/data', 'particles_array'):
        data = common.load(f'particle_detection/data/particles_array{file}.npz')

        particles_array = data['particles_array']
        num_repeats = particles_array.shape[0]

        Fs = []
        F_uncs = []
        ks = None

        t0 = time.time()

        for i in tqdm.trange(num_repeats, desc='array'):
            particles = particles_array[i, :, :]

            warnings.warn('is d_frames now d_times?')

            particles_at_frame, times, d_frames, min_K, max_K, data = isf.calc_both.setup(file, particles, data, d_frames)

            F, F_unc, ks, F_unbinned, F_unc_unbinned, k_unbinned, k_x, k_y, d_frames = scattering_functions.intermediate_scattering(
                'F', num_k_bins, max_time_origins, d_frames, 
                particles_at_frame, times, max_K, min_K, cores=16,
            )

            Fs.append(F)
            F_uncs.append(F)

        t1 = time.time()
            
        filename = f"isf/data/F_array{file}"

        to_save = dict(F_array=np.stack(Fs), F_unc_array=np.stack(F_uncs), k=ks, 
            t=d_frames*data['time_step'], min_K=min_K,
            num_k_bins=num_k_bins, max_time_origins=max_time_origins, computation_time=t1-t0,
            particle_diameter=data['particle_diameter'],
            pixel_size=data.get('pixel_size'), pack_frac_given=data.get('pack_frac_given'), pack_frac=data.get('pack_frac'),
            window_size_x=data.get('window_size_x'), window_size_y=data.get('window_size_y'),
            NAME=data.get('NAME'), channel=data.get('channel'),
            num_timesteps=times, max_time_hours=data.get('max_time_hours'), density=data.get('density'),
        )

        common.save_data(filename, **to_save)
