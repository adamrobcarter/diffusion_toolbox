import scattering_functions.calc_both
import scattering_functions.scattering_functions_nonumba
import common
import pickle
import time
import numpy as np

d_frames = [1, 2, 3, 4, 5, 6, 7, 8]

if __name__ == '__main__':
    for file in common.files_from_argv('particle_detection/data', 'particles_'):

        particles_at_frame, num_timesteps, d_frames, min_K, max_K, data = scattering_functions.calc_both.setup('F', file, d_frames)

        all_data = []

        for max_time_origins in np.logspace(1, 5, num=25):
        # for max_time_origins in [1, 8, 64, 512]:
        # for max_time_origins in [1, 64]:
            print('doing', max_time_origins)

            t0 = time.time()

            Fs, F_unc, ks, F_unbinned, F_unc_unbinned, k_unbinned, k_x, k_y = scattering_functions.scattering_functions_nonumba.intermediate_scattering(
                False, 'F', 50, max_time_origins, d_frames,
                particles_at_frame, num_timesteps, max_K, min_K, cores=16, use_zero=False, use_big_k=False, linear_log_crossover_k=0.2
            )

            t1 = time.time()

            ret = dict(F=Fs, F_unc=F_unc, k=ks, k_x=k_x, k_y=k_y,
                F_unbinned=F_unbinned, k_unbinned=k_unbinned, F_unc_unbinned=F_unc_unbinned,
                t=d_frames*data['time_step'], min_K=min_K,
                max_time_origins=max_time_origins, computation_time=t1-t0,
                particle_diameter=data['particle_diameter'],
                pixel_size=data['pixel_size'], pack_frac_given=data.get('pack_frac_given'),
                window_size_x=data.get('window_size_x'), window_size_y=data.get('window_size_y'),
                NAME=data.get('NAME'), channel=data.get('channel'),
            )

            all_data.append(ret)

            print(f'done in {t1-t0:.0f}s')

            if t1-t0 > 200:
                break
        
        with open(f'isf/data/F_quantify_noise_{file}.pickle', "wb" ) as f:
            print('saving')
            pickle.dump(all_data, f)
            print('saved')
