from countoscope import Countoscope as CountoscopeNew
import countoscope_old
import common
import numpy as np

for file in common.files_from_argv('particle_detection/data', 'particles_'):
    data = common.load(f'particle_detection/data/particles_{file}.npz')
    time_step                = data['time_step']
    pixel_size               = data.get('pixel_size', 1)
    num_timesteps            = data['num_timesteps']
    particle_diameter        = data.get('particle_diameter', np.nan)
    particle_diameter_calced = data.get('particle_diameter_calced')
    depth_of_field           = data.get('depth_of_field')
    window_size_x            = data.get('window_size_x')
    window_size_y            = data.get('window_size_y')
    particles                = data['particles']

    box_size = 8

    results_old = countoscope_old.calculate_nmsd(data=particles,
                                        #  window_size_x=window_size_x, window_size_y=window_size_y, 
                                         box_sizes=[box_size],
                                         sep_sizes=[0],
                                         return_counts=True
    )
    old_nmsd = results_old.nmsd
    print(old_nmsd.shape)

    countoscope_new = CountoscopeNew(
        box_size          = np.array([box_size, box_size]),
        trajectory_array  = particles,
        trajectory_labels = ['x', 'y', 't'],
    )
    countoscope_new.count()
    new_nmsd = countoscope_new.evaluate_deltaN2()
    print(new_nmsd.shape)

    common.save_data(f'box_counting/data/new_vs_old_{file}.npz', new=new_nmsd, old=old_nmsd)