import countoscope_old as countoscope
import numpy as np
import time
import sys
import common

def calc_and_save(file, box_sizes_x, box_sizes_y, sep_sizes_x, sep_sizes_y, data, output_file_name,
                  save_counts=False, extra_to_save={},
                  save_data=True, skip_processing=False):
    t0 = time.time()

    if not save_data:
        print('WARNING: I am not going to save the data')

    particles                = data['particles']
    time_step                = data['time_step']
    particle_diameter        = data.get('particle_diameter', np.nan)
    particle_diameter_calced = data.get('particle_diameter_calced')
    window_size_x            = data.get('window_size_x')
    window_size_y            = data.get('window_size_y')
    # assert 'density' in data

    time_values = np.unique(particles[:, 2])
    assert time_values[0] == 0
    assert time_values[1] == 1

    results = countoscope.calculate_N1N2(
        data=particles,
        window_size_x=window_size_x, window_size_y=window_size_y, 
        box_sizes_x=box_sizes_x, box_sizes_y=box_sizes_y,
        sep_sizes_x=sep_sizes_x, sep_sizes_y=sep_sizes_y,
        # skip_processing=skip_processing,
        return_counts=save_counts,
        # return_histograms=False,
    )

    output = dict(
        filename=output_file_name,
        box_sizes_x=box_sizes_x, box_sizes_y=box_sizes_y,
        sep_sizes_x=sep_sizes_x, sep_sizes_y=sep_sizes_y,
        time_step=time_step, particle_diameter=particle_diameter,
        particle_diameter_calced=particle_diameter_calced, computation_time=time.time()-t0,
        pack_frac=data.get('pack_frac'), density=data.get('density', common.calc_density(particles, window_size_x, window_size_y)),
        pack_frac_given=data.get('pack_frac_given'), max_time_hours=data.get('max_time_hours'),
        window_size_x=window_size_x, window_size_y=window_size_y, pixel_size=data.get('pixel_size'),
        **results._asdict(),
        **extra_to_save
    )

    if save_data:
        common.save_data(**output)
    else:
        print('not saving data')

    return output


def go(file):
    data = common.load(f'particle_detection/data/particles_{file}.npz')

    output_filename = f'box_counting/data/pnv_{file}'
        # output_filename += '_moreoverlap'
    output_filename += '.npz'

    # box_sizes_x = np.array([2, 5, 10, 20])
    # box_sizes_y = np.array([2, 5, 10, 20])
    # spacing_x = np.array([2, 2, 2, 2])
    # spacing_y = np.array([2, 2, 2, 2])
    box_sizes_x = np.logspace(np.log10(0.5), np.log10(20), num=10)[:9]
    box_sizes_y = box_sizes_x
    spacing_x = np.full_like(box_sizes_x, 4)
    spacing_y = spacing_x
    sep_sizes_x = spacing_x - box_sizes_x
    sep_sizes_y = spacing_y - box_sizes_y

    t0 = time.time()
    calc_and_save(
        file=file, data=data,
        box_sizes_x=box_sizes_x, box_sizes_y=box_sizes_y,
        sep_sizes_x=sep_sizes_x, sep_sizes_y=sep_sizes_y,
        output_file_name=output_filename, save_counts=False,
        save_data=True, extra_to_save=dict(
            v_profile = data['v_profile'],
        ),
    )
    print(f'took {time.time()-t0:.0f}s')

if __name__ == '__main__':
    for file in common.files_from_argv('particle_detection/data', 'particles_'):
        go(file)