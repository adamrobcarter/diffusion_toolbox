import countoscope_old as countoscope
import numpy as np
import time
import common
from box_counting.calc_pnv import calc_and_save

def go(file, frame_deltas=[0, 1], aspect=1):
    data = common.load(f'particle_detection/data/particles_{file}.npz')

    output_filename = f'box_counting/data/pnv_{file}_aspect{aspect}_special'
        # output_filename += '_moreoverlap'
    output_filename += '.npz'


    a = 0.5
    b = 1.25992105/2
    c = 1.58740105/2
    box_sizes_base = np.array([#a,b,c, 
                               2*a,2*b,2*c, 4*a,4*b,4*c, 8*a,8*b,8*c, 16*a,16*b,16*c, 32*a,32*b,32*c])
    box_sizes_x = np.full(3*box_sizes_base.size, np.nan)
    box_sizes_x[0::3] = box_sizes_base / 1.1
    box_sizes_x[1::3] = box_sizes_base
    box_sizes_x[2::3] = box_sizes_base * 1.1

    # do the aspect ratio to get y
    box_sizes_y = np.copy(box_sizes_x)
    box_sizes_x *= np.sqrt(aspect)
    box_sizes_y /= np.sqrt(aspect)

    spacing_x = np.full_like(box_sizes_x, 4)
    spacing_y = spacing_x
    sep_sizes_x = spacing_x - box_sizes_x
    sep_sizes_y = spacing_y - box_sizes_y

    t0 = time.time()
    calc_and_save(
        file=file, data=data, frame_deltas=frame_deltas,
        box_sizes_x=box_sizes_x, box_sizes_y=box_sizes_y,
        sep_sizes_x=sep_sizes_x, sep_sizes_y=sep_sizes_y,
        output_file_name=output_filename, save_counts=True,
        save_data=True, extra_to_save=dict(
            v_profile = data.get('v_profile'),
            velocity_multiplier = data.get('velocity_multiplier'),
        ),
    )
    print(f'took {time.time()-t0:.0f}s')

if __name__ == '__main__':
    for file in common.files_from_argv('particle_detection/data', 'particles_'):
        go(file)