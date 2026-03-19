import common
import numpy as np

if __name__ == '__main__':
    file = 'raw_data/carlos/Project_Test1_ch00.tif'
    stack = common.load_tif(file)

    time_step = 0.2
    max_time_hours = round(stack.shape[0]*time_step/60/60, 2)
    pixel_size = 0.324355
    particle_diameter = None

    common.save_data('preprocessing/data/stack_carlos03.npz',
        stack=stack.astype(np.float32),
        time_step=time_step, pixel_size=pixel_size, particle_diameter=particle_diameter,
        window_size_x=stack.shape[1]*pixel_size, window_size_y=stack.shape[2]*pixel_size,
        max_time_hours=max_time_hours,
        source_file=file,
    )