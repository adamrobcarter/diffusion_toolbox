import common
import tifffile
import numpy as np


if __name__ == '__main__':
    file = 'raw_data/carlos/carlos01.tif'
    stack = common.load_tif(file)

    time_step = 1
    max_time_hours = round(stack.shape[0]*time_step/60/60, 2)
    pixel_size = 0.25935

    common.save_data('preprocessing/data/stack_carlos01.npz',
        stack=stack.astype(np.float32),
        time_step=time_step, pixel_size=pixel_size, particle_diameter=2.5, # 3.0 is my guess
        window_size_x=stack.shape[1]*pixel_size, window_size_y=stack.shape[2]*pixel_size,
        max_time_hours=max_time_hours,
        source_file=file,
    )