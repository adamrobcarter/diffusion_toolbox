import common
import numpy as np

for directory in ['sim', 'simdownsampled', 'exp']:
    print(f'doing {directory}')
    tifs = common.get_directory_files(f'raw_data/pierre_{directory}', 'tif')

    data = [common.load_tif(tif) for tif in tifs]
    stack = np.stack(data, axis=0)
    print(stack.shape)
    if directory == 'exp':
        pixel_size = 0.6
        time_step = 0.015
        particle_diameter = 2.4
    elif directory == 'sim':
        pixel_size = 0.1
        time_step = 0.012
        particle_diameter = 0.6
    elif directory == 'simdownsampled':
        pixel_size = 0.4
        time_step = 0.012
        particle_diameter = 0.6
    else:
        raise Exception('need to provide pixel size')

    print('saving')
    np.savez(f'preprocessing/data/stack_pierre_{directory}.npz', stack=stack, 
             pixel_size=pixel_size, particle_diameter=particle_diameter, time_step=time_step)
    
    print()