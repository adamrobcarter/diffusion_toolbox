import DDM.phiDM
import common
import numpy as np
import time

background_removed = True

if __name__ == '__main__':
    for file in common.files_from_argv('preprocessing/data', 'stack_'):
        data = common.load(f'preprocessing/data/stack_{file}.npz')
        stack      = data['stack']
        pixel_size = data['pixel_size']
        time_step  = data['time_step']

        print(stack.shape)
        stack = common.add_drift_intensity(stack, 2)
        print('added drift')
        print(stack.shape)

        t0 = time.time()
        Rs = DDM.phiDM.calc(stack, pixel_size)
        t1 = time.time()

        common.save_data(f'DDM/data/phiDM_{file}', R=Rs)