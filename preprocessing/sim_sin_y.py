import common
import numpy as np

if __name__ == '__main__':
    x = np.arange(0, 1024)
    y = np.arange(0, 1024)

    xx, yy = np.meshgrid(x, y, indexing='ij')

    f = 0.02

    frame = np.sin(2 * np.pi * xx * f)
    frame += 0.01 * np.random.normal(size=frame.shape)
    frame += 2

    stack = np.zeros((1, *frame.shape))
    stack[0] = frame

    common.save_data(f'preprocessing/data/stack_sin_y_{f}.npz', stack=stack, pixel_size=1.0, time_step=1.0)