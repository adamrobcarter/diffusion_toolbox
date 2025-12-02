import numpy as np
import cv2
import matplotlib.pyplot as plt
import tqdm
import common

if __name__ == '__main__':
    # read video
    video = cv2.VideoCapture('raw_data/Pilier seul.avi')

    num_frames = int(video.get(cv2.CAP_PROP_FRAME_COUNT))
    window_size_x = int(video.get(cv2.CAP_PROP_FRAME_WIDTH))
    window_size_y = int(video.get(cv2.CAP_PROP_FRAME_HEIGHT))
    frame_rate = video.get(cv2.CAP_PROP_FPS)

    stack = np.full((num_frames, window_size_y, window_size_x), np.nan)
    print('stack size', common.arraysize(stack))

    for frame_i in tqdm.trange(num_frames):
        success, frame = video.read()
        assert success
        frame_bw = frame.mean(axis=2) # convert to grayscale by averaging RGB channels
        stack[frame_i, :, :] = frame_bw

    stack = np.swapaxes(stack, 1, 2)
    stack = stack[:, :, ::-1]

    print('histogram: stack')
    common.term_hist(stack)

    # stack[stack < 100] = stack.mean()

    common.save_data('preprocessing/data/stack_sophie1.npz',
        stack=stack,
        window_size_x=window_size_x,
        window_size_y=window_size_y,
        time_step=1/frame_rate,
        pixel_size=1
    )