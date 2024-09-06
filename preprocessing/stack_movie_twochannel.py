# import matplotlib.pyplot as plt
# from matplotlib.colors import LinearSegmentedColormap
# import numpy as np

# colors = [(0, 0, 0), (1, 0, 0)] # first color is black, last is red
# cm = LinearSegmentedColormap.from_list(
#         "Custom", colors, N=20)
# mat = np.indices((10,10))[1]
# plt.imshow(mat, cmap=cm)
# plt.show()


import numpy as np
import matplotlib.pyplot as plt
import common
import sys
import matplotlib.cm
import warnings
from preprocessing.stack_movie import save_array_movie, TWOCHANNEL

def go(file):
    data1 = common.load(f'preprocessing/data/stack_{file}r.npz')
    stack1 = data1['stack']
    pixel_size1 = data1['pixel_size']
    time_step1 = data1['time_step']

    data2 = common.load(f'preprocessing/data/stack_{file}g.npz')
    stack2 = data2['stack']
    pixel_size2 = data2['pixel_size']
    time_step2 = data2['time_step']

    assert np.all(stack1.shape == stack2.shape)
    assert pixel_size1 == pixel_size2
    assert time_step1 == time_step2


    filename = f'stack_movie_{file}'
    
    save_array_movie(stack=None, pixel_size=pixel_size1, time_step=time_step1, file=file,
                        outputfilename=f"preprocessing/figures_png/{filename}.gif",
                        stacks=[stack1, stack2], stackcolors=['red', 'green'],
                        method=TWOCHANNEL)
    # save_array_movie(stack, pixel_size, time_step, file, f"/home/acarter/presentations/cmd31/figures/{filename}.mp4",
    #                  nth_frame=data.get('nth_frame', 1), show_func=show_twochannel_image)
    # save_array_movie(stack_copy, pixel_size, time_step, file, f"/home/acarter/presentations/cin_first/figures/{filename}.mp4")



if __name__ == '__main__':
    go('marine2_C0')
    go('marine2_C1')
    go('marine2_C2')
    go('marine2_C3')
    go('marine2_C4')
    go('marine2_D0')
    go('marine2_D1')
    go('marine2_E0')
    go('marine2_E1')
    go('marine2_E2')
    go('marine2_E3')