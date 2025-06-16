import numpy as np
import matplotlib.pyplot as plt
import common
import sys
import matplotlib.cm
import warnings
import preprocessing.stack_movie

DISPLAY_SMALL = False
INVERSE_COLORS = False
SHOW_TIMESTEP = True
HIDE_ANNOTATIONS = False

SHOW_HISTOGRAMS = False

HIGHLIGHTS = False # displays 50 frames evenly throughout the stack instead of the first 50
# HIGHLIGHTS = True
# BACKWARDS = True
BACKWARDS = False # plays the stack backwards. DIFF_WITH_ZERO is now compared to the last frame

# possible display methods
REMOVE_BACKGROUND = 1
DIFF_WITH_ZERO = 2 # frame - frame1
DIFF_WITH_PREVIOUS = 3
DIFF_WITH_TEN = 6
NONE = 4
TWOCHANNEL = 5

METHOD_NAME = {
    REMOVE_BACKGROUND: 'i(t) - mean(i)',
    DIFF_WITH_ZERO: 'i(t) - i(0)',
    DIFF_WITH_PREVIOUS: 'i(t) - i(t-1)',
    DIFF_WITH_TEN: 'i(t) - i(t-10)',
    NONE: '',
    TWOCHANNEL: '',
}

# select your method here
# METHOD = DIFF_WITH_PREVIOUS
# METHOD = DIFF_WITH_ZERO
# METHODS = [NONE, DIFF_WITH_PREVIOUS]
# METHODS = [DIFF_WITH_PREVIOUS]
# METHODS = [NONE, DIFF_WITH_TEN,REMOVE_BACKGROUND, DIFF_WITH_ZERO, DIFF_WITH_PREVIOUS]
# METHODS = [DIFF_WITH_TEN]
# METHODS = [NONE]
METHODS = [NONE, REMOVE_BACKGROUND]
# METHODS = [DIFF_WITH_ZERO]
# METHOD = REMOVE_BACKGROUND

ADD_DRIFT = False

if __name__ == '__main__':
    for file in common.files_from_argv('preprocessing/data/', 'stack_'):
        # num = int(file[6:9])
        # if num < 114:
        #     continue
        
        data = common.load(f'preprocessing/data/stack_{file}.npz')
        
        # try:
        if True:
            for METHOD in METHODS:
                for crop in [False]:
                # for crop in [False, True]:
            
                    filename = f'preprocessing/figures_png/stack_movie_{file}'

                    if crop:
                        filename += '_crop'

                    if METHOD == REMOVE_BACKGROUND:
                        filename += '_bkgrem'
                    elif METHOD == DIFF_WITH_ZERO:
                        filename += '_diffframe1'
                    elif METHOD == DIFF_WITH_PREVIOUS:
                        filename += '_diffprev'
                    elif METHOD == DIFF_WITH_TEN:
                        filename += '_diff10'
                    if HIGHLIGHTS:
                        filename += '_highlights'
                    if BACKWARDS:
                        filename += '_backwards'
                    if ADD_DRIFT:
                        filename += '_drifted'
                    
                    filename += '.mp4'

                    preprocessing.stack_movie.go(
                        file,
                        data,
                        outputfilename=filename,
                        add_drift=ADD_DRIFT,
                        method=METHOD,
                        display_small=crop,
                        highlights=HIGHLIGHTS,
                        backward=BACKWARDS
                    )
        # except Exception as err:
        #     print()
        #     print(err)
        #     print()