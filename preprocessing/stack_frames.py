import common
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


# select your method here
# METHOD = DIFF_WITH_PREVIOUS
# METHOD = DIFF_WITH_ZERO
# METHODS = [NONE, DIFF_WITH_PREVIOUS]
# METHODS = [DIFF_WITH_PREVIOUS]
# METHODS = [NONE, DIFF_WITH_TEN,REMOVE_BACKGROUND, DIFF_WITH_ZERO, DIFF_WITH_PREVIOUS]
# METHODS = [DIFF_WITH_TEN]
METHODS = [preprocessing.stack_movie.NONE]
# METHODS = [NONE, DIFF_WITH_PREVIOUS, REMOVE_BACKGROUND]
# METHODS = [DIFF_WITH_PREVIOUS]
# METHODS = [DIFF_WITH_ZERO]
METHODS = [preprocessing.stack_movie.REMOVE_BACKGROUND]

ADD_DRIFT = False

        
if __name__ == '__main__':
    for file in common.files_from_argv('preprocessing/data/', 'stack_'):
        
        data = common.load(f'preprocessing/data/stack_{file}.npz')
        
        # try:
        if True:
            for METHOD in METHODS:
                for crop in [False]:
                # for crop in [False, True]:
            
                    filename = f'preprocessing/figures_png/stack_movie_{file}'

                    if crop:
                        filename += '_crop'

                    if METHOD == preprocessing.stack_movie.REMOVE_BACKGROUND:
                        filename += '_bkgrem'
                    elif METHOD == preprocessing.stack_movie.DIFF_WITH_ZERO:
                        filename += '_diffframe1'
                    elif METHOD == preprocessing.stack_movie.DIFF_WITH_PREVIOUS:
                        filename += '_diffprev'
                    elif METHOD == preprocessing.stack_movie.DIFF_WITH_TEN:
                        filename += '_diff10'
                    if HIGHLIGHTS:
                        filename += '_highlights'
                    if BACKWARDS:
                        filename += '_backwards'
                    if ADD_DRIFT:
                        filename += '_drifted'
                    
                    filename += '.gif'

                    preprocessing.stack_movie.go(
                        file,
                        data,
                        outputfilename=filename,
                        add_drift=ADD_DRIFT,
                        method=METHOD,
                        display_small=crop,
                        highlights=HIGHLIGHTS,
                        backward=BACKWARDS,
                        output_type='frames',
                        figsize_mult = 1/3,
                    )
        # except Exception as err:
        #     print()
        #     print(err)
        #     print()