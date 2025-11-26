import common
import preprocessing.stack_movie

path = '/home/acarter/presentations/soil_25_nov/figures'

preprocessing.stack_movie.go(
    'faxtor030b_small',
    outputfilename=f'{path}/faxtor030b_small_highlights',
    # method='frames',
    # display_small=crop,
    highlights=True,
    output_type='frames',
    figsize_mult = 1/2,
    annotation_color='black',
)
preprocessing.stack_movie.go(
    'faxtor006a',
    outputfilename=f'{path}/faxtor006a_highlights',
    highlights=True,
    output_type='frames',
    # display_small=True,
    figsize_mult = 1/2,
    annotation_color='black',
)
preprocessing.stack_movie.go(
    'faxtor006a',
    outputfilename=f'{path}/faxtor006a',
    output_type='frames',
    figsize_mult = 1/2,
    annotation_color='black',
    # display_small=True,
)
# preprocessing.stack_movie.go(
#     'faxtor017b_small',
#     outputfilename=f'{path}/faxtor017b_small_highlights',
#     highlights=True,
#     output_type='frames',
# )