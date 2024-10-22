import preprocessing.stack_movie
import visualisation.Dc_literature
import visualisation.Ds_overlapped

path = '/home/acarter/presentations/followupcommittee/figures/'

figsize_Ds_overlapped_big = (3.5, 3)

visualisation.Ds_overlapped.go(
    'eleanorlong034',
    ['boxcounting_collective', 'MSD_first', 'f_short'],
    export_destination=path+'Ds_overlapped_eleanorlong034_boxcounting_MSD_f.pdf',
    figsize=figsize_Ds_overlapped_big,
)
visualisation.Ds_overlapped.go(
    'eleanorlong034',
    ['boxcounting_collective'],
    export_destination=path+'Ds_overlapped_eleanorlong034_boxcounting_collective.pdf',
    figsize=figsize_Ds_overlapped_big,
)
visualisation.Ds_overlapped.go(
    'eleanorlong034',
    ['boxcounting_collective', 'MSD_first'],
    export_destination=path+'Ds_overlapped_eleanorlong034_boxcounting_collective_MSD.pdf',
    figsize=figsize_Ds_overlapped_big,
)

visualisation.Ds_overlapped.go(
    'eleanorlong034',
    ['f_short'],
    export_destination=path+'Ds_overlapped_eleanorlong034_f_short.pdf',
    figsize=figsize_Ds_overlapped_big,
    PLOT_AGAINST_K=True,
)
visualisation.Ds_overlapped.go(
    'eleanorlong034',
    ['f_short', 'MSD_first'],
    export_destination=path+'Ds_overlapped_eleanorlong034_f_short_MSD.pdf',
    figsize=figsize_Ds_overlapped_big,
    PLOT_AGAINST_K=True,
)

preprocessing.stack_movie.go(
    'pierre_exp',
    path+'pierre_exp.mp4',
    display_small=True,
    flip_y=True,
)

visualisation.Dc_literature.go(path+'Dc_literature.pdf')