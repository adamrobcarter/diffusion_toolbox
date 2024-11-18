import preprocessing.stack_movie
import visualisation.Dc_literature
import visualisation.Ds_overlapped
import visualisation.Ds_overlapped_mult
import box_counting.D_of_L
import box_counting.msd_single

path = 'tmp'

########################## fig 3 #########################
fig3_size = (8/3, 2)
# box_counting.D_of_L.go(
#     'eleanorlong001',
#     'var',
#     export_destination=f'{path}/fig3a.pdf',
#     figsize=fig3_size,
#     save_data=False,
#     title='',
#     show_nofit_cutoff=False
# )
# visualisation.Ds_overlapped_mult.go(
#     ['eleanorlong001', 'eleanorlong010'],
#     sources=['MSD_first', 'timescaleint', 'timescaleint_nofit'],
#     export_destination=f'{path}/fig3b.pdf',
#     figsize=fig3_size,
# )
# visualisation.Ds_overlapped_mult.go(
#     ['eleanorlong010', 'sim_nohydro_010_L544', 'brennan_hydro_010_L544'],
#     sources=['MSD_first', 'timescaleint', 'timescaleint_nofit'],
#     export_destination=f'{path}/fig3c.pdf',
#     figsize=fig3_size,
# )


########################## fig 5 #########################
visualisation.Ds_overlapped_mult.go(
    ['eleanorlong001_crop1.0', 'eleanorlong001_crop0.5', 'eleanorlong001_crop0.25', 'eleanorlong001_crop0.125', 'eleanorlong001_crop0.0625'],
    sources=['MSD_first', 'f_first_first'],
    export_destination=f'{path}/fig5d.pdf',
    figsize=fig3_size,
    plot_against_k=True
)


########################## fig 6 #########################
fig6_size = (8/4, 2)
box_counting.msd_single.go(
    
)

visualisation.Ds_overlapped.go(
    'eleanorlong010',
    sources=['MSD_first', 'timescaleint', 'timescaleint_replacement'],
    export_destination=f'{path}/fig6b.pdf',
    figsize=fig6_size,
)