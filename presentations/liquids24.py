path = '/home/acarter/presentations/liquids24/figures/'

import visualisation.Ds_overlapped
# visualisation.Ds_overlapped.go('eleanorlong034',
#     # ['MSD_short', 'D0Sk', 'f', 'f_short',],
#     ['MSD_short', 'D0Sk_theory', 'f_short',],
#     PLOT_AGAINST_K=True, TWO_PI=True, logarithmic_y=True,
#     show_window=False, show_pixel=False, show_pack_frac_plateau=True,
#     ylim=(0.7, 7),
#     export_destination=path+'Ds_f.pdf',
# )
# visualisation.Ds_overlapped.go('eleanor0.01',
#     # ['MSD_short', 'D0Sk', 'f', 'f_short',],
#     ['MSD_short', 'D0Sk_theory', 'f_short',],
#     PLOT_AGAINST_K=True, TWO_PI=True, logarithmic_y=True,
#     show_window=False, show_pixel=False, show_pack_frac_plateau=True,
#     export_destination=path+'Ds_f_lowdensity.pdf', PRESENT_SMALL=True,
# )
# visualisation.Ds_overlapped.go('eleanorlong034',
#     ['MSD_short', 'D_of_L_theory', 'boxcounting_shorttime', 'boxcounting_collective'],#, 'timescaleint_nofit', 'timescaleint'],
#     PLOT_AGAINST_K=False, TWO_PI=True, logarithmic_y=True,
#     show_window=False, show_pixel=False, show_pack_frac_plateau=True,
#     label_pack_frac=True, ylim=(0.5, 36),
#     export_destination=path+'Ds_boxcounting.pdf',
# )
# visualisation.Ds_overlapped.go('eleanor0.01',
#     ['MSD_short', 'boxcounting_shorttime', 'boxcounting_collective', 'D_of_L_theory'],#, 'timescaleint_nofit', 'timescaleint'],
#     PLOT_AGAINST_K=False, TWO_PI=True, logarithmic_y=True,
#     show_window=False, show_pixel=False, show_pack_frac_plateau=True,
#     export_destination=path+'Ds_boxcounting_lowdensity.pdf', PRESENT_SMALL=True,
# )
# visualisation.Ds_overlapped.go('eleanorlong034',
#     ['MSD_short', 'boxcounting_collective', 'f_short'],#, 'timescaleint_nofit', 'timescaleint'],
#     PLOT_AGAINST_K=False, TWO_PI=True, logarithmic_y=True, ylim=(0.5, 36),
#     show_window=False, show_pixel=False, show_pack_frac_plateau=False, label_k_scaling=True,
#     label_pack_frac=False,
#     export_destination=path+'Ds_boxcounting_and_f.pdf',
# )
visualisation.Ds_overlapped.go('eleanorlong034',
    ['MSD_short', 'boxcounting_collective', 'f_short', 'boxcounting_shorttime'],#, 'timescaleint_nofit', 'timescaleint'],
    plot_against_k=False, TWO_PI=True, logarithmic_y=True, ylim=(0.5, 36),
    show_window=False, show_pixel=False, show_pack_frac_plateau=False, label_k_scaling=True,
    label_pack_frac=False,
    export_destination=path+'Ds_boxcounting_and_f.png',
)
# visualisation.Ds_overlapped.go('eleanorlong034',
#     ['MSD_short', 'boxcounting_collective', 'f_short'],#, 'timescaleint_nofit', 'timescaleint'],
#     PLOT_AGAINST_K=False, TWO_PI=False, logarithmic_y=True, ylim=(0.5, 36),
#     show_window=False, show_pixel=False, show_pack_frac_plateau=False, label_k_scaling=True,
#     export_destination=path+'Ds_boxcounting_and_f_oneoverk.pdf', PRESENT_SMALL=True,
# )

# import particle_linking.show_static
# particle_linking.show_static.go('eleanor0.01', path+'linked.png', path+'positions.png', path+'intensities.png')


# import scattering_functions.show_Fs_overlayed
# scattering_functions.show_Fs_overlayed.go('eleanor0.01', SHOW_FIT=True, export_destination=path+'fkt_decay.pdf')

# import box_counting.example
# box_counting.example.go('eleanor0.34', 6.6, 10, path+'boxcounting_example_small.pdf', label_boxes=False)
# box_counting.example.go('eleanor0.34', 38, 10, path+'boxcounting_example_big.pdf', label_boxes=False)


# import visualisation.Ds_overlapped_mult
# visualisation.Ds_overlapped_mult.go(['eleanorlong010', 'eleanorlong034', 'eleanorlong066'], export_destination=path+'Ds_mult.pdf')