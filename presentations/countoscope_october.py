import visualisation.Ds_overlapped
# import visualisation.Ds_overlapped
import scattering_functions.show_decaytime

path = '/home/acarter/presentations/countoscope_october/figures/'


visualisation.Ds_overlapped.go(
    'eleanorlong034',
    ['MSD_short', 'timescaleint_nofit_cropped'],
    export_destination=path+'Ds_dominiguez_counting1.pdf',
    figsize=(5, 4),
    show_pixel=False, show_window=False,
    hide_msd=True,
)

visualisation.Ds_overlapped.go(
    'eleanorlong034',
    ['MSD_short', 'timescaleint_nofit_cropped'],
    export_destination=path+'Ds_dominiguez_counting1.pdf',
    figsize=(5, 4),
    show_pixel=False, show_window=False,
    hide_msd=True,
)
visualisation.Ds_overlapped.go(
    'eleanorlong034',
    ['MSD_short', 'timescaleint_nofit_cropped', 'D_of_L_theory'],
    export_destination=path+'Ds_dominiguez_counting2.pdf',
    figsize=(5, 4),
    show_pixel=False, show_window=False,
    hide_msd=True,
)
visualisation.Ds_overlapped.go(
    'eleanorlong034',
    ['MSD_short', 'timescaleint_nofit_cropped', 'D_of_L_theory', 'D_of_L_theory_Lh'],
    export_destination=path+'Ds_dominiguez_counting3.pdf',
    figsize=(5, 4),
    show_pixel=False, show_window=False,
    hide_msd=True,
)

visualisation.Ds_overlapped.go(
    'eleanorlong034',
    ['MSD_short', 'f_short'],
    export_destination=path+'Ds_dominiguez_fkt1.pdf',
    figsize=(5, 4),
    show_pixel=False, show_window=False,
    hide_msd=True,
)
visualisation.Ds_overlapped.go(
    'eleanorlong034',
    ['MSD_short', 'f_short', 'D0Sk_theory'],
    export_destination=path+'Ds_dominiguez_fkt2.pdf',
    figsize=(5, 4),
    show_pixel=False, show_window=False,
    hide_msd=True,
)
visualisation.Ds_overlapped.go(
    'eleanorlong034',
    ['MSD_short', 'f_short', 'D0Sk_theory', 'dominiguez_theory'],
    export_destination=path+'Ds_dominiguez_fkt3.pdf',
    figsize=(5, 4),
    show_pixel=False, show_window=False,
    hide_msd=True,
)

scattering_functions.show_decaytime.go(
    'eleanorlong034',
    num_Ts=1,
    export_destination=path+'f_decaytime.pdf'
)

visualisation.Ds_overlapped.go(
    'eleanorlong034',
    ['boxcounting_collective', 'timescaleint'],
    export_destination=path+'Ds_timescaleint_boxcountingcollective.pdf',
    figsize=(5, 4),
)

visualisation.Ds_overlapped.go(
    'eleanorlong034',
    ['timescaleint'],
    export_destination=path+'Ds_timescaleint.pdf',
    figsize=(5, 4),
)

visualisation.Ds_overlapped.go(
    'eleanorlong034',
    ['timescaleint_nofit_cropped', 'timescaleint'],
    export_destination=path+'Ds_timescaleint_timescaleintnofit.pdf',
    figsize=(5, 4),
)