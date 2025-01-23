import common
import visualisation.Ds_overlapped_mult
import matplotlib.pyplot as plt
import cmocean

SINGLE_FIGSIZE = (4.5, 4)
path = '/home/acarter/presentations/countoscope_24_december/figures'

########################################################
fig, ax = plt.subplots(1, 1, figsize=SINGLE_FIGSIZE)

kwargs = dict(
    files = ['sim_nohydro_011_L160_longer_mergedD', 'sim_nohydro_011_L320_longer_mergedD', 'sim_nohydro_011_L640_longer_mergedD'],
    file_labels = ['$L_x=160\mathrm{\mu m}$', '$L_x=320\mathrm{\mu m}$', '$L_x=640\mathrm{\mu m}$'],
    colors = [[cmocean.cm.ice(0.75)], [cmocean.cm.ice(0.5)], [cmocean.cm.ice(0.25)]]
)

visualisation.Ds_overlapped_mult.go(
    sources = ['timescaleint_nofit_cropped_var'],
    ax      = ax,
    markers = 'o',
    **kwargs
)
ax.set_ylim(0.87, 2.5)

common.save_fig(fig, f'{path}/periodic_counting.pdf', hide_metadata=True)

########################################################
fig, ax = plt.subplots(1, 1, figsize=SINGLE_FIGSIZE)

visualisation.Ds_overlapped_mult.go(
    sources = ['f_first_first'],
    ax      = ax,
    markers = 'd',
    disable_ylabel=True,
    **kwargs
)
ax.set_ylim(0.87, 2.5)

common.save_fig(fig, f'{path}/periodic_fkt.pdf', hide_metadata=True)


########################################################
fig, ax = plt.subplots(1, 1, figsize=SINGLE_FIGSIZE)

visualisation.Ds_overlapped_mult.go(
        ['sim_nohydro_002_L320', 'sim_nohydro_002_L640_crop320'],
        sources=[
            'f_first_first'
        ],
        ax=ax,
        plot_against_k=True,
        markers='d',
        file_labels=['periodic boundary', 'non periodic boundary'],
        source_labels=[''],
        colors=[['olivedrab'], ['darkgreen']]
    )
ax.set_ylim(0, 30)
common.save_fig(fig, f'{path}/periodic_non_fkt.pdf', hide_metadata=True)


########################################################
fig, ax = plt.subplots(1, 1, figsize=SINGLE_FIGSIZE)

visualisation.Ds_overlapped_mult.go(
        ['eleanorlong001'],
        sources=[
            'f_first_first'
        ],
        ax=ax,
        plot_against_k=True,
        markers='d',
        file_labels=['experiment'],
        source_labels=[''],
        colors=[['olivedrab']]
    )
ax.set_ylim(0, 30)
common.save_fig(fig, f'{path}/periodic_fkt_exp.pdf', hide_metadata=True)


########################################################
fig, ax = plt.subplots(1, 1, figsize=SINGLE_FIGSIZE)

visualisation.Ds_overlapped_mult.go(
        ['sim_nohydro_002_L320'],
        sources=[
            'f_first_first'
        ],
        ax=ax,
        plot_against_k=True,
        markers='d',
        file_labels=['simulation'],
        source_labels=[''],
        colors=[['olivedrab']]
    )
ax.set_ylim(0, 30)
common.save_fig(fig, f'{path}/periodic_fkt_single.pdf', hide_metadata=True)


########################################################
fig, ax = plt.subplots(1, 1, figsize=SINGLE_FIGSIZE)

visualisation.Ds_overlapped_mult.go(
        ['sim_nohydro_002_L320', 'sim_nohydro_002_L640_crop320'],
        sources=[
            'timescaleint_nofit_cropped_var'
        ],
        ax=ax,
        markers='o',
        file_labels=['periodic boundary', 'non periodic boundary'],
        source_labels=[''],
        colors=[['olivedrab'], ['darkgreen']]
    )
ax.set_ylim(0, 30)
common.save_fig(fig, f'{path}/periodic_non_counting.pdf', hide_metadata=True)