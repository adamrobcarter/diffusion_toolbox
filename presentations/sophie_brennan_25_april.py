import matplotlib.pyplot as plt
import common
import visualisation.Ds_overlapped_mult

PATH = '/home/acarter/presentations/sophie_brennan_25_april/figures'

fig, (ax_boxcounting, ax_fkt) = plt.subplots(1, 2, figsize=(6.5, 3.5))
ylim = (0.7, 7)
visualisation.Ds_overlapped_mult.go(
    [
        ('sim_hydro_011_L640_pot',       'timescaleint_nofit_cropped_var'),
    ],
    ax_boxcounting,
    labels = ['single wall']
)
visualisation.Ds_overlapped_mult.go(
    [
        ('sim_hydro_011_L640_pot',       'f_t128'),
    ],
    ax_fkt,
    labels = ['single wall']
)
ax_boxcounting.set_title('Counting')
ax_fkt        .set_title('$f(k, t=128\mathrm{s})$')
ax_boxcounting.set_ylim(*ylim)
ax_fkt        .set_ylim(*ylim)
common.save_fig(fig, f'{PATH}/zconf.png')

fig, (ax_boxcounting, ax_fkt) = plt.subplots(1, 2, figsize=(6.5, 3.5))
ylim = (0.7, 7)
visualisation.Ds_overlapped_mult.go(
    [
        ('sim_hydro_011_L640_pot',       'timescaleint_nofit_cropped_var'),
        ('sim_hydro_011_L640_pot_zconf', 'timescaleint_nofit_cropped_var'),
    ],
    ax_boxcounting,
    labels = ['single wall', 'z trap']
)
visualisation.Ds_overlapped_mult.go(
    [
        ('sim_hydro_011_L640_pot',       'f_t128'),
        ('sim_hydro_011_L640_pot_zconf', 'f_t128'),
    ],
    ax_fkt,
    labels = ['single wall', 'z trap']
)
ax_boxcounting.set_title('Counting')
ax_fkt        .set_title('$f(k, t=128\mathrm{s})$')
ax_boxcounting.set_ylim(*ylim)
ax_fkt        .set_ylim(*ylim)
common.save_fig(fig, f'{PATH}/zconf2.png')