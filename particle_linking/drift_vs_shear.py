import common
import matplotlib.pyplot as plt
import numpy as np
import scipy

fig, ax = plt.subplots(1, 1)
ax.set_xlabel('shear $S$ ($v_z = Sz$)')
ax.set_ylabel('drift x ($\mathrm{\mu m/s}$)')

shears = []
drifts = []
zs = []

for file in common.files_from_argv('particle_linking/data', 'trajs_'):
    data = common.load(f'particle_linking/data/trajs_{file}.npz')
    particle_diameter = data['particle_diameter']   
    drift_xyz, drift_xys_unc = common.find_drift(data['particles'], data.get('dimension', 2))
    z_avg = data['particles'][:, 2].mean()

    # print('drift', drift_xyz)
    shears.append(data['shear'])
    drifts.append(drift_xyz[0]) # um/s
    zs.append(z_avg)

zs = np.array(zs)
shears = np.array(shears)

target_v = 0.10356

func = lambda x, m, c: m*x + c
popt, pcov = common.curve_fit(func, drifts, shears)
print(f'S = {func(target_v, *popt)}')


ax.scatter(shears, drifts, label='simulation nblobs=642')
ax.legend()
common.save_fig(fig, f'particle_linking/figures_png/drift_vs_shear_{file}_0.png')

# func = lambda x, m, c: m*x+c
# popt, unc = common.curve_fit(func, shears, drifts)
# ax.plot(shears, func(shears, *popt), label=f'fit y = {popt[0]:.2f}x = {popt[1]:.2f}', color='black')

shear_theory = np.linspace(min(shears), max(shears))
shear_theory = shears
S = shear_theory
a = particle_diameter / 2
h = zs
delta = h - a
print('d/a', delta/a)

U = h * S * 0.7431 / (0.6376 - 0.200 * np.log(delta/a)) # https://www.sciencedirect.com/science/article/pii/0009250967800484 eq 4.11
ax.plot(shear_theory, U, label='Goldman 1976b')
print(U)

ax.legend()
common.save_fig(fig, f'particle_linking/figures_png/drift_vs_shear_{file}_1.png')