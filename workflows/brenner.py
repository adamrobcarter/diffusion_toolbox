import common
import matplotlib.pyplot as plt
import matplotlib.gridspec
import numpy as np
from visualisation.Ds_overlapped import get_D0
import scipy

fig = plt.figure(figsize=(5, 5))

gs = matplotlib.gridspec.GridSpec(2, 1, height_ratios=[3, 1])
ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1], sharex=ax1)

# hide x-axis labels on top plot
plt.setp(ax1.get_xticklabels(), visible=False)

z_over_a_Brenner = np.linspace(1.2, 1.6, num=100)
z_over_a = np.linspace(1, 1.6, num=500)

D_MSD_Eleanor, particle_diameter, _ = get_D0('eleanorlong001')
# D_SE = common.stokes_einstein_D(particle_diameter)
# ax1.hlines(D_SE, height_radii.min(), height_radii.max(), label='Stokes-Einstein', color=common.tab_color(2))
D_sim = 0.04095378722612345

# A = particle_diameter / (particle_diameter + 2*h) # this is radius / z_centre btw
a_over_z = 1 / z_over_a
a_over_z_Brenner = 1 / z_over_a_Brenner
Brenner_translational_correction = (1 - 9/16*a_over_z_Brenner + 1/8*a_over_z_Brenner**3 - 45/256*a_over_z_Brenner**4 - 1/16*a_over_z_Brenner**5 )
z_over_a = np.linspace(1, 1.6, num=500)
Perkin_correction = 1 - 8/15*np.log(1 - a_over_z) + 0.029*a_over_z + 0.04973*a_over_z**2 - 0.1249*a_over_z**3
Perkin_correction = 1/Perkin_correction
# from 111 years of brownian motion
print('BP', Brenner_translational_correction, Perkin_correction)
# 111 years of brownian motion says accurate for z >> a
# Faxen_translational_correction = (1 - 9/16*A + 1/8*A**3 )

data = common.load('workflows/emrfd.npz')
z = data['data'][:, 2]
print('z/a mean',z.mean()/(data['particle_diameter']/2))

particle_diameter = data['particle_diameter'] # um
particle_diameter_m = particle_diameter * 1e-6 # m

stokes_einstein_D_m = scipy.constants.k * data['T'] / (6 * np.pi * data['eta'] * particle_diameter_m/2)
stokes_einstein_D = stokes_einstein_D_m * 1e12  # um^2/s

ax2.set_xlabel('$z/a$')
ax1.set_ylabel('$D$')
ax1.set_ylim(0, stokes_einstein_D*1.05)

ax1.hlines(stokes_einstein_D, z_over_a.min(), z_over_a.max(), label='Stokes-Einstein', color=common.tab_color(2), linestyle='dashed')
ax1.plot(z_over_a_Brenner, stokes_einstein_D * Brenner_translational_correction, label='Brenner/Faxen', color=common.tab_color(2))
ax1.plot(z_over_a, stokes_einstein_D * Perkin_correction, label='Perkin & Jones', color=common.tab_color(2), linestyle='dotted')
ax1.legend(loc='lower right')
common.save_fig(fig, 'workflows/figures/brenner0.png')

ax1.hlines(D_MSD_Eleanor, z_over_a.min(), z_over_a.max(), label='experiment', color=common.tab_color(1))
ax1.legend(loc='lower right')
common.save_fig(fig, 'workflows/figures/brenner1.png')

ax1.hlines(D_sim, z_over_a.min(), z_over_a.max(), label='simulation', color=common.tab_color(0))
ax2.hist(z/(data['particle_diameter']/2), bins=z_over_a, color=common.tab_color(0))
ax1.legend(loc='lower right')
common.save_fig(fig, 'workflows/figures/brenner2.png')

ax1.vlines(z.mean()/(data['particle_diameter']/2), *ax1.get_ylim(), color=common.tab_color(0), alpha=0.5)
ax2.vlines(z.mean()/(data['particle_diameter']/2), *ax2.get_ylim(), color=common.tab_color(0), alpha=0.5)
ax1.legend(loc='lower right')
common.save_fig(fig, 'workflows/figures/brenner3.png')
