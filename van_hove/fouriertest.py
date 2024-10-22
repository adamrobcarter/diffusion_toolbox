import common
import countoscope_theory.structure_factor
import numpy as np
import matplotlib.pyplot as plt

k = np.linspace(0, 80/2.8, 800)
phi = 0.34
sigma = 2.8
S = countoscope_theory.structure_factor.hard_spheres_2d(k, phi, sigma)
print(S)
print('nanfracs', common.nanfrac(k[1:]), common.nanfrac(S[1:]))
r, g = common.inverse_fourier(k[1:], S[1:] - 1)

print(r.shape, g.shape)
print('g nanfrac', common.nanfrac(g))
# print(r)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(4, 6))

ax1.scatter(k, S)
ax2.plot(np.abs(g))

print(g[::100])

common.save_fig(fig, f'van_hove/figures_png/fouriertest.png')