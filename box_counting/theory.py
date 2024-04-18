import sDFT_interactions
import matplotlib.pyplot as plt
import numpy as np
import common

t = np.logspace(-3, 3)
phi = 0.01
D0 = 0.04
sigma = 3
plt.plot(t, sDFT_interactions.sDFT_interactions(1, t, phi, D0, sigma), color='tab:blue')
plt.plot(t, sDFT_interactions.sDFT_interactions(2, t, phi, D0, sigma), color='tab:blue')
plt.plot(t, sDFT_interactions.sDFT_interactions(4, t, phi, D0, sigma), color='tab:blue')
phi = 0.66
plt.plot(t, sDFT_interactions.sDFT_interactions(1, t, phi, D0, sigma), color='tab:orange')
plt.plot(t, sDFT_interactions.sDFT_interactions(2, t, phi, D0, sigma), color='tab:orange')
plt.plot(t, sDFT_interactions.sDFT_interactions(4, t, phi, D0, sigma), color='tab:orange')
plt.plot(t, common.N2_nointer_2D(t, D0, 1, 1)/2, color='tab:green')
plt.plot(t, common.N2_nointer_2D(t, D0, 1, 2)/2, color='tab:green')
plt.plot(t, common.N2_nointer_2D(t, D0, 1, 4)/2, color='tab:green')
plt.loglog()
common.save_fig(plt.gcf(), 'box_counting/figures_png/theory.png')