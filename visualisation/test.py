import common
import countoscope_theory.timescaleint
import matplotlib.pyplot as plt
import numpy as np

L = np.logspace(-2, 3)
D = countoscope_theory.timescaleint.D_of_L(L, 0.01, 0.34, 2.8)
print(D)
plt.loglog(L, -np.array(D))
common.save_fig(plt.gcf(), 'visualisation/figures_png/test.png')