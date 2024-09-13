from countoscope import Countoscope as CountoscopeNew
import countoscope_old
import common
import numpy as np
import matplotlib.pyplot as plt

for file in common.files_from_argv('particle_detection/data', 'particles_'):
    data = common.load(f'box_counting/data/new_vs_old_{file}.npz')

    plt.figure(figsize=(4, 3.5))
    plt.plot(data['new'], label='new')
    plt.plot(data['old'][0, :], label='old')
    plt.legend()
    plt.loglog()
    plt.xlabel('timestep')
    plt.ylabel('nmsd')
    common.save_fig(plt.gcf(), f'box_counting/figures_png/new_vs_old_{file}.png')