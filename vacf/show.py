import common
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate

LOG = False

for file in common.files_from_argv('vacf/data', 'vacf_'):
    data = common.load(f'vacf/data/vacf_{file}.npz')
    vacf = data['vacf']
    time_step = data['time_step']

    print(vacf[0:10])

    fig, ax = plt.subplots()

    ax.plot(time_step * np.arange(vacf.size), vacf/vacf[0], marker='o')

    if 'particle_mass' in data:
        t = np.linspace(0, time_step*100, num=1000)
        m = data['particle_mass'] / 1e12 # comes in pN
        friction = 6*np.pi * data['eta'] * (data['particle_diameter']/1e6) / 2
        # eta is kg s^-1 m–1, so this is kg / s
        ax.plot(t, np.exp(-friction * t / m), label='theory')

    ax.set_xlabel('$t$')
    ax.set_ylabel('$C(t)/C(0)$')
    if LOG:
        ax.semilogy()
        ax.semilogx()
        ax.set_ylim(1e-4, 1.5)
    else:
        ax.set_xlim(0, 6)
    ax.legend()
    ax.hlines(0, *ax.get_xlim(), colors='gray', alpha=0.5)

    common.save_fig(fig, f'vacf/figures/vacf_{file}.png')

    D = scipy.integrate.trapezoid(vacf[:100], dx=time_step)
    print(f'D = {D:.3g} um^2/s')
    