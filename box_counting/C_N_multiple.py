import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize
import common
import scipy.integrate, scipy.special, scipy.optimize, scipy.signal
import sDFT_interactions
import tqdm
import sys

"""
Countoscope appendix:

Given the above estimates for h N2(t)i and hN2i hNi2, the timescale integral in equation (8) 
is calculated using trapezoidal summation. When the observation box size is very small, the 
integrand in Eq. (8) decays very quickly, and much of its support may fall before the first 
frame. When the observation box is very large, the decay is slow and we may not have enough 
data to calculate the integral accurately. Both of these issues can be addressed by fitting 
our available data at short and long times to theory-informed functional forms and calculating 
the missing contributions to the timescale integral. This procedure is described in detail in 
Supp. Mat. Sec. 1.5 which also includes a discussion of how the errorbars were calculated for 
Fig. 5.
"""
f = lambda tau: np.sqrt(tau / np.pi) * ( np.exp(-1/tau) - 1) + scipy.special.erf(np.sqrt(1/tau)) # countoscope eq. 2

# import warnings
# warnings.filterwarnings('ignore')

integrand_fig, (integrand_ax_old, mid_ax, extra_ax) = plt.subplots(3, 1, figsize=(6, 15))
integ_fig, integ_axs = plt.subplots(1, 1, figsize=(6, 5))

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
color_index = 0

for file in sys.argv[1:]:

    

    # ax.set_prop_cycle(None) # reset colour cycle
    color = colors[color_index]
    color_index += 1

    # D0 = { # countoscope paper, table 1
    #     0.02: 0.0416,
    #     0.34: 0.0310,
    #     0.66: 0.0175
    # }[phi]

    # data = common.load(f'box_counting/data/counted_dense_{file}.npz')
    data = common.load(f'box_counting/data/counted_{file}.npz')
    N2_mean = data['N2_mean']
    N2_std  = data['N2_std']
    N_stats = data['N_stats']
    sigma = data['particle_diameter']
    phi = data['pack_frac']
    time_step    = data['time_step']

    box_sizes = N_stats[:, 0]
    sep_sizes = data['sep_sizes']
    N_mean    = N_stats[:, 1]
    N_var     = N_stats[:, 2]
    num_boxes_used = N_stats[:, 5]

    num_timesteps = N2_mean.shape[1]
    num_boxes     = N2_mean.shape[0]
    t = np.arange(0, num_timesteps) * time_step

    D_of_Ls = np.full((num_boxes), np.nan)
    T_of_Ls = np.full((num_boxes), np.nan)

    # sDFT_data = common.load(f'data/sDFT_{phi}.npz')
    # late_integrals    = sDFT_data['late_integrals']
    # integrand_theorys = sDFT_data['integrand_theorys']
    # assert np.all(sDFT_data['box_sizes'] == box_sizes)

    for box_size_index in tqdm.trange(num_boxes):
        # T_of_Ls[box_size_index] = common.T_of_L(N2_mean[box_size_index, :], N_var[box_size_index], t)
        L = box_sizes[box_size_index]
        # D_of_Ls[box_size_index] = common.D_of_L(N2_mean[box_size_index, :], N_var[box_size_index], t, L) / D0

        T_integrand = (1 - 0.5 * N2_mean[box_size_index, :] / N_var[box_size_index] )**2 # countoscope paper eq. 8
        C_N = N_var[box_size_index] - 0.5 * N2_mean[box_size_index, :] # countoscope paper eq. 1
        interval_end = (C_N[1] / N_var[box_size_index])**2
        interval_end_2 = T_integrand[1]
        # print(t[1], interval_end, interval_end_2)
        

        # scipy.integrate.quad(early_func, 0, t[1])
        # T_of_Ls[box_size_index] = 2 * (scipy.integrate.quad(fit_func, 0, t[1])[0] + scipy.integrate.trapezoid(T_int[1:], t[1:])) # countoscope eq. 8
        
    
        L_2 = L

        # popt_full, pcov_full = scipy.optimize.curve_fit(integrand_full, t[:M2_index], T_integrand[:M2_index], maxfev=2000)

    
        if box_size_index%1 == 0:

            integrand_ax = integ_axs

            # integrand_ax.set_title(rf'$\phi={phi}$')

            # N_var_theory = N_mean[box_size_index] * sDFT_interactions.sDFT_interactions(L=L, t=np.array([np.inf]), phi=phi, D0=D0, sigma=sigma)[0]
            # print('var', N_var_theory, N_var[box_size_index], N_var_theory / N_var[box_size_index])
            #              ^^^^ should be density * box volume?

            # print(f'early data late {early_integral:.1f}, {data_integral:.1f}, {late_integral:.1f}')
            # print('ratio', late_integral/late_integral1)
            N2_func_full = lambda t, D0: 2 * N_mean[box_size_index] * (1 - f(4*D0*t/L**2) * f(4*D0*t/L**2)) # countoscope eq. 2, countoscope overleaf doc
            t_C_N_theory = np.logspace(-2, 5, 200)
            # N2_theory = N2_func_full(t_C_N_theory, D0)
            # integrand_theory = (1 - 0.5 * N2_theory / N_var_theory)**2
            # integrand_theory2 = (1 - 0.5 * N2_theory / N_mean[box_size_index])**2
            # integrand_ax.plot(t_C_N_theory, integrand_theory, color='grey', linewidth=0.8, zorder=5)
            # integrand_ax.plot(t_C_N_theory, integrand_theory2, color='black', linewidth=0.5, zorder=5)
            # integrand_rescaled_ax.plot(t_C_N_theory/L**2, integrand_theory, color='black', linewidth=0.5, zorder=5, alpha=1)

            label = rf'$L = {L}\mathrm{{\mu m}}$'
            # label += f', D={D0:.3f}'
            label += fr', $sep = {sep_sizes[box_size_index]:.1f}\mathrm{{\mu m}}$'
            label += f', $n = {num_boxes_used[box_size_index]:.0f}$'
            line = integrand_ax.plot(t[1:], T_integrand[1:], label=label, color=color, zorder=5)
            
            
integrand_ax.legend(fontsize=8)
integrand_ax.loglog()
integrand_ax.set_xlabel('$t$')

integrand_ax.set_ylim(1e-5, 1e0)
# integrand_ax.set_xlim(0.1, 1e5) 

integrand_ax.set_xlabel('$t$')
integrand_ax.set_ylabel('$C_N(t)^2 / C_N(0)^2$')

# integrand_ax.seintegrand_xlabel('$L/\sigma$')
# integrand_ax.seintegrand_ylabel(r'$integrand(L)$')
# integrand_ax.set_ylim(1e-4, 1.1e0)
# integrand_fig.savefig('figures_png/integrand.png', dpi=300)
integrand_ax.set_title(fr'{file} timescale integrand, $\phi={phi:.3f}$, $\sigma={sigma}$')
integ_fig.savefig(f'box_counting/figures_png/C_N_multiple_{file}.png', dpi=300)
integ_fig.savefig(f'box_counting/figures_eps/C_N_multiple_{file}.eps', dpi=300)