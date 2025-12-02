import scipy.optimize
import common
import matplotlib.pyplot as plt
import scipy
import pickle
import numpy as np

if __name__ == '__main__':
    for file in common.files_from_argv('particle_detection/data', 'particles_'):
        
        with open(f'isf/data/F_quantify_noise_{file}.pickle', "rb") as f:
            all_data = pickle.load(f)

        num_time_origins = []
        error = []
        comptimes = []

        for data in all_data:
            
            t         = data["t"]
            F_all     = data["F"]
            F_unc_all = data['F_unc']
            k_all     = data["k"]

            f_all = F_all / F_all[0, :]

            errors_inside = []

            for k_index in range(f_all.shape[1]):
                f = f_all[:, k_index]

                if np.isnan(f).all():
                    continue

                func = lambda x, m, c: m*x + c
                (m, c), pcov = scipy.optimize.curve_fit(func, t[1:], f[1:])
                assert 0.8 < c < 1.2, f'c={c:.3f}'

                pred_vals = func(t, m, c)
                # err = np.sum((pred_vals[1:] - f[1:])**2)
                err = np.sqrt(pcov[0, 0]) / np.abs(m)
                print(err)
                errors_inside.append(err)

            num_time_origins.append(data['max_time_origins'])
            error.append(np.mean(errors_inside))
            comptimes.append(data['computation_time'])

        fig, ax = plt.subplots(1, 1)

        ax.scatter(num_time_origins, error)
        ax.semilogx()
        ax.semilogy()
        ax.set_xlabel('max time origins')
        ax.set_ylabel('error on gradient')
        
        ax.yaxis.label.set_color('tab:blue')
        # par1.spines["right"].set_edgecolor(p2.get_color())
        ax.tick_params(axis='y', colors='tab:blue')

        ax2 = ax.twinx()
        ax2.scatter(num_time_origins, comptimes, color='tab:orange', marker='.') # can confirm comptime is linear in num_time_origins
        ax2.semilogy()
        ax2.set_ylabel('computation time (s)')
        
        ax2.yaxis.label.set_color('tab:orange')
        # par1.spines["right"].set_edgecolor(p2.get_color())
        ax2.tick_params(axis='y', colors='tab:orange')


        common.save_fig(plt.gcf(), f'isf/data/f_quantify_noise_{file}.png')

        