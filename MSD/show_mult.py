import matplotlib.pyplot as plt
import common
import numpy as np
import matplotlib.cm

if __name__ == '__main__':
    fig, ax = plt.subplots(1, 1, figsize=(5.5, 5))

    files = common.files_from_argv('MSD/data', 'msd_')

    for i, file in enumerate(files):
        data = common.load(f'MSD/data/msd_{file}.npz')
        msd = data['msd']
        msd_unc = data['msd_unc']
        
        t_indexes = np.arange(0, msd.size)
        n = (msd.size-1) / t_indexes
        msd_unc = msd_unc / np.sqrt(n)
        # ^^ this is based on thinking about how many independent measurements there are


        t = data['t']
        print(t)

        color =  matplotlib.cm.afmhot(np.interp(i, (0, len(files)), (0.2, 0.75)))
                
        # ax.errorbar(t[1:], msd[1:], msd_unc[1:], linestyle='none', marker='none', color='lightskyblue')
        ax.plot(t[1:], msd[1:], marker='.', markersize=8, linestyle='none', color=color, label=file)
        ax.fill_between(t[1:], msd[1:]-msd_unc[1:], msd[1:]+msd_unc[1:], alpha=0.2, color=color)
        
    ax.loglog()
    ax.set_ylim(msd[1:].min()*0.6, msd.max()/0.8)
    # ax.set_xlim(t[1]*0.8, t[-1]/0.8)

    ax.set_ylabel(r'$\langle r(t)^2 \rangle$ ($\mathrm{\mu m}$)')
    ax.set_xlabel('$t$ (s)')

    # common.save_fig(fig, f'/home/acarter/presentations/cin_first/figures/msd_nofit_{file}.pdf', hide_metadata=True)

    ax.legend(fontsize=8)

    filename = '_'.join(files)
    # common.save_fig(fig, f'/home/acarter/presentations/cin_first/figures/msd_{file}.pdf', hide_metadata=True)
    common.save_fig(fig, f'MSD/figures_png/msd_combined_{filename}.png')
