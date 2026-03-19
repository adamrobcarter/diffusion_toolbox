import common
import visualisation.Ds_overlapped_mult
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 1)
visualisation.Ds_overlapped_mult.go(
    [
        dict(
            file = 'sim_nointer_0.1_L207_t1m_dt100ms_sigma0.5_precropL828',
            source = 'f_first_first'
        ),
        dict(
            file = 'sim_nointer_0.1_L207_t24h_dt4m_sigma0.5_precropL828',
            source = 'f_t1024',
        )
    ],
    ax = ax,
)
common.save_fig(fig, 'workflows/figures/el.png')