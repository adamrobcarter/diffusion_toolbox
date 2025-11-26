import visualisation.Ds_overlapped
import MSD.calc
import particle_detection.show_movie

nohydro = 'ld_self_stripe_0.07_nowalls_ztrap_w0p5a_h0a_L500_t1h_12s_dt50_nolub'
hydro = 'ld_dpstokes_stripe_0.07_nowalls_ztrap_w0p5a_h0a_L500_t1h_20s_dt100_nolub'

# MSD.calc.go(f'{hydro}_unwrap')
# MSD.calc.go(f'{nohydro}_unwrap')

D0_hydro, _, _ = visualisation.Ds_overlapped.get_D0(hydro)
D0_nohydro, _, _ = visualisation.Ds_overlapped.get_D0(nohydro)

print(D0_hydro, D0_nohydro)


particle_detection.show_movie.go(
    hydro,
    outfile = f'workflows/figures/movie_hydro.gif',
    every_nth_frame=10,
    max_num_frames = 15
)
particle_detection.show_movie.go(
    nohydro,
    outfile = f'workflows/figures/movie_nohydro.gif',
    every_nth_frame=10,
    max_num_frames = 15
)