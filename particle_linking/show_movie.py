import common
import particle_detection.show_movie

for file in common.files_from_argv('preprocessing/data', 'stack_'):
    particle_detection.show_movie.go(
        file,
        stack_file=f'preprocessing/data/stack_{file}.npz',
        particles_file=f'particle_linking/data/trajs_{file}.npz',
        output_file=f"particle_linking/figures_png/movie_{file}.gif")