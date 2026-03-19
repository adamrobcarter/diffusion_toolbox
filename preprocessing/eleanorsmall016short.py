import numpy as np
import common
import pandas

if __name__ == '__main__':
    df = pandas.read_pickle('raw_data/eleanor/eleanor016small/trajs_section.p')
    particles = df.to_numpy(dtype=np.float32)

    print(particles.shape)

    pixel_size = 0.288
    particles[:, [0, 1]] *= pixel_size
    # particles[:, 2] /= 40

    num_timesteps = int(particles[:, 2].max()) + 1
    window_size_x = particles[:, 1].max()
    window_size_y = particles[:, 0].max()
    print(particles[:, 1].min(), particles[:, 1].max(), particles[:, 0].min(), particles[:, 0].max())

    EDGE_CROP = 4
    particles = common.crop_particles(particles, window_size_y-EDGE_CROP, window_size_x-EDGE_CROP, EDGE_CROP, EDGE_CROP)

    metadata = dict(
        time_step = 0.25,
        # pack_frac_given = 0.342,
        particle_diameter = 2.01,
    )

    common.save_trajs(
        'eleanorsmall016short',
        particles = particles,
        particles_labels = ['y', 'x', 't', 'id'],
        **metadata
    )
    common.save_particles(
        'eleanorsmall016short',
        particles = particles[:, [0, 1, 2]],
        particles_labels = ['y', 'x', 't'],
        **metadata
    )

    end_timestep = num_timesteps // 8
    particles = particles[particles[:, 2] < end_timestep, :]

    common.save_trajs(
        'eleanorsmall016short_div8',
        particles = particles,
        particles_labels = ['y', 'x', 't', 'id'],
        **metadata
    )
    common.save_particles(
        'eleanorsmall016short_div8',
        particles = particles[:, [0, 1, 2]],
        particles_labels = ['y', 'x', 't'],
        **metadata
    )



"""
python pickle file of a pandas DataFrame

columns: 'y', 'x', 'frame', 'particle'
y and x are coordinates in pixels. to convert to um use 0.288(+/-0.003)um/px
frame is the time in frames. to convert to s use 4fps
	NB: frames are only included at 0.1fps. So the frame column goes 0,40,80...
particle is an integer ID label assigned to distinct trajectories

particles are 1.9um diameter
I get a packing fraction of 0.163 and a short-time self diffusion coefficient of 0.133um2/s

My reference for this experiment is 251218/5
    """