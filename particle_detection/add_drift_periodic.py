import common
import tqdm
import functools
import numpy as np

# V_PROFILE = 'cos'
# V_PROFILE = 'cos2'
# V_PROFILE = 'step'
V_PROFILE = 'const'
# V_PROFILE = 'step2'

def get_velocity_cos(window_size_x, window_size_y, x, y):
    v = np.full((*x.shape, 2), np.nan)
    print(x.shape)
    
    v[..., 0] = np.cos(2*np.pi * y / window_size_y)
    v[..., 1] = 0
    # print(v.shape)
    # print(v)
    assert np.isfinite(v).all()
    return v

def get_velocity_cos2(window_size_x, window_size_y, x, y):
    v = np.full((*x.shape, 2), np.nan)
    print(x.shape)
    
    v[..., 0] = 1 + 0.5 * np.cos(2*np.pi * y / window_size_y)
    v[..., 1] = 0
    # print(v.shape)
    # print(v)
    assert np.isfinite(v).all()
    return v

def get_velocity_step(window_size_x, window_size_y, x, y):
    v = np.full((*x.shape, 2), np.nan)
    
    v[y <=   window_size_y,   :] = np.array([-1/2, 0])
    v[y <  2*window_size_y/3, :] = np.array([1, 0])
    v[y <    window_size_y/3, :] = np.array([0, 0])
    # if y < window_size_y/3:
    #     return
    # elif y < 2*window_size_y/3:
    #     return np.array([1, 0])
    # else:
    #     return np.array([-0.5, 0])
    assert np.isfinite(v).all()
    return v

def get_velocity_const(window_size_x, window_size_y, x, y):
    v = np.full((*x.shape, 2), 1)
    v[..., 1] = 0
    
    assert np.isfinite(v).all()
    return v

def get_velocity_step2(window_size_x, window_size_y, x, y):
    v = np.full((*x.shape, 2), np.nan)
    
    v[y <= 10*window_size_y/10, :] = np.array([-1.0, 0])
    v[y <   9*window_size_y/10, :] = np.array([-0.8, 0])
    v[y <   8*window_size_y/10, :] = np.array([-0.6, 0])
    v[y <   7*window_size_y/10, :] = np.array([-0.4, 0])
    v[y <   6*window_size_y/10, :] = np.array([-0.2, 0])
    v[y <   5*window_size_y/10, :] = np.array([ 0.0, 0])
    v[y <   4*window_size_y/10, :] = np.array([ 0.2, 0])
    v[y <   3*window_size_y/10, :] = np.array([ 0.4, 0])
    v[y <   2*window_size_y/10, :] = np.array([ 0.6, 0])
    v[y <   1*window_size_y/10, :] = np.array([ 0.8, 0])
    # if y < window_size_y/3:
    #     return
    # elif y < 2*window_size_y/3:
    #     return np.array([1, 0])
    # else:
    #     return np.array([-0.5, 0])
    assert np.isfinite(v).all()
    return v

get_v = {
    'const': get_velocity_const,
    'cos': get_velocity_cos,
    'cos2': get_velocity_cos2,
    'step': get_velocity_step,
    'step2': get_velocity_step2
}

def go_internal(infile, outfile, v_profile, velocity_multiplier=1):
    data = common.load(infile)

    particles = data['particles']
    window_size_x = data['window_size_x']
    window_size_y = data['window_size_y']
    time_step     = data['time_step']

    assert np.all(particles[:, 0] <= window_size_x), 'needed for get_velocity'
    assert np.all(particles[:, 1] <= window_size_y), 'needed for get_velocity'

    v = functools.partial(get_v[v_profile], window_size_x, window_size_y)

        # add velocity
    x = particles[:, 0]
    y = particles[:, 1]
    frame = particles[:, 2]
    assert np.unique(frame)[1] == 1 # checks that this is frame number not timestep

    new_xy = v(x, y) * velocity_multiplier * frame[:, np.newaxis] * time_step
    particles[:, 0] += new_xy[:, 0]
    particles[:, 1] += new_xy[:, 1]

        # rewrap into periodic domain
    particles[:, 0] = particles[:, 0] % window_size_x
    particles[:, 1] = particles[:, 1] % window_size_y

    assert np.isfinite(particles).all()

    newdata = dict(data)
    newdata['particles'] = particles
    newdata['v_profile'] = v_profile
    newdata['velocity_multiplier'] = velocity_multiplier
    common.save_data(outfile, **newdata)

def go(file, v_profile, velocity_multiplier=1):
    go_internal(
        infile = f'particle_detection/data/particles_{file}.npz',
        outfile = f'particle_detection/data/particles_{file}_drifted_{v_profile}_v{velocity_multiplier}.npz',
        v_profile = v_profile,
        velocity_multiplier = velocity_multiplier,
    )

if __name__ == '__main__':
    for file in common.files_from_argv('particle_detection/data', 'particles_'):
        go(
            file = file,
            v_profile = V_PROFILE,
        )