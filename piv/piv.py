import openpiv.pyprocess, openpiv.validation
import numpy as np
import tqdm

def go(stack, time_step, interrogation_window):
    overlap = [interrogation_window/2, interrogation_window/2]

    u = []
    v = []
    signal_to_noise = []

    for frame in tqdm.trange(1, stack.shape[0]):
        us, vs, signal_to_noises = openpiv.pyprocess.extended_search_area_piv(
            stack[frame-1, :, :], stack[frame, :, :],
            window_size=interrogation_window, overlap=overlap, dt=time_step
        )

        u.append(us)
        v.append(vs)
        signal_to_noise.append(signal_to_noises)

    u = np.array(u)
    v = np.array(v)
    signal_to_noise = np.array(signal_to_noise)

    return u, v, signal_to_noise

