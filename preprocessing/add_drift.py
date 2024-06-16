import common

for file in common.files_from_argv('preprocessing/data', 'stack'):
    data = common.load(f'preprocessing/data/stack_{file}.npz')
    stack = data['stack']


    # drift_width_per_stack = 0.1
    # drift_px_per_stack = drift_width_per_stack * stack.shape[1]
    # drift_px_per_frame = drift_px_per_stack / stack.shape[0]
    # print('drift px per frame', drift_px_per_frame)
    drift_px_per_frame = 1

    newstack = common.add_drift_intensity(stack, drift_px_per_frame)

    common.save_data(f'preprocessing/data/stack_{file}_drifted.npz',
                     drift=drift_px_per_frame,
                     stack=newstack,
                     pixel_size=data['pixel_size'], time_step=data['time_step'],
                     we_processed=data.get('we_processed'), nth_frame=data.get('nth_frame'))