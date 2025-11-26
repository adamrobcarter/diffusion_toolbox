from preprocessing.libmobility_diffusion import rsync_and_load, get_name
import numpy as np
import common

keys_that_dont_need_to_match = [
    'job_started_at',
    'source_file',
    'extra_source_file'
]

def load_and_save(source_dirs, nums, skip_rsync=False):
    dirs = [source_dirs.format(num) for num in nums]

    datas, metadatas = [], []

    for dir in dirs:
        data, metadata = rsync_and_load(dir, skip_rsync=skip_rsync)
        datas.append(data)
        metadatas.append(metadata)

    for key in metadatas[0]:
        if key in keys_that_dont_need_to_match:
            continue

        for i in range(1, len(metadatas)):
            assert metadatas[i][key] == metadatas[0][key], f'key={key}, {metadatas[i][key]} != {metadatas[0][key]}'

    new_metadata = dict(metadatas[0])
    del new_metadata['job_started_at']
    new_metadata['jobs_started_at'] = [metadata['job_started_at'] for metadata in metadatas]
    del new_metadata['source_file']
    new_metadata['source_files'] = [metadata['source_file'] for metadata in metadatas]
    del new_metadata['extra_source_file']
    new_metadata['extra_source_files'] = [metadata['extra_source_file'] for metadata in metadatas]

    new_data = np.stack(datas, axis=0)
    
    output_name = get_name(metadatas[0], prefix=f'array{len(nums)}_')
    
    common.save_data(f'particle_detection/data/particles_{output_name}.npz',
        particles_array=new_data, **new_metadata
    )

# load_and_save(
#     '/store/cartera/libmobility_diffusion/solver_NBody_N_7638_L_640_single_wall_dt_125_t_3600_32_run_{}/',
#     range(8),
#     skip_rsync=True
# )
# load_and_save(
#     '/store/cartera/libmobility_diffusion/solver_NBody_N_7638_L_640_single_wall_dt_125_t_450_32_run_{}/',
#     range(64),
#     skip_rsync=True
# )
# load_and_save(
#     '/store/cartera/libmobility_diffusion/solver_NBody_N_30552_L_1280_single_wall_dt_125_t_3600_32_run_{}/',
#     range(8),
#     skip_rsync=True
# )
# load_and_save(
#     '/store/cartera/libmobility_diffusion/solver_NBody_N_7638_L_640_single_wall_dt_125_t_1025_32_run_{}/',
#     range(28),
# )
load_and_save(
    '/store/cartera/libmobility_diffusion/solver_NBody_N_30552_L_1280_single_wall_dt_125_t_1025_32_run_{}/',
    range(27),
)