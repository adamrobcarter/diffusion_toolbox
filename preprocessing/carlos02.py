import common
import tifffile
import numpy as np

file = 'raw_data/carlos/PS_Tiff_Stack_maybe2fps.tif'
stack = common.load_tif(file)

# with tifffile.TiffFile(file) as tif:
#     stack = tif.asarray()
#     # axes = tif.series[0].axes 'ZYX'
#     imagej_metadata = tif.imagej_metadata
#     print(tif.metaseries_metadata)
#     print(tif.eer_metadata)
#     # print(tif.pilatus_metadata)
#     # for key in tif:
#     #     print(key)
#     # print(dir(tif))
#     for prop in dir(tif):
#         print(prop)
#         print(getattr(tif, prop))
#         print()

# print(stack.shape)

time_step = 0.5
max_time_hours = round(stack.shape[0]*time_step/60/60, 2)
pixel_size = 0.214

common.save_data('preprocessing/data/stack_carlos02.npz',
    stack=stack.astype(np.float32),
    time_step=time_step, pixel_size=pixel_size, particle_diameter=2,
    window_size_x=stack.shape[1]*pixel_size, window_size_y=stack.shape[2]*pixel_size,
    max_time_hours=max_time_hours,
    source_file=file,
)