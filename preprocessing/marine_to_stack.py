import czifile
import common
import numpy as np
import json


for i, file in enumerate(common.get_directory_files('raw_data/marine', 'czi')):
    # image = czifile.imread(file)
    
    with czifile.CziFile(file) as czi:
        image = czi.asarray()
        stack = image.squeeze()

        metadata = czi.metadata(raw=False)["ImageDocument"]["Metadata"]
        pixel_x = metadata["Scaling"]["Items"]["Distance"][0]["Value"]
        pixel_y = metadata["Scaling"]["Items"]["Distance"][0]["Value"]
        assert pixel_x == pixel_y
        pixel_x *= 1e6

        time_step = metadata["Information"]["Image"]["Dimensions"]["T"]["Positions"]["Interval"]["Increment"]

        print(f'{file} is marine{i}, px={pixel_x:.4f}, t={time_step:.3f}')
        np.savez(f'preprocessing/data/stack_marine{i}', stack=stack, pixel_size=pixel_x, time_step=time_step, particle_diameter=np.nan)
    

    # with open('czi.json', 'w') as o:
    #     json.dump(metadata, o, sort_keys=False, indent=4)
    # break