import czifile
import common
import numpy as np
import json

def save(file, name, pixel, depth_of_field, stack, time_step):
    print(f'{file} is {name}, px={pixel:.4f}, t={time_step:.3f}, dof={depth_of_field:.3f}')
    np.savez(f'preprocessing/data/stack_{name}', stack=stack, pixel_size=pixel, time_step=time_step,
            particle_diameter=np.nan
            #  , channels=channels)
            #, depth_of_field=depth_of_field
            )

def do_file(i, file, used_suffix):
    with czifile.CziFile(file) as czi:
        image = czi.asarray()
        stack = image.squeeze()

        metadata = czi.metadata(raw=False)["ImageDocument"]["Metadata"]
        pixel_x = metadata["Scaling"]["Items"]["Distance"][0]["Value"]
        pixel_y = metadata["Scaling"]["Items"]["Distance"][0]["Value"]
        assert pixel_x == pixel_y
        pixel_x *= 1e6

        lensNA = metadata["Information"]["Instrument"]["Objectives"]["Objective"]["LensNA"]
        depth_of_field = 0.55 / lensNA**2 # https://www.physics.hmc.edu/~gerbode/wppriv/wp-content/uploads/2012/06/DS_ZISRAW-FileFormat_Rel_ZEN2011.pdf
        # print(depth_of_field)

        time_step = metadata["Information"]["Image"]["Dimensions"]["T"]["Positions"]["Interval"]["Increment"]

        channel_datas = metadata["DisplaySetting"]["Channels"]["Channel"]
        # [print(channel) for channel in channel_datas]
        # channels = [{'color': channel['color'][:7]} for channel in channel_datas]
        # print(channels)
        print(type(channel_datas))
        if multiple_channels := type(channel_datas) is list: # more than one channel
            # channels = [{'color': channel['Color'][:7]} for channel in channel_datas]
            channels = [{'color': channel['Color'][:]} for channel in channel_datas]

        if multiple_channels:
            for channel in range(stack.shape[1]):
                channel_desc = {0:'r', 1:'g'}[channel]
                name = f'marine{used_suffix}{i}{channel_desc}'
                save(file, name, pixel_x, depth_of_field, stack[:, channel, :, :], time_step)

        else:
            name = f'marine{used_suffix}{i}'
            save(file, name, pixel_x, depth_of_field, stack, time_step)

    with open('preprocessing/data/czi.json', 'w') as o:
        json.dump(metadata, o, sort_keys=False, indent=4)
    # break
        
if __name__ == '__main__':
    for i, file in enumerate(common.get_directory_files('raw_data/marine', 'czi', 'DDTP_TL')):
        do_file(i, file, 'A')
            
    for i, file in enumerate(common.get_directory_files('raw_data/marine', 'czi', 'PPR_TL')):
        do_file(i, file, 'B')